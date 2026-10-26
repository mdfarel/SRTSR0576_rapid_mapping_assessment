"""Static layers per site: Copernicus DEM, ESA WorldCover 2021 and JRC Global Surface Water.

Each layer is mosaicked from the Planetary Computer tiles that touch the site and
resampled onto the 10 m site grid. Writes data/<site>/static/dem.tif,
worldcover.tif and jrc_<layer>.tif; a layer already on disk is skipped.

    python -m pipeline.fetch_static [site ...]
"""
import sys

import numpy as np
from rasterio.enums import Resampling

import config
from pipeline.common import bbox_wgs84, catalog, read_to_grid, site_dir, site_grid, write_tif

JRC_LAYERS = ["occurrence", "seasonality", "recurrence", "extent", "transitions"]
BBOX_PAD_DEG = 0.01  # search a little beyond the grid so edge tiles are not missed
DEM_NODATA = -32767.0
WORLDCOVER_YEAR = "2021"


def mosaic(items, asset, grid, resampling, dtype, nodata):
    """Read one asset of every tile onto the grid; where tiles overlap, the first tile read wins."""
    out = None
    for it in items:
        if asset not in it.assets:
            continue
        arr, _ = read_to_grid(it.assets[asset].href, grid, resampling, dtype, nodata)
        if out is None:
            out = np.full(arr.shape, nodata, dtype)
        fill = (arr != nodata) & (out == nodata)
        out[fill] = arr[fill]
    return out


def is_worldcover_year(item):
    """True for a WorldCover tile of WORLDCOVER_YEAR (the collection also holds the 2020 map)."""
    date = item.datetime or item.properties.get("start_datetime", "")
    return str(date).startswith(WORLDCOVER_YEAR)


def fetch(site_key):
    """Download whichever static layers are not yet on disk for one site."""
    cat = catalog()
    g10 = site_grid(site_key, 10)
    bbox = bbox_wgs84(g10, BBOX_PAD_DEG)
    d = site_dir(site_key, "static")

    if not (d / "dem.tif").exists():
        items = list(cat.search(collections=["cop-dem-glo-30"], bbox=bbox).items())
        dem = mosaic(items, "data", g10, Resampling.bilinear, np.float32, DEM_NODATA)
        write_tif(d / "dem.tif", dem, g10, nodata=DEM_NODATA)
        print(site_key, "dem", len(items), "tiles")

    if not (d / "worldcover.tif").exists():
        all_items = cat.search(collections=["esa-worldcover"], bbox=bbox).items()
        items = [it for it in all_items if is_worldcover_year(it)]
        wc = mosaic(items, "map", g10, Resampling.nearest, np.uint8, 0)
        write_tif(d / "worldcover.tif", wc, g10, nodata=0)
        print(site_key, "worldcover 2021", len(items), "tiles")

    items = list(cat.search(collections=["jrc-gsw"], bbox=bbox).items())
    for layer in JRC_LAYERS:
        if (d / f"jrc_{layer}.tif").exists():
            continue
        arr = mosaic(items, layer, g10, Resampling.nearest, np.uint8, 255)
        write_tif(d / f"jrc_{layer}.tif", arr, g10, nodata=255)
        print(site_key, "jrc", layer, len(items), "tiles")


if __name__ == "__main__":
    for key in sys.argv[1:] or list(config.SITES):
        fetch(key)
