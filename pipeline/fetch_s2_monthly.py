"""Monthly Sentinel-2 NDVI for the analysis year, for crop phenology.

For each month, the clearest scenes of every MGRS tile touching the study area
are read at 20 m (red, near infrared, scene classification), clouds and shadows
are removed, and the per-pixel median NDVI is kept with its count of clear
observations. Months without a clear view stay empty and are interpolated later
in pipeline.phenology.

Writes data/<site>/s2/ndvi_monthly_<year>.tif (12 bands, January to December)
and ndvi_monthly_<year>_count.tif. Scene bands share the cache of pipeline.fetch_s2.

    python -m pipeline.fetch_s2_monthly [site ...]
"""
import sys

import numpy as np

import config
from pipeline.common import bbox_wgs84, catalog, fill_cache, site_dir, site_grid, write_tif
from pipeline.fetch_s2 import QUANTIFICATION, boa_offset, clearest_per_tile, read_band
from pipeline.layers import ndi

PER_TILE = 4  # least cloudy scenes kept per MGRS tile and month
MAX_CLOUD = 80  # scene-level cloud cover (%) accepted
DOWNLOAD_THREADS = 8


def month_scenes(grid, month):
    """The clearest scenes of each tile for one month of config.YEAR."""
    y = config.YEAR
    start = f"{y}-{month:02d}-01"
    end = f"{y}-{month + 1:02d}-01" if month < 12 else f"{y + 1}-01-01"
    items = catalog().search(collections=["sentinel-2-l2a"], bbox=bbox_wgs84(grid),
                             datetime=f"{start}/{end}",
                             query={"eo:cloud_cover": {"lt": MAX_CLOUD}}).items()
    # the search end date is inclusive, so drop scenes from the first day of the next month
    in_month = [it for it in items if it.datetime.month == month]
    return clearest_per_tile(in_month, PER_TILE)


def scene_ndvi(item, grid, cache):
    """NDVI of one scene on the grid; NaN where the pixel is not clear."""
    scl = read_band(item, "SCL", grid, cache)
    red = read_band(item, "B04", grid, cache).astype(np.float32)
    nir = read_band(item, "B08", grid, cache).astype(np.float32)
    ok = np.isin(scl, config.S2_VALID_SCL) & (red > 0) & (nir > 0)
    off = boa_offset(item)
    r, n = (red + off) / QUANTIFICATION, (nir + off) / QUANTIFICATION
    return np.where(ok, ndi(n, r), np.nan).astype(np.float32)


def run(key):
    """Write the 12 monthly NDVI medians and their clear-observation counts for one site."""
    out = site_dir(key, "s2") / f"ndvi_monthly_{config.YEAR}.tif"
    if out.exists():
        return
    grid = site_grid(key, 20)
    cache = site_dir(key, "cache", "s2")
    ndvi = np.full((12, grid.height, grid.width), np.nan, np.float32)
    count = np.zeros((12, grid.height, grid.width), np.uint8)
    for m in range(1, 13):
        items = month_scenes(grid, m)
        jobs = [(it, band, grid, cache) for it in items for band in ("B04", "B08", "SCL")]
        fill_cache(read_band, jobs, DOWNLOAD_THREADS)
        if items:
            s = np.stack([scene_ndvi(it, grid, cache) for it in items])
            with np.errstate(all="ignore"):
                ndvi[m - 1] = np.nanmedian(s, axis=0)
            count[m - 1] = np.isfinite(s).sum(0)
        print(f"{key} month {m:02d}: {len(items)} scenes, clear median {int(np.median(count[m - 1]))}", flush=True)
    write_tif(out, ndvi, grid, nodata=np.nan, dtype="float32")
    write_tif(out.with_name(f"ndvi_monthly_{config.YEAR}_count.tif"), count, grid, dtype="uint8")


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
