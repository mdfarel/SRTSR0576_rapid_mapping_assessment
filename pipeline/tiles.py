"""Web map layers: styled RGBA rasters packed as PMTiles, plus per-site summary JSON.

Each product is coloured as in the static maps (maps.py), warped to Web Mercator, cut into an
MBTiles pyramid with GDAL and repacked as one PMTiles archive, which the browser reads with HTTP
range requests. Needs the GDAL command-line tools (gdalwarp, gdal_translate, gdaladdo).

Inputs: the products in data/<site>/products/ and the cached Sentinel-2 dry-season composite.
Outputs: webgis/tiles/<site>/<layer>.pmtiles, webgis/data/<site>.json (the numbers the page shows),
webgis/data/<site>_water_features.csv, webgis/data/sites.json (the sites the page offers) and a copy of
the comparative CSV, webgis/data/waterfields_comparative_<year>.csv, which run.py refreshes once the
comparative dataset has been rebuilt.

    python -m pipeline.tiles [site ...]
"""
import json
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np
import rasterio
from matplotlib.colors import to_rgb
from pmtiles.convert import mbtiles_to_pmtiles
from rasterio.enums import ColorInterp
from rasterio.transform import from_origin

import config
from pipeline.common import Grid, write_tif
from pipeline.layers import Site, upsample
from pipeline.maps import (ACCENT, CHANGE_COLORS, CYCLE_COLORS, DENSITY_CLASSES, FIELD_COLORS, GRAIN_COLORS,
                           GRAIN_MIN_COHERENCE, JRC_CLASSES, MONTH_CLASSES, NDVI_CLASSES, SEASONAL, WATER_COLORS,
                           ndvi_change_class, radar_cannot_see_water, read_json)
from pipeline.water import RECTILINEAR

WEB = config.ROOT / "webgis"
MAX_ZOOM_PMTILES = 22  # highest zoom the PMTiles header may announce; the pyramid itself stops earlier

# The value layers are drawn in classes, each in the exact colour of its row in the web legend (webgis/app.js),
# so that the map, the legend and the identify tool, which names a pixel by its legend colour, agree. The
# classes (MONTH_CLASSES, DENSITY_CLASSES, JRC_CLASSES and NDVI_CLASSES) are those of the static maps, which
# maps.py defines.
COMPARATIVE_CSV = config.OUT / "dataset" / f"waterfields_comparative_{config.YEAR}.csv"


def rgba(shape):
    """An empty (fully transparent) RGBA image of the given rows x columns, bands first."""
    return np.zeros((4,) + shape, np.uint8)


def put(img, mask, colour, alpha=255):
    """Set the pixels under mask to one colour and opacity (0 to 255)."""
    r, g, b = (np.array(to_rgb(colour)) * 255).round().astype(np.uint8)
    img[0][mask], img[1][mask], img[2][mask], img[3][mask] = r, g, b, alpha


def put_classes(img, values, classes, alpha):
    """Colour the pixels whose value lies in each ((lowest, highest), colour, label) class; NaN takes no class."""
    for (lo, hi), colour, _ in classes:
        put(img, (values >= lo) & (values <= hi), colour, alpha)


# One function per layer. Each returns the image, its grid and the tile format: PNG keeps class
# colours exact, WEBP is smaller for the true-colour image.

def layer_truecolor(site):
    d = site.s2_composite("dry", 20)
    tc = np.stack([d["B04"], d["B03"], d["B02"]])
    img = rgba(tc.shape[1:])
    ok = np.isfinite(tc).all(0)
    img[:3] = (np.clip(np.nan_to_num(tc) / 0.16, 0, 1) ** 0.8 * 255).astype(np.uint8)
    img[3] = np.where(ok, 255, 0)
    return img, site.g20, "WEBP"


def layer_water(site):
    cls = site.load("water_class_10m")
    img = rgba(cls.shape)
    if not radar_cannot_see_water(site):  # seasonal inundation is radar water, left out where radar cannot see it
        seasonal = upsample(site.load("water_change_20m") == 2, 2)[: cls.shape[0], : cls.shape[1]]
        put(img, seasonal & (cls == 0), SEASONAL, 220)
    for k, c in WATER_COLORS.items():
        put(img, cls == k, c)
    return img, site.g10, "PNG"


def layer_cardinal(site):
    cls = site.load("water_class_10m")
    card = site.load("water_cardinal_10m") == 1
    rect = np.isin(cls, RECTILINEAR)
    img = rgba(cls.shape)
    put(img, (cls > 0) & ~rect, "#9fb8cf", 200)
    put(img, rect & ~card, "#3a5f7d")
    put(img, rect & card, ACCENT)
    return img, site.g10, "PNG"


def layer_fields(site):
    cls = site.load("field_class_20m")
    img = rgba(cls.shape)
    for k, c in FIELD_COLORS.items():
        put(img, cls == k, c, 230)
    return img, site.g20, "PNG"


def layer_water_change(site):
    wc = site.load("water_change_20m")
    img = rgba(wc.shape)
    for k, c in CHANGE_COLORS.items():
        put(img, wc == k, c)
    return img, site.g20, "PNG"


def layer_ndvi_change(site):
    """Greener and browner are changes beyond NDVI_CHANGE, as the statistics count them (maps.ndvi_change_class)."""
    ndvi = ndvi_change_class(site.load("ndvi_change_20m"))  # -1 where a season had no clear view
    img = rgba(ndvi.shape)
    for k, (colour, _) in enumerate(NDVI_CLASSES):
        put(img, ndvi == k, colour, 230)
    return img, site.g20, "PNG"


def layer_grain(site):
    dev = site.load("field_grain_deviation_100m")
    coh = site.load("field_grain_coherence_100m")
    img = rgba(dev.shape)
    has = np.isfinite(dev) & (np.nan_to_num(coh) >= GRAIN_MIN_COHERENCE)
    for (lo, hi), c in zip([(0, 7.5), (7.5, 15), (15, 30), (30, 46)], GRAIN_COLORS):
        put(img, has & (dev >= lo) & (dev < hi), c, 200)
    g = site.g10
    grid100 = Grid(g.epsg, from_origin(g.transform.c, g.transform.f, 100, 100), dev.shape[1], dev.shape[0], 100)
    return img, grid100, "PNG"


def layer_jrc(site):
    occ = site.static10("jrc_occurrence").astype(np.float32)  # 0 to 100%, 255 where there is no data
    img = rgba(occ.shape)
    put_classes(img, occ, JRC_CLASSES, 220)
    return img, site.g10, "PNG"


def layer_cycles(site):
    cyc = site.load("crop_cycles_20m")
    ag = site.load("field_class_20m") > 0
    img = rgba(cyc.shape)
    for k, c in CYCLE_COLORS.items():
        put(img, ag & (cyc == k), c, 230)
    return img, site.g20, "PNG"


def layer_months(site):
    months = site.load("s1_water_months_20m").astype(np.float32)
    img = rgba(months.shape)
    put_classes(img, months, MONTH_CLASSES, 230)
    return img, site.g20, "PNG"


def layer_density(site):
    dens = np.round(site.load("pond_density_20m"))
    img = rgba(dens.shape)
    put_classes(img, dens, DENSITY_CLASSES, 210)
    return img, site.g20, "PNG"


LAYERS = {
    "truecolor": layer_truecolor, "water": layer_water, "cardinal": layer_cardinal, "density": layer_density,
    "months": layer_months, "fields": layer_fields, "cycles": layer_cycles, "water_change": layer_water_change,
    "ndvi_change": layer_ndvi_change, "grain": layer_grain, "jrc": layer_jrc,
}


def clip_to_aoi(site, img, grid):
    """Make everything outside the study polygon transparent (the image passed in is not changed)."""
    step = int(round(grid.res / 10))
    inside = site.mask(10)[::step, ::step][: img.shape[1], : img.shape[2]]
    img = img.copy()
    img[3][: inside.shape[0], : inside.shape[1]][~inside] = 0
    return img


def sh(*args):
    """Run a command-line tool, stopping on any error."""
    subprocess.run([str(a) for a in args], check=True, capture_output=True)


def pack(img, grid, fmt, dest):
    """Warp an RGBA image to Web Mercator and write it as a PMTiles archive of fmt (PNG or WEBP) tiles."""
    # Class layers (PNG) are resampled with nearest neighbour so their colours stay exact.
    nearest = fmt == "PNG"
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        write_tif(tmp / "utm.tif", img, grid, dtype="uint8")
        with rasterio.open(tmp / "utm.tif", "r+") as ds:
            ds.colorinterp = [ColorInterp.red, ColorInterp.green, ColorInterp.blue, ColorInterp.alpha]
        sh("gdalwarp", "-q", "-t_srs", "EPSG:3857", "-r", "near" if nearest else "bilinear",
           "-co", "COMPRESS=DEFLATE", tmp / "utm.tif", tmp / "merc.tif")
        opts = ["-co", f"TILE_FORMAT={fmt}", "-co", "ZOOM_LEVEL_STRATEGY=UPPER",
                "-co", "RESAMPLING=NEAREST" if nearest else "RESAMPLING=BILINEAR"]
        if fmt == "WEBP":
            opts += ["-co", "QUALITY=82"]
        sh("gdal_translate", "-q", "-of", "MBTiles", *opts, tmp / "merc.tif", tmp / "t.mbtiles")
        sh("gdaladdo", "-q", "-r", "nearest" if nearest else "average", tmp / "t.mbtiles", "2", "4", "8", "16", "32")
        dest.parent.mkdir(parents=True, exist_ok=True)
        mbtiles_to_pmtiles(str(tmp / "t.mbtiles"), str(dest), MAX_ZOOM_PMTILES)


def rounded(points):
    """(lon, lat) points rounded to 5 decimals (about 1 m), as lists for JSON."""
    return [[round(lon, 5), round(lat, 5)] for lon, lat in points]


def summary(site):
    """Everything the web page shows for one site: its name, the study-area outline and the statistics files."""
    products = site.products
    return {
        "key": site.key, "name": site.meta["name"], "country": site.meta["country"],
        "centre": [site.meta["lon"], site.meta["lat"]],
        "boundary": rounded(site.aoi_lonlat),
        "aoi_source": site.meta["aoi_source"],
        "ftw": read_json(products / "ftw_stats.json"),
        "accuracy": read_json(products / "accuracy_stats.json"),
        "seasons": site.meta["seasons"], "year": config.YEAR,
        "water": site.load_json("water_stats"), "fields": site.load_json("field_stats"),
        "change": site.load_json("change_stats"),
    }


def write_manifest():
    """List the sites whose web data is complete, so the page only offers those."""
    data = WEB / "data"
    keys = [k for k in config.SITES if (data / f"{k}.json").exists() and (WEB / "tiles" / k / "water.pmtiles").exists()]
    (data / "sites.json").write_text(json.dumps(keys))
    return keys


def write_summary(key):
    """Rewrite only the per-site summary JSON (for example after new accuracy checks), without re-tiling."""
    data = WEB / "data"
    data.mkdir(parents=True, exist_ok=True)
    (data / f"{key}.json").write_text(json.dumps(summary(Site(key)), separators=(",", ":")))


def copy_comparative():
    """Copy the comparative CSV of the four sites into webgis/data/, where the statistics pane links to it."""
    if COMPARATIVE_CSV.exists():
        data = WEB / "data"
        data.mkdir(parents=True, exist_ok=True)
        (data / COMPARATIVE_CSV.name).write_bytes(COMPARATIVE_CSV.read_bytes())


def run(key):
    """All web layers of one site, its summary JSON and water-feature table, then the site list.

    The comparative CSV is copied too when it exists; run.py copies it again after rebuilding it.
    """
    site = Site(key)
    for name, make_layer in LAYERS.items():
        dest = WEB / "tiles" / key / f"{name}.pmtiles"
        img, grid, fmt = make_layer(site)
        pack(clip_to_aoi(site, img, grid), grid, fmt, dest)
        print(key, name, f"{dest.stat().st_size / 1e6:.1f} MB", flush=True)
    write_summary(key)
    src = site.products / "water_features.csv"
    (WEB / "data" / f"{key}_water_features.csv").write_bytes(src.read_bytes())
    copy_comparative()
    print("sites on the web map:", write_manifest(), flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
