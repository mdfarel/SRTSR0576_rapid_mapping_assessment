"""Seasonal Sentinel-2 L2A median composites per site.

For each season the clearest scenes of every MGRS tile that touches the study
area are read onto the site grid, clouds, shadows and cirrus are removed with
the scene classification layer (SCL), and the per-pixel median of the remaining
observations is kept. The number of clear observations per pixel is written
alongside, so thin coverage stays visible instead of hidden.

Writes data/<site>/s2/<season>_<res>m.tif (surface reflectance, one band per
entry in config.S2_BANDS_20 or S2_BANDS_10) and <season>_<res>m_count.tif.
Scene bands are cached in data/<site>/cache/s2, so an interrupted run resumes.

    python -m pipeline.fetch_s2 [site ...]
"""
import sys
from collections import defaultdict

import numpy as np
from rasterio.enums import Resampling

import config
from pipeline.common import bbox_wgs84, cached, catalog, fill_cache, read_to_grid, site_dir, site_grid, write_tif

MAX_CLOUD = {"dry": 60, "wet": 90}  # scene-level cloud cover (%) accepted; the wet season needs more leeway
SCENES_PER_TILE = {("dry", 20): 8, ("wet", 20): 10, ("dry", 10): 6}  # (season, resolution in m): scenes kept
NATIVE_10M = {"B02", "B03", "B04", "B08"}  # bands stored at 10 m; read from their 20 m overview for a 20 m grid
QUANTIFICATION = 10000.0  # L2A digital numbers per unit of reflectance
BOA_OFFSET = -1000  # added to every digital number from processing baseline 04.00 on
OFFSET_BASELINE = 4.0
DOWNLOAD_THREADS = 6
ROW_BLOCK = 512  # the median is taken a block of rows at a time to limit memory


def cloud_cover(item):
    return item.properties["eo:cloud_cover"]


def clearest_per_tile(items, per_tile):
    """The `per_tile` least cloudy days of each MGRS tile, keeping the latest reprocessing of a tile-day."""
    by_tile = defaultdict(dict)
    for it in items:
        tile = it.properties.get("s2:mgrs_tile")
        day = it.datetime.date().isoformat()
        prev = by_tile[tile].get(day)
        if prev is None or it.id > prev.id:  # a later processing date sorts later in the id
            by_tile[tile][day] = it
    chosen = []
    for days in by_tile.values():
        ranked = sorted(days.values(), key=cloud_cover)
        chosen += ranked[:per_tile]
    return chosen


def pick_scenes(site_key, season, grid, per_tile):
    """Scenes for one season, in date order."""
    start, end = config.SITES[site_key]["seasons"][season]
    items = catalog().search(collections=["sentinel-2-l2a"], bbox=bbox_wgs84(grid),
                             datetime=f"{start}/{end}",
                             query={"eo:cloud_cover": {"lt": MAX_CLOUD[season]}}).items()
    chosen = clearest_per_tile(items, per_tile)
    return sorted(chosen, key=lambda i: i.datetime)


def boa_offset(item):
    """Planetary Computer serves L2A as delivered by ESA, so baseline 04.00+ scenes carry the -1000 DN
    offset. Checked on clear water: with the offset removed, a large clear reservoir reads B08 0.038 and
    B11 0.016 (dry season 2024); removing it twice would push water reflectance below zero."""
    baseline = float(item.properties.get("s2:processing_baseline", "0") or 0)
    return BOA_OFFSET if baseline >= OFFSET_BASELINE else 0


def read_band(item, band, grid, cache_dir):
    """One band of one scene on the grid as raw digital numbers (0 = nodata), cached on disk."""
    resampling = Resampling.nearest if band == "SCL" else Resampling.average  # SCL holds class codes
    overview = 0 if (grid.res >= 20 and band in NATIVE_10M) else None
    path = cache_dir / f"{item.id}_{band}_{int(grid.res)}m.npz"

    def download():
        arr, _ = read_to_grid(item.assets[band].href, grid, resampling, np.uint16, 0, overview)
        return arr
    return cached(path, download)


def clear_masks(items, grid, cache):
    """(scenes, rows, cols) stack, True where the scene classification marks the pixel as clear."""
    masks = []
    for it in items:
        scl = read_band(it, "SCL", grid, cache)
        masks.append(np.isin(scl, config.S2_VALID_SCL))
    return np.stack(masks)


def median_reflectance(items, band, valid, grid, cache):
    """Per-pixel median reflectance of one band over the clear observations; NaN where none."""
    stack = np.stack([read_band(it, band, grid, cache) for it in items])
    offsets = np.array([boa_offset(it) for it in items], np.float32)[:, None, None]
    out = np.full((grid.height, grid.width), np.nan, np.float32)
    for r0 in range(0, grid.height, ROW_BLOCK):
        rows = slice(r0, r0 + ROW_BLOCK)
        dn = stack[:, rows].astype(np.float32)
        ok = valid[:, rows] & (dn > 0)
        refl = np.where(ok, (dn + offsets) / QUANTIFICATION, np.nan)
        with np.errstate(all="ignore"):
            out[rows] = np.nanmedian(refl, axis=0)
    return out


def composite(site_key, season, res, bands):
    """Write the median composite and clear-observation count for one season and resolution."""
    out_dir = site_dir(site_key, "s2")
    out = out_dir / f"{season}_{res}m.tif"
    if out.exists():
        return
    grid = site_grid(site_key, res)
    g20 = site_grid(site_key, 20)
    cache = site_dir(site_key, "cache", "s2")
    items = pick_scenes(site_key, season, g20, SCENES_PER_TILE[(season, res)])
    print(f"{site_key} {season} {res}m: {len(items)} scenes", flush=True)

    jobs = [(it, band, grid, cache) for it in items for band in bands]
    jobs += [(it, "SCL", grid, cache) for it in items]
    fill_cache(read_band, jobs, DOWNLOAD_THREADS)

    valid = clear_masks(items, grid, cache)
    count = valid.sum(0).astype(np.uint8)
    result = np.full((len(bands), grid.height, grid.width), np.nan, np.float32)
    for bi, band in enumerate(bands):
        result[bi] = median_reflectance(items, band, valid, grid, cache)
        print(f"  {band} done", flush=True)

    write_tif(out, result, grid, nodata=np.nan, dtype="float32",
              tags={"bands": ",".join(bands), "scenes": ",".join(i.id for i in items)})
    write_tif(out_dir / f"{season}_{res}m_count.tif", count, grid, nodata=255, dtype="uint8")


def run(site_key):
    """Dry and wet composites at 20 m, and a dry composite at 10 m."""
    for season in ("dry", "wet"):
        composite(site_key, season, 20, config.S2_BANDS_20)
    composite(site_key, "dry", 10, config.S2_BANDS_10)


if __name__ == "__main__":
    for key in sys.argv[1:] or list(config.SITES):
        run(key)
