"""Monthly Sentinel-1 RTC backscatter (VV, VH) per site for the analysis year.

Radar sees through cloud, so it carries the wet-season water mapping and the
rice-flooding signal where optical images are mostly cloud. Uses the
radiometrically terrain-corrected (gamma0) product, read onto the 20 m site
grid and summarised as monthly medians in dB.

Writes data/<site>/s1/monthly_<year>.tif with 24 bands: vv_01 to vv_12, then
vh_01 to vh_12. Every acquisition of the year is used. Scenes are cached in
data/<site>/cache/s1, so an interrupted run resumes.

    python -m pipeline.fetch_s1 [site ...]
"""
import sys

import numpy as np
from rasterio.enums import Resampling

import config
from pipeline.common import bbox_wgs84, cached, catalog, fill_cache, read_to_grid, site_dir, site_grid, write_tif

POLS = ("vv", "vh")
S1_NODATA = -32768.0
DOWNLOAD_THREADS = 6


def scenes(grid):
    """Every Sentinel-1 RTC acquisition over the grid in config.YEAR, in date order."""
    y = config.YEAR
    items = catalog().search(collections=["sentinel-1-rtc"], bbox=bbox_wgs84(grid),
                             datetime=f"{y}-01-01/{y}-12-31").items()
    return sorted(items, key=lambda i: i.datetime)


def read_pol(item, pol, grid, cache_dir):
    """One polarisation of one scene on the grid as linear gamma0, NaN where there is no valid value; cached."""
    path = cache_dir / f"{item.id}_{pol}.npz"

    def download():
        arr, nd = read_to_grid(item.assets[pol].href, grid, Resampling.average, np.float32, S1_NODATA,
                               overview_level=0)
        arr[(arr == nd) | ~np.isfinite(arr) | (arr <= 0)] = np.nan
        return arr
    return cached(path, download)


def monthly_median_db(items, pol, month, grid, cache):
    """Per-pixel median backscatter of one month in dB; NaN everywhere if the month has no acquisition."""
    month_items = [it for it in items if it.datetime.month == month]
    if not month_items:
        return np.full((grid.height, grid.width), np.nan, np.float32)
    stack = np.stack([read_pol(it, pol, grid, cache) for it in month_items])
    with np.errstate(all="ignore"):
        db = 10 * np.log10(np.nanmedian(stack, axis=0))
    return db.astype(np.float32)


def run(site_key):
    """Write the monthly VV and VH medians for one site."""
    out_dir = site_dir(site_key, "s1")
    out = out_dir / f"monthly_{config.YEAR}.tif"
    if out.exists():
        return
    grid = site_grid(site_key, 20)
    cache = site_dir(site_key, "cache", "s1")
    items = scenes(grid)
    print(f"{site_key} S1: {len(items)} scenes", flush=True)
    jobs = [(it, pol, grid, cache) for it in items for pol in POLS]
    fill_cache(read_pol, jobs, DOWNLOAD_THREADS)

    bands, names = [], []
    for pol in POLS:
        for m in range(1, 13):
            bands.append(monthly_median_db(items, pol, m, grid, cache))
            names.append(f"{pol}_{m:02d}")
        print(f"  {pol} done", flush=True)
    write_tif(out, np.stack(bands), grid, nodata=np.nan, dtype="float32",
              tags={"bands": ",".join(names), "scenes": ",".join(i.id for i in items)})


if __name__ == "__main__":
    for key in sys.argv[1:] or list(config.SITES):
        run(key)
