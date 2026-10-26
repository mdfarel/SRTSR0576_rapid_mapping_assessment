"""Crop phenology from monthly Sentinel-2 NDVI: cropping cycles per year.

Months without a clear view are filled by linear interpolation between the
nearest clear months, the series is lightly smoothed, and peaks are counted on
the year treated as a circle (a crop that peaks across New Year counts once).
A peak counts as a cropping cycle when NDVI reaches 0.45 and rises at least
0.15 above both neighbouring troughs (the lowest NDVI in the three months
before and in the three months after the peak). In a single-monsoon climate,
two or more cycles a year generally mean irrigation or residual-moisture
cropping, which is the signal radar loses over Java's small, asynchronously
planted plots.

Inputs: the monthly NDVI stack and its clear-observation counts from
pipeline.fetch_s2_monthly (data/<site>/s2/ndvi_monthly_<year>.tif and _count.tif).
Outputs (20 m GeoTIFFs in data/<site>/products): crop_cycles_20m (255 where not
judged), ndvi_peak_month_20m, ndvi_amplitude_20m and ndvi_valid_months_20m.

    python -m pipeline.phenology [site ...]
"""
import sys

import numpy as np

import config
from pipeline.common import read_tif
from pipeline.layers import Site

PEAK_MIN = 0.45        # NDVI a peak must reach to count as a crop
PROMINENCE = 0.15      # rise of the peak above both neighbouring troughs
WINDOW = 3             # months either side when looking for troughs
MIN_CLEAR_MONTHS = 6   # fewer clear months than this and the pixel is not judged
NOT_JUDGED = 255       # nodata value of the crop-cycle raster


def fill_gaps(v):
    """Linear interpolation over months along axis 0, holding the ends."""
    v = v.copy()
    n = v.shape[0]
    ok = np.isfinite(v)
    idx = np.arange(n)[:, None, None]
    # for every month, the index of the last clear month at or before it (-1 if none)
    # and of the first clear month at or after it (n if none)
    prev_i = np.where(ok, idx, -1)
    next_i = np.where(ok, idx, n)
    prev_i = np.maximum.accumulate(prev_i, axis=0)
    next_i = np.minimum.accumulate(next_i[::-1], axis=0)[::-1]
    prev_v = np.take_along_axis(v, np.clip(prev_i, 0, n - 1), axis=0)
    next_v = np.take_along_axis(v, np.clip(next_i, 0, n - 1), axis=0)
    has_p, has_n = prev_i >= 0, next_i < n
    # between two clear months interpolate; before the first or after the last, hold the nearest value
    w = np.where(has_p & has_n, (idx - prev_i) / np.maximum(next_i - prev_i, 1), 0.0)
    filled = np.where(has_p & has_n, prev_v + w * (next_v - prev_v), np.where(has_p, prev_v, next_v))
    return np.where(ok, v, filled).astype(np.float32)


def cycles(ndvi):
    """Number of crop peaks in the year at each pixel, and the smoothed monthly NDVI.

    The months wrap round (December is next to January), so a peak across New Year counts once.
    """
    # 1-2-1 smoothing over neighbouring months
    s = 0.25 * np.roll(ndvi, 1, 0) + 0.5 * ndvi + 0.25 * np.roll(ndvi, -1, 0)
    n = s.shape[0]
    count = np.zeros(s.shape[1:], np.uint8)
    for i in range(n):
        v = s[i]
        is_max = (v >= s[(i - 1) % n]) & (v > s[(i + 1) % n])
        trough_before = np.min([s[(i - k) % n] for k in range(1, WINDOW + 1)], axis=0)
        trough_after = np.min([s[(i + k) % n] for k in range(1, WINDOW + 1)], axis=0)
        prom = v - np.maximum(trough_before, trough_after)
        count += (is_max & (v >= PEAK_MIN) & (prom >= PROMINENCE)).astype(np.uint8)
    return count, s


def run(key):
    """Crop-cycle count, NDVI peak month, amplitude and clear-month count for one site."""
    site = Site(key)
    ndvi = read_tif(site.s2 / f"ndvi_monthly_{config.YEAR}.tif").astype(np.float32)
    counts = read_tif(site.s2 / f"ndvi_monthly_{config.YEAR}_count.tif")
    valid_months = (counts > 0).sum(0).astype(np.uint8)
    filled = fill_gaps(ndvi)
    n_cycles, smooth = cycles(np.nan_to_num(filled, nan=0.0))
    inside = site.mask(20)
    n_cycles[(valid_months < MIN_CLEAR_MONTHS) | ~inside] = NOT_JUDGED  # too few clear months, or outside the area
    peak_month = (np.argmax(smooth, axis=0) + 1).astype(np.uint8)
    amplitude = (smooth.max(0) - smooth.min(0)).astype(np.float32)
    site.save("crop_cycles_20m", n_cycles, site.g20, nodata=NOT_JUDGED, dtype="uint8")
    site.save("ndvi_peak_month_20m", np.where(inside, peak_month, 0), site.g20, nodata=0, dtype="uint8")
    site.save("ndvi_amplitude_20m", np.where(inside, amplitude, np.nan), site.g20, nodata=np.nan, dtype="float32")
    site.save("ndvi_valid_months_20m", valid_months, site.g20, dtype="uint8")
    judged = n_cycles != NOT_JUDGED
    print(key, "clear months median", int(np.median(valid_months[inside])), "judged",
          f"{100 * judged.sum() / inside.sum():.0f}%",
          {c: int((n_cycles == c).sum()) for c in range(4)}, flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
