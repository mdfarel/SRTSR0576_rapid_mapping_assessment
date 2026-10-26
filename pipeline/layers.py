"""Loading inputs onto the analysis grids, and small raster helpers.

Site gives every analysis step the same view of one site: its 10 m and 20 m
grids, the study-area mask, the cached inputs under data/<site>/ and the
products folder that the steps write to. The functions above it are the shared
raster arithmetic: normalised differences, Otsu thresholds, block averages and
slope. Nothing here downloads data; run the fetch steps first.
"""
import json

import numpy as np
from pyproj import CRS, Proj
from scipy import ndimage

import config
from pipeline.common import aoi_lonlat, aoi_mask, read_tif, site_dir, site_grid, write_tif

# Normalised differences are left undefined where the two reflectances sum to this or less (see ndi).
MIN_REFLECTANCE_SUM = 0.02
# Fewer finite values than this are too few for a stable threshold.
MIN_SAMPLE = 1000
# Below this the edge sample has no clear second peak, so Otsu is not trusted.
MIN_BIMODALITY = 0.3
# Width, in pixels, of the band either side of a water-land edge sampled by edge_otsu.
EDGE_BUFFER_PX = 2


def ndi(a, b):
    """Normalised difference of two reflectance bands.

    Undefined (NaN) where the two reflectances sum to 0.02 or less: after the
    processing-baseline offset is removed, very dark pixels can reach zero or
    below, and dividing by that sum gives values in the millions that would
    dominate any edge or orientation measure. Results are clipped to [-1, 1].
    """
    with np.errstate(all="ignore"):
        s = a + b
        out = np.where(s > MIN_REFLECTANCE_SUM, (a - b) / s, np.nan)
    return np.clip(out, -1.0, 1.0).astype(np.float32)


def otsu(values, lo, hi):
    """Otsu threshold of finite values, clipped to [lo, hi] so a unimodal scene cannot drift.

    The value that best splits a 256-bin histogram into two classes (largest
    between-class variance). The extreme 0.5 % at each end is clipped first so a
    few outliers cannot stretch the bins.
    """
    v = values[np.isfinite(values)]
    if v.size < MIN_SAMPLE:
        return (lo + hi) / 2
    v = np.clip(v, np.percentile(v, 0.5), np.percentile(v, 99.5))
    hist, edges = np.histogram(v, bins=256)
    centres = (edges[:-1] + edges[1:]) / 2
    # For every candidate threshold: pixel count and mean on each side.
    w0 = np.cumsum(hist)
    w1 = w0[-1] - w0
    m0 = np.cumsum(hist * centres)
    mu0 = m0 / np.maximum(w0, 1)
    mu1 = (m0[-1] - m0) / np.maximum(w1, 1)
    between = w0 * w1 * (mu0 - mu1) ** 2
    return float(np.clip(centres[np.argmax(between)], lo, hi))


def bimodality(values, t):
    """Between-class share of total variance at threshold t (Otsu's eta): near 1 is clearly two-peaked."""
    v = values[np.isfinite(values)]
    if v.size < MIN_SAMPLE:
        return float("nan")
    a, b = v[v <= t], v[v > t]
    if a.size == 0 or b.size == 0:
        return 0.0
    between = a.size * b.size / v.size ** 2 * (a.mean() - b.mean()) ** 2
    return float(between / max(v.var(), 1e-12))


def edge_otsu(index, inside, initial, lo, hi, fallback, water_high=True):
    """Otsu on the pixels around water-land edges (after Donchyts et al. 2016; Markert et al. 2020).

    A first split at `initial` locates the edges; the threshold is then set from the
    values within EDGE_BUFFER_PX of those edges, where water and land are sampled in
    similar numbers even when water covers only a few percent of the area.
    water_high=False is for indices where water is the low side (radar backscatter).
    Returns the threshold and the bimodality of the edge sample. Where the sample is
    not bimodal (little or no water, as in the Prambanan region) the conventional
    `fallback` value is returned instead, and the low bimodality is reported.
    """
    valid = inside & np.isfinite(index)
    if water_high:
        water = (index > initial) & valid
    else:
        water = (index < initial) & valid
    edge = water ^ ndimage.binary_erosion(water)
    zone = ndimage.binary_dilation(edge, iterations=EDGE_BUFFER_PX) & valid
    sample = index[zone]
    if sample.size < MIN_SAMPLE:
        sample = index[inside]
    t = otsu(sample, lo, hi)
    eta = round(bimodality(sample, t), 3)
    two_peaks = eta >= MIN_BIMODALITY  # False also when eta is NaN (too few pixels to judge)
    if not two_peaks:
        return float(fallback), eta
    return t, eta


def block_mean(arr, k):
    """Average k x k blocks (10 m to 20 m with k=2), ignoring NaN."""
    h, w = arr.shape
    a = arr[: h // k * k, : w // k * k].reshape(h // k, k, w // k, k).astype(np.float32)
    with np.errstate(all="ignore"):
        return np.nanmean(a, axis=(1, 3))


def upsample(arr, k):
    """Repeat each pixel k x k times (20 m to 10 m with k=2)."""
    return np.repeat(np.repeat(arr, k, axis=0), k, axis=1)


def slope_deg(dem, res):
    """Terrain slope in degrees from a DEM with square pixels of res metres."""
    gy, gx = np.gradient(dem.astype(np.float32), res)
    return np.degrees(np.arctan(np.hypot(gx, gy)))


def grid_convergence_deg(site_key):
    """Angle between UTM grid north and true north at the site centre."""
    s = config.SITES[site_key]
    factors = Proj(CRS.from_epsg(s["epsg"])).get_factors(s["lon"], s["lat"])
    return float(factors.meridian_convergence)


class Site:
    """Lazy access to every input layer of one site."""

    def __init__(self, key):
        self.key = key
        self.meta = config.SITES[key]
        self.g10 = site_grid(key, 10)
        self.g20 = site_grid(key, 20)
        self.static = site_dir(key, "static")
        self.s2 = site_dir(key, "s2")
        self.s1 = site_dir(key, "s1")
        self.products = site_dir(key, "products")
        self._masks = {}

    def mask(self, res):
        """True inside the study-area polygon."""
        if res not in self._masks:
            self._masks[res] = aoi_mask(self.key, res)
        return self._masks[res]

    @property
    def aoi_lonlat(self):
        return aoi_lonlat(self.key)

    def static10(self, name):
        """A static layer on the 10 m grid: dem, worldcover or jrc_<layer>."""
        return read_tif(self.static / f"{name}.tif")

    def s2_composite(self, season, res):
        """Seasonal Sentinel-2 median reflectance as {band name: 2-D array}."""
        arr = read_tif(self.s2 / f"{season}_{res}m.tif")
        bands = config.S2_BANDS_20 if res == 20 else config.S2_BANDS_10
        return {b: arr[i] for i, b in enumerate(bands)}

    def s2_count(self, season, res):
        """Number of clear Sentinel-2 observations behind each composite pixel."""
        return read_tif(self.s2 / f"{season}_{res}m_count.tif")

    def s1_monthly(self):
        """Monthly Sentinel-1 VV and VH medians in dB, each a (12, rows, cols) stack for January to December."""
        arr = read_tif(self.s1 / f"monthly_{config.YEAR}.tif")
        return arr[:12], arr[12:]

    def season_months(self, season):
        """Calendar months (1 to 12) covered by a season, e.g. [2, 3, 4] for the mainland dry season."""
        start, end = self.meta["seasons"][season]
        m0, m1 = int(start[5:7]), int(end[5:7])  # dates are "YYYY-MM-DD"
        return list(range(m0, m1 + 1))

    def save(self, name, arr, grid, nodata=None, dtype=None):
        write_tif(self.products / f"{name}.tif", arr, grid, nodata=nodata, dtype=dtype)

    def load(self, name):
        return read_tif(self.products / f"{name}.tif")

    def save_json(self, name, obj):
        (self.products / f"{name}.json").write_text(json.dumps(obj, indent=2))

    def load_json(self, name):
        return json.loads((self.products / f"{name}.json").read_text())
