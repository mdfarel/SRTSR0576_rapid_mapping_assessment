"""Grid, STAC access and cached raster reads shared by every step.

Every layer is resampled onto one square grid per site, so all products line up
pixel for pixel. Each remote read is cached to disk, so an interrupted run
resumes where it stopped.

Inputs: the study areas in config.py and Cloud-Optimised GeoTIFFs on the
Microsoft Planetary Computer. This module writes nothing by itself; the fetch
and analysis steps call its helpers.
"""
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import planetary_computer as pc
import pystac_client
import rasterio
from pyproj import Transformer
from rasterio.enums import Resampling
from rasterio.features import rasterize
from rasterio.transform import from_origin
from rasterio.vrt import WarpedVRT

import config

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"

# Grid edges are snapped to multiples of 20 m, so the 10 m grid nests exactly inside the 20 m grid.
SNAP_M = 20

# GDAL settings for reading remote files over HTTP: fetch only the byte ranges needed, and retry
# slow or failed requests instead of stopping the run.
GDAL_ENV = dict(
    GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR",
    CPL_VSIL_CURL_ALLOWED_EXTENSIONS=".tif,.TIF,.tiff,.zip,/content",
    GDAL_HTTP_MULTIRANGE="YES",
    GDAL_HTTP_MERGE_CONSECUTIVE_RANGES="YES",
    GDAL_HTTP_MAX_RETRY="5",
    GDAL_HTTP_RETRY_DELAY="2",
    GDAL_HTTP_TIMEOUT="90",
    GDAL_HTTP_CONNECTTIMEOUT="20",
    VSI_CACHE="TRUE",
)


@dataclass(frozen=True)
class Grid:
    """A north-up UTM pixel grid: EPSG code, affine transform, size in pixels and pixel size in metres."""
    epsg: int
    transform: object
    width: int
    height: int
    res: float

    @property
    def crs(self):
        return f"EPSG:{self.epsg}"

    @property
    def bounds(self):
        """(left, bottom, right, top) in metres."""
        x0, y1 = self.transform.c, self.transform.f
        return x0, y1 - self.height * self.res, x0 + self.width * self.res, y1

    @property
    def pixel_area_km2(self):
        return self.res * self.res / 1e6


def aoi_lonlat(site_key):
    """Study-area polygon as a closed lon/lat ring."""
    s = config.SITES[site_key]
    if s.get("aoi"):
        ring = list(s["aoi"])
    else:
        ring = square_lonlat(s)
    return ring + [ring[0]]


def square_lonlat(s):
    """Corners of the config.SIDE_M square centred on the site's lat/lon, in lon/lat (used for Angkor)."""
    to_utm = Transformer.from_crs("EPSG:4326", f"EPSG:{s['epsg']}", always_xy=True)
    to_ll = Transformer.from_crs(f"EPSG:{s['epsg']}", "EPSG:4326", always_xy=True)
    cx, cy = to_utm.transform(s["lon"], s["lat"])
    cx, cy = round(cx / SNAP_M) * SNAP_M, round(cy / SNAP_M) * SNAP_M
    h = config.SIDE_M / 2
    corners = [(cx - h, cy + h), (cx - h, cy - h), (cx + h, cy - h), (cx + h, cy + h)]
    return [to_ll.transform(x, y) for x, y in corners]


def aoi_utm(site_key):
    """Study-area polygon as a closed ring in the site's UTM zone (metres)."""
    s = config.SITES[site_key]
    to_utm = Transformer.from_crs("EPSG:4326", f"EPSG:{s['epsg']}", always_xy=True)
    return [to_utm.transform(lon, lat) for lon, lat in aoi_lonlat(site_key)]


def snap_down(v):
    """Round a UTM coordinate down to the 20 m lattice.

    Rounding to the millimetre first means a lon/lat round trip cannot shift a cell.
    """
    return np.floor(round(v, 3) / SNAP_M) * SNAP_M


def snap_up(v):
    """Round a UTM coordinate up to the 20 m lattice (see snap_down)."""
    return np.ceil(round(v, 3) / SNAP_M) * SNAP_M


def site_grid(site_key, res):
    """UTM grid enclosing the study area, snapped to 20 m so the 10 m grid is exactly twice as fine."""
    s = config.SITES[site_key]
    xs, ys = zip(*aoi_utm(site_key))
    x0, x1 = snap_down(min(xs)), snap_up(max(xs))
    y0, y1 = snap_down(min(ys)), snap_up(max(ys))
    width = int(round((x1 - x0) / res))
    height = int(round((y1 - y0) / res))
    return Grid(s["epsg"], from_origin(x0, y1, res, res), width, height, res)


def aoi_mask(site_key, res):
    """True inside the study-area polygon on the site grid."""
    grid = site_grid(site_key, res)
    shape = {"type": "Polygon", "coordinates": [aoi_utm(site_key)]}
    return rasterize([(shape, 1)], out_shape=(grid.height, grid.width), transform=grid.transform,
                     fill=0, dtype="uint8").astype(bool)


def bbox_wgs84(grid, pad_deg=0.0):
    """Lon/lat box [west, south, east, north] around a grid, widened by pad_deg, for catalogue searches."""
    to_ll = Transformer.from_crs(grid.crs, "EPSG:4326", always_xy=True)
    x0, y0, x1, y1 = grid.bounds
    xs, ys = to_ll.transform([x0, x1, x0, x1], [y0, y0, y1, y1])
    return [min(xs) - pad_deg, min(ys) - pad_deg, max(xs) + pad_deg, max(ys) + pad_deg]


def catalog():
    """The Planetary Computer STAC catalogue."""
    return pystac_client.Client.open(STAC_URL)


def read_to_grid(href, grid, resampling=Resampling.nearest, dtype=None, nodata=None, overview_level=None, sign=True):
    """Warp one remote COG onto the site grid. Pixels outside the source are nodata.

    Returns the array and the nodata value used. overview_level=0 reads the first
    internal overview (20 m for a 10 m band), which moves a quarter of the bytes
    when the target grid is coarser than the source. sign=False is for files
    outside the Planetary Computer, which need no access token.
    """
    url = pc.sign(href) if sign else href
    open_options = {}
    if overview_level is not None:
        open_options["overview_level"] = overview_level
    with rasterio.Env(**GDAL_ENV):
        with rasterio.open(url, **open_options) as src:
            nd = src.nodata if nodata is None else nodata
            with WarpedVRT(src, crs=grid.crs, transform=grid.transform, width=grid.width,
                           height=grid.height, resampling=resampling,
                           src_nodata=nd, nodata=nd) as vrt:
                arr = vrt.read(1)
    if dtype is not None:
        arr = arr.astype(dtype, copy=False)
    return arr, nd


def cached(path, fn):
    """Return the array stored at path, or compute it with fn() and store it."""
    path = Path(path)
    if path.exists():
        return np.load(path)["a"]
    arr = fn()
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.stem + ".partial.npz")  # written first, so a killed run leaves no broken cache file
    np.savez_compressed(tmp, a=arr)
    tmp.replace(path)
    return arr


def fill_cache(read, jobs, threads):
    """Call read(*job) for every job on several threads, only to fill the disk cache.

    Each array is dropped as soon as it is read, so memory stays low; the caller
    then reads the arrays back from the cache one band at a time.
    """
    def read_and_drop(job):
        read(*job)

    with ThreadPoolExecutor(threads) as pool:
        for _ in pool.map(read_and_drop, jobs):
            pass  # wait for every read; an error in any of them stops the run here


def write_tif(path, arr, grid, nodata=None, dtype=None, tags=None):
    """Write a 2-D array, or a (bands, rows, cols) stack, as a compressed tiled GeoTIFF on the grid."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = arr[None]
    dtype = dtype or arr.dtype
    is_integer = np.issubdtype(np.dtype(dtype), np.integer)
    profile = dict(driver="GTiff", width=grid.width, height=grid.height, count=arr.shape[0],
                   dtype=dtype, crs=grid.crs, transform=grid.transform, nodata=nodata,
                   compress="deflate", tiled=True, blockxsize=512, blockysize=512,
                   predictor=2 if is_integer else 3)  # 2 suits integers, 3 suits floating point
    with rasterio.open(path, "w", **profile) as dst:
        dst.write(arr.astype(dtype, copy=False))
        if tags:
            dst.update_tags(**tags)


def read_tif(path):
    """Read a GeoTIFF: a 2-D array for one band, a (bands, rows, cols) stack otherwise."""
    with rasterio.open(path) as src:
        return src.read(1) if src.count == 1 else src.read()


def site_dir(site_key, *parts):
    """The folder data/<site>/<parts...>, created if missing."""
    p = config.DATA.joinpath(site_key, *parts)
    p.mkdir(parents=True, exist_ok=True)
    return p
