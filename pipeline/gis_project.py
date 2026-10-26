"""A QGIS project per study area: our interpretation, with the open data it was made from.

Each site gets a folder out/gis/<site>/ holding <site>.qgz, its layers and README.txt, and the folder
is zipped to out/gis_zip/<site>_waterfields_qgis.zip. The project refers to its layers by relative
paths, so the unzipped folder opens in QGIS on any computer. The layers sit in two groups:

    Interpretation (results)  the study area and the water, field and change products, in the colours
                              and with the legends of the static maps (pipeline/maps.py)
    Open data (inputs)        the datasets named in the brief, clipped to the study area: Sentinel-2
                              true colour and Sentinel-1 VV for the dry and the wet season, JRC water
                              occurrence, ESA WorldCover 2021 and a Copernicus DEM hillshade, with
                              Fields of The World 2024 as the reference field map

Everything outside the study-area polygon is no data. The rasters are compressed so that each zip
stays small: DEFLATE for classes and values, some values stored as small integers with a scale
factor that QGIS applies on reading (NDVI change to 0.01, radar to 0.25 dB), and JPEG for the true
colour and the radar images, which are for viewing. The layers alone are kept under MAX_ZIP_MIB:
where the Fields of The World polygons would take them over, their outlines go in as a 10 m raster
instead, and README.txt says where the polygons can be read.

The project also holds one A3 print layout per analysis (SHEETS): the map over a faded dry-season
Sentinel-2 true colour, with a graticule, the study-area outline, a legend whose labels carry the
areas or counts of the products, a north arrow, a scale bar, a short analysis, the data sources and
their dates, the coordinate reference system and a credit line. Each layout draws a map theme of
the same name, which can also be shown in the map canvas; the layouts reach the styles they need
(the faded background, single-colour water) as named styles of the project's layers, so the layers
keep their own styles in the Layers panel. The layouts are exported to layouts/NN_<name>.pdf in the
site folder, with the raster images at LAYOUT_DPI, and to one <site>_map_atlas.pdf, and the zip
holds them too. The PDFs can take a zip over MAX_ZIP_MIB; the build then says by how much.

Inputs: data/<site>/products/, the cached data in data/<site>/s2/, s1/ and static/, and the
comparative CSV in out/dataset/ when it exists (the layouts read some of their numbers from it).
Needs QGIS's own Python (PyQGIS), which the rest of the pipeline does not. PyMuPDF (fitz), when
importable, recompresses the images in the PDFs; it, pypdf or PyPDF2 joins them into the atlas.
The layout text is set in Segoe UI where its files can be found, and wrapped and boxed so that it
also fits in DejaVu Sans, which QGIS takes instead on a computer without Segoe UI.

    python3 -m pipeline.gis_project [site ...]
"""
import csv
import json
import os
import re
import shutil
import sys
import warnings
import zipfile
from dataclasses import dataclass
from datetime import date

import numpy as np
import rasterio
import shapely
from qgis.core import (Qgis, QgsApplication, QgsCategorizedSymbolRenderer, QgsColorRampLegendNodeSettings,
                       QgsColorRampShader, QgsContrastEnhancement, QgsCoordinateReferenceSystem,
                       QgsCoordinateTransformContext, QgsFeature, QgsField, QgsGeometry, QgsHillshadeRenderer,
                       QgsLayerTreeLayer, QgsLayoutExporter, QgsLayoutItemLabel, QgsLayoutItemLegend,
                       QgsLayoutItemMap, QgsLayoutItemMapGrid, QgsLayoutItemPicture, QgsLayoutItemScaleBar,
                       QgsLayoutMeasurement, QgsLayoutPoint, QgsLayoutRenderContext, QgsLayoutSize, QgsLayoutUtils,
                       QgsLegendRenderer, QgsLegendStyle, QgsLineSymbol, QgsMapLayerLegendUtils,
                       QgsMapThemeCollection, QgsMultiBandColorRenderer, QgsPalettedRasterRenderer, QgsPrintLayout,
                       QgsProject, QgsProjectMetadata, QgsRasterLayer, QgsRasterRange, QgsRasterShader, QgsRectangle,
                       QgsReferencedRectangle, QgsRenderContext, QgsRendererCategory, QgsSingleBandGrayRenderer,
                       QgsSingleBandPseudoColorRenderer, QgsSingleSymbolRenderer, QgsSymbol, QgsSymbolLayerUtils,
                       QgsTextFormat, QgsTextRenderer, QgsVectorFileWriter, QgsVectorLayer)
from qgis.PyQt.QtCore import QMetaType, QRectF, Qt
from qgis.PyQt.QtGui import QColor, QFont, QFontDatabase
from rasterio.features import rasterize, shapes
from scipy import ndimage

import config
from pipeline import ftw
from pipeline.common import Grid, aoi_lonlat, aoi_utm, read_tif
from pipeline.dataset import passes_rotation_test
from pipeline.layers import Site
from pipeline.change import NDVI_CHANGE
from pipeline.maps import (ACCENT, AG_SHARE_CLASSES, CHANGE_COLORS, CHANGE_LABELS, CONCENTRATED, CYCLE_COLORS,
                           CYCLE_LABELS, DENSITY_CLASSES, FIELD_COLORS, FIELD_LABELS, GRAIN_COLORS,
                           GRAIN_MIN_COHERENCE, INK, JRC_CLASSES, MIN_CLASS_KM2, MONTH_CLASSES, MUTED, NDVI_CLASSES,
                           NDVI_STRONG, SPARSE_KM2, SWATCH_OUTLINE, TRAJECTORY, TRAJECTORY_CODES, WATER_COLORS,
                           WATER_FOR_REFERENCE, WATER_LABELS, field_classes_shown, km2_text, lower_first,
                           ndvi_change_class, radar_cannot_see_water, read_json)
from pipeline.water import CLASSES, RECTILINEAR

ZIP_DIR = config.OUT / "gis_zip"
MAX_ZIP_MIB = 25            # size limit of a zip, so that each site can be sent as one file
RESULTS = "Interpretation (results)"
INPUTS = "Open data (inputs)"

# Sentinel-2 true colour: reflectance / WHITE_REFLECTANCE, raised to TRUE_COLOUR_GAMMA, as on the maps
WHITE_REFLECTANCE = 0.18
TRUE_COLOUR_GAMMA = 0.8
TRUE_COLOUR_QUALITY = 85    # JPEG quality
# Sentinel-1 VV, stored as 1 to 255 in steps of VV_STEP_DB above VV_FLOOR_DB, drawn black to white over VV_SHOWN_DB
VV_FLOOR_DB = -40.0
VV_STEP_DB = 0.25
VV_SHOWN_DB = (-25.0, 0.0)
RADAR_QUALITY = 85          # JPEG quality
NDVI_STEP = 0.01            # NDVI change stored in hundredths
OPEN_END = 1e6              # top of a class with no highest value, above any value a layer holds
DEM_Z_FACTOR = 2            # hillshade vertical exaggeration, as on the maps
FTW_PRECISION_M = 1.0       # Fields of The World vertices rounded to 1 m; they follow 10 m pixel edges
FTW_COLOUR = "#2f6f3a"
# Crop cycles in the package: 0 to 3 cycles, TOO_FEW_MONTHS where there were too few clear months to
# count them, NOT_SHOWN off agricultural land and outside the study area
TOO_FEW_MONTHS = 4
NOT_SHOWN = 255

# ESA WorldCover 2021 classes and their official colours
WORLDCOVER = {
    10: ("#006400", "Tree cover"), 20: ("#ffbb22", "Shrubland"), 30: ("#ffff4c", "Grassland"),
    40: ("#f096ff", "Cropland"), 50: ("#fa0000", "Built-up"), 60: ("#b4b4b4", "Bare or sparse vegetation"),
    70: ("#f0f0f0", "Snow and ice"), 80: ("#0064c8", "Permanent water bodies"),
    90: ("#0096a0", "Herbaceous wetland"), 95: ("#00cf75", "Mangroves"), 100: ("#fae6a0", "Moss and lichen"),
}
# Products copied as they are, apart from no data outside the study area: (name, resolution, no-data value)
COPIED = [
    ("water_class_10m", 10, 0), ("water_cardinal_10m", 10, 0), ("water_change_20m", 20, 0),
    ("s1_water_months_20m", 20, 255), ("field_class_20m", 20, 0),
]
# Columns of water_features.csv carried by the water-body polygons, with their types (lon and lat are left out)
WATER_FIELDS = [
    ("id", int), ("class", str), ("area_ha", float), ("length_m", int), ("mean_width_m", float),
    ("elongation", float), ("rectangularity", float), ("rectilinearity", float),
    ("deviation_from_cardinal_deg", float), ("encloses_land_ratio", float), ("azimuth_deg", float), ("cardinal", int),
]
QGIS_TYPES = {int: QMetaType.Type.Int, float: QMetaType.Type.Double, str: QMetaType.Type.QString}


# Files for the package

def write_raster(path, arr, grid, dtype, nodata=None, mask=None, jpeg_quality=None, scale=None, offset=0.0,
                 nbits=None):
    """Write a 2-D array or a (bands, rows, cols) stack as a tiled GeoTIFF on a site grid.

    DEFLATE unless jpeg_quality is given. JPEG cannot keep a no-data value exact, so a JPEG file
    keeps `mask` (True where there is data) as its internal mask instead. With scale, readers see
    value * scale + offset.
    """
    arr = np.asarray(arr)
    if arr.ndim == 2:
        arr = arr[None]
    profile = dict(driver="GTiff", width=grid.width, height=grid.height, count=arr.shape[0], dtype=dtype,
                   crs=grid.crs, transform=grid.transform, nodata=nodata, tiled=True, blockxsize=512,
                   blockysize=512)
    if jpeg_quality:
        profile.update(compress="jpeg", jpeg_quality=jpeg_quality)
        if arr.shape[0] == 3:
            profile["photometric"] = "ycbcr"
    else:
        is_integer = np.issubdtype(np.dtype(dtype), np.integer)
        profile.update(compress="deflate", zlevel=9, predictor=2 if is_integer else 3)
    if nbits:
        profile.update(nbits=nbits, predictor=1)
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(path, "w", **profile) as dst:
            dst.write(arr.astype(dtype, copy=False))
            if mask is not None:
                dst.write_mask(mask.astype(np.uint8) * 255)
            if scale:
                dst.scales = (scale,) * arr.shape[0]
                dst.offsets = (offset,) * arr.shape[0]


def season_months(site, season):
    """Array indices (0 to 11) of the months of a season."""
    return [m - 1 for m in site.season_months(season)]


def season_text(site, season):
    """The dates of a season in words, for example '1 February to 30 April 2024'."""
    start, end = (date.fromisoformat(d) for d in site.meta["seasons"][season])
    return f"{start.day} {start:%B} to {end.day} {end:%B} {end.year}"


def write_inputs(site, out):
    """The brief's datasets on the site grids, with no data outside the study area."""
    in10, in20 = site.mask(10), site.mask(20)
    for season in ("dry", "wet"):
        d = site.s2_composite(season, 20)
        rgb = np.stack([d["B04"], d["B03"], d["B02"]])
        ok = np.isfinite(rgb).all(axis=0) & in20
        rgb8 = np.round(np.clip(np.nan_to_num(rgb) / WHITE_REFLECTANCE, 0, 1) ** TRUE_COLOUR_GAMMA * 255)
        write_raster(out / f"s2_truecolour_{season}_20m.tif", rgb8, site.g20, "uint8", mask=ok,
                     jpeg_quality=TRUE_COLOUR_QUALITY)

    # radar for each season as the water analysis reads it (pipeline/water.py): the median of the
    # dry-season monthly medians, and the lowest wet-season monthly median, which is the wettest state
    vv, _ = site.s1_monthly()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # pixels never observed stay NaN
        radar = {"dry": np.nanmedian(vv[season_months(site, "dry")], axis=0),
                 "wet": np.nanmin(vv[season_months(site, "wet")], axis=0)}
    for season, db in radar.items():
        code = np.clip(np.round((np.nan_to_num(db) - VV_FLOOR_DB) / VV_STEP_DB), 1, 255)
        write_raster(out / f"s1_vv_{season}_20m.tif", code, site.g20, "uint8", mask=np.isfinite(db) & in20,
                     jpeg_quality=RADAR_QUALITY, scale=VV_STEP_DB, offset=VV_FLOOR_DB)

    # JRC occurrence is a 30 m product, so keeping every other cell of its 10 m copy loses nothing that matters
    jrc = site.static10("jrc_occurrence")[::2, ::2]
    write_raster(out / "jrc_occurrence_20m.tif", np.where(in20, jrc, 255), site.g20, "uint8", nodata=255)
    # the JRC transition classes (0 to 10), the same 30 m product kept in the same way
    transitions = site.static10("jrc_transitions")[::2, ::2]
    write_raster(out / "jrc_transitions_20m.tif", np.where(in20, transitions, 255), site.g20, "uint8", nodata=255)
    wc = site.static10("worldcover")
    write_raster(out / "worldcover_2021_10m.tif", np.where(in10, wc, 0), site.g10, "uint8", nodata=0)
    dem = site.static10("dem")[::2, ::2]
    write_raster(out / "dem_20m.tif", np.where(in20 & (dem > -1000), np.round(dem), -32768), site.g20, "int16",
                 nodata=-32768)


def write_ftw_polygons(site, out, fields):
    """Fields of The World polygons inside the study area (chosen as in pipeline/ftw.py), with area and confidence."""
    utm, conf = fields
    geoms = shapely.set_precision(utm, FTW_PRECISION_M)
    keep = np.isin(shapely.get_type_id(geoms), (3, 6))  # polygons and multipolygons; drops any that collapsed
    rows = [(g, [round(float(a) / 1e4, 3), None if np.isnan(c) else round(float(c), 1)])
            for g, a, c in zip(geoms[keep], shapely.area(utm[keep]), conf[keep])]
    write_polygons(out / "ftw_fields_2024.gpkg", "ftw_fields_2024", site.g10.crs,
                   [("area_ha", float), ("confidence", float)], rows)


def write_ftw_outlines(site, out, fields):
    """The outlines of the same Fields of The World polygons as a 1-bit raster on the 10 m grid."""
    g = site.g10
    outline = rasterize(((line, 1) for line in shapely.boundary(fields[0])), out_shape=(g.height, g.width),
                        transform=g.transform, fill=0, dtype="uint8")
    write_raster(out / "ftw_field_outlines_10m.tif", outline * site.mask(10), g, "uint8", nodata=0, nbits=1)


def ftw_source(key):
    """The folder in which Fields of The World keeps the polygons of a site's country (the file pipeline/ftw.py
    reads is the one GeoParquet file in it), named by its country code only."""
    cc = ftw.COUNTRY[key].split("/")[0]
    return f"{ftw.BASE}/admin:country_code={cc}/"


def in_steps(values, step, classify):
    """Values as whole numbers of `step`, each kept in the class that classify (as the maps classify) gives it.

    A value is rounded to the nearest step, unless that would carry it across a class break; it then
    goes to the step on its own side of the break, at most one step from the computed value. In this
    way every pixel of a stored copy takes the class it has on the maps. NaN stays NaN.
    """
    per_unit = round(1 / step)
    code = np.round(values * per_unit)
    crossed = classify(code / per_unit) != classify(values)
    return np.where(crossed, code + np.sign(values * per_unit - code), code)


def ag_share_class(share):
    """The AG_SHARE_CLASSES class of each share (1 to 4), 0 below the lowest class or where there is no value."""
    return np.digitize(np.nan_to_num(share, nan=-1.0), [lo for lo, _, _ in AG_SHARE_CLASSES])


def write_results(site, out):
    """The products, with no data outside the study area; scaled or masked copies where noted."""
    masks = {10: site.mask(10), 20: site.mask(20)}
    for name, res, nodata in COPIED:
        arr = site.load(name)
        write_raster(out / f"{name}.tif", np.where(masks[res], arr, nodata), site.g10 if res == 10 else site.g20,
                     arr.dtype, nodata=nodata)
    g20, in20 = site.g20, masks[20]

    dens = site.load("pond_density_20m")
    write_raster(out / "pond_density_20m.tif", np.where(in20, dens, np.nan), g20, "float32", nodata=np.nan)
    # crop cycles on agricultural land only, as on the map
    ag = (site.load("field_class_20m") > 0) & in20
    cyc = site.load("crop_cycles_20m")
    cyc = np.where(cyc > 3, TOO_FEW_MONTHS, cyc)
    write_raster(out / "crop_cycles_20m.tif", np.where(ag, cyc, NOT_SHOWN), g20, "uint8", nodata=NOT_SHOWN)
    # share of agricultural land in 1 km, as a whole percentage
    share = site.load("ag_intensity_20m")
    pct = in_steps(share, 0.01, ag_share_class)
    write_raster(out / "ag_share_1km_20m.tif", np.where(np.isfinite(share) & in20, pct, 255), g20, "uint8",
                 nodata=255)
    dn = site.load("ndvi_change_20m")
    code = np.where(np.isfinite(dn) & in20, np.clip(in_steps(dn, NDVI_STEP, ndvi_change_class), -100, 100), -128)
    write_raster(out / "ndvi_change_20m.tif", code, g20, "int8", nodata=-128, scale=NDVI_STEP)

    # dominant field-boundary direction, only where the boundaries are directional enough, as on the map
    with rasterio.open(site.products / "field_grain_deviation_100m.tif") as src:
        g100 = Grid(site.meta["epsg"], src.transform, src.width, src.height, 100)
    study = {"type": "Polygon", "coordinates": [aoi_utm(site.key)]}
    inside = rasterize([(study, 1)], out_shape=(g100.height, g100.width), transform=g100.transform, fill=0,
                       dtype="uint8").astype(bool)
    dev = site.load("field_grain_deviation_100m")
    show = np.isfinite(dev) & (site.load("field_grain_coherence_100m") >= GRAIN_MIN_COHERENCE) & inside
    write_raster(out / "field_edge_direction_100m.tif", np.where(show, dev, np.nan), g100, "float32", nodata=np.nan)

    write_water_bodies(site, out)
    shutil.copy(site.products / "water_features.csv", out / "water_features.csv")
    area = {"type": "Feature", "geometry": {"type": "Polygon", "coordinates": [aoi_lonlat(site.key)]},
            "properties": {"name": site.meta["name"], "source": site.meta["aoi_source"]}}
    (out / "study_area.geojson").write_text(json.dumps({"type": "FeatureCollection", "features": [area]}))


def write_water_bodies(site, out):
    """One polygon per water body with its row of water_features.csv.

    The bodies are labelled as in pipeline/water.py (8-connected groups of dry-season water), so a
    body's label is its id in the CSV; the check on areas makes sure that the two still agree.
    """
    lab, n = ndimage.label(site.load("water_dry_10m") > 0, np.ones((3, 3)))
    with open(site.products / "water_features.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    pixels = np.bincount(lab.ravel())[1:]
    area_ha = np.array([float(r["area_ha"]) for r in rows])
    assert len(rows) == n and np.allclose(pixels * 0.01, area_ha, atol=0.006), "water_features.csv does not match"
    code = {name: k for k, name in CLASSES.items()}
    features = []
    for geom, value in shapes(lab.astype(np.int32), mask=lab > 0, connectivity=8, transform=site.g10.transform):
        row = rows[int(value) - 1]
        values = [row[name] if kind is str else kind(float(row[name])) for name, kind in WATER_FIELDS]
        features.append((shapely.geometry.shape(geom), [code[row["class"]]] + values))
    write_polygons(out / "water_bodies.gpkg", "water_bodies", site.g10.crs, [("class_code", int)] + WATER_FIELDS,
                   features)


def write_polygons(path, name, crs, fields, rows):
    """Write (shapely polygon, attribute values) rows to a GeoPackage layer of multipolygons."""
    layer = QgsVectorLayer(f"MultiPolygon?crs={crs}", name, "memory")
    layer.dataProvider().addAttributes([QgsField(field, QGIS_TYPES[kind]) for field, kind in fields])
    layer.updateFields()
    features = []
    for wkb, (_, values) in zip(shapely.to_wkb([g for g, _ in rows]), rows):
        geometry = QgsGeometry()
        geometry.fromWkb(wkb)
        geometry.convertToMultiType()
        f = QgsFeature(layer.fields())
        f.setGeometry(geometry)
        f.setAttributes(values)
        features.append(f)
    layer.dataProvider().addFeatures(features)
    options = QgsVectorFileWriter.SaveVectorOptions()
    options.driverName = "GPKG"
    options.layerName = name
    result = QgsVectorFileWriter.writeAsVectorFormatV3(layer, str(path), QgsCoordinateTransformContext(), options)
    assert result[0] == QgsVectorFileWriter.NoError, result


# Layer styles

def values_in(layer):
    """The distinct values of a one-band class raster."""
    return set(np.unique(read_tif(layer.source())).tolist())


def paletted(layer, colours, labels):
    """One colour and legend label per class value, for the classes present in the layer."""
    present = values_in(layer)
    classes = [QgsPalettedRasterRenderer.Class(k, QColor(colours[k]), labels[k]) for k in colours if k in present]
    layer.setRenderer(QgsPalettedRasterRenderer(layer.dataProvider(), 1, classes))


def classed(layer, classes, opacity=1.0, hide_below=None):
    """Draw a value layer in classes, given as (highest value of the class, colour, legend label), lowest first.

    Each pixel takes the colour of the first class whose highest value it does not exceed (a discrete
    colour ramp), so the layer shows only the colours of its legend, one legend entry per class.
    Values up to hide_below are treated as no data, as the maps leave them out. QGIS saves this range
    with its upper end included, so hide_below lies halfway between two values that the layer holds.
    """
    fn = QgsColorRampShader(classes[0][0], classes[-1][0])
    fn.setColorRampType(QgsColorRampShader.Discrete)
    fn.setColorRampItemList([QgsColorRampShader.ColorRampItem(v, QColor(c), label) for v, c, label in classes])
    legend = QgsColorRampLegendNodeSettings()
    legend.setUseContinuousLegend(False)
    fn.setLegendSettings(legend)
    shader = QgsRasterShader()
    shader.setRasterShaderFunction(fn)
    renderer = QgsSingleBandPseudoColorRenderer(layer.dataProvider(), 1, shader)
    renderer.setClassificationMin(classes[0][0])
    renderer.setClassificationMax(classes[-1][0])
    renderer.setOpacity(opacity)
    layer.setRenderer(renderer)
    if hide_below is not None:
        layer.dataProvider().setUserNoDataValue(1, [QgsRasterRange(-1e12, hide_below)])


def whole_classes(classes):
    """Classes of maps.py on whole values, ((lowest, highest), colour, label), in the form classed() takes.

    A class ends half a unit above its highest value, so that each value takes the class of the whole
    number it rounds to, as on the maps; an open-ended class (highest value inf) takes all values above.
    """
    return [(hi + 0.5 if np.isfinite(hi) else OPEN_END, colour, label) for (_, hi), colour, label in classes]


def ndvi_classes():
    """The NDVI change classes of maps.py (maps.ndvi_change_class) for the values stored in hundredths, lowest first.

    A stored value on a break takes the class the maps give it: 0.15 is little change and 0.3 is
    strongly greener. The stored values keep the class of the computed values (in_steps).
    """
    half = NDVI_STEP / 2
    tops = [OPEN_END, NDVI_STRONG - half, NDVI_CHANGE + half, -NDVI_CHANGE - half, -NDVI_STRONG + half]
    return [(top, colour, label) for top, (colour, label) in reversed(list(zip(tops, NDVI_CLASSES)))]


def share_classes():
    """The agricultural-share classes of maps.py for the stored whole percentages, lowest first.

    A class runs from its lowest share up to the next class's lowest share, so it ends half a per cent
    below that share.
    """
    tops = [round(lo * 100) - 0.5 for lo, _, _ in AG_SHARE_CLASSES[1:]] + [OPEN_END]
    return [(top, colour, label) for top, (_, colour, label) in zip(tops, AG_SHARE_CLASSES)]


def alpha_band(layer):
    """The band that QGIS makes of a file's internal mask, or -1 if there is none."""
    n = layer.bandCount()
    return n if layer.dataProvider().colorInterpretation(n) == Qgis.RasterColorInterpretation.AlphaBand else -1


def grey(layer, lo, hi):
    """Black at lo to white at hi, with masked pixels left transparent."""
    renderer = QgsSingleBandGrayRenderer(layer.dataProvider(), 1)
    renderer.setAlphaBand(alpha_band(layer))
    stretch = QgsContrastEnhancement(layer.dataProvider().dataType(1))
    stretch.setContrastEnhancementAlgorithm(QgsContrastEnhancement.StretchToMinimumMaximum)
    stretch.setMinimumValue(lo)
    stretch.setMaximumValue(hi)
    renderer.setContrastEnhancement(stretch)
    layer.setRenderer(renderer)


def outline(layer, colour, width, dashed=False):
    """Draw a polygon layer as outlines without fill (width in mm)."""
    line = layer.renderer().symbol().symbolLayer(0)
    line.setFillColor(QColor(0, 0, 0, 0))
    line.setStrokeColor(QColor(colour))
    line.setStrokeWidth(width)
    if dashed:
        line.setStrokeStyle(Qt.PenStyle.DashLine)


def by_class(layer, field, colours, labels, width=0.1):
    """Polygons filled and outlined (width in mm) in their class colour, for the classes present.

    field may also be an expression that gives each polygon its class.
    """
    present = set(layer.uniqueValues(layer.fields().indexOf(field))) if field in layer.fields().names() else None
    categories = []
    for k, colour in colours.items():
        if present is not None and k not in present:
            continue
        categories.append(QgsRendererCategory(k, filled(layer, colour, width), labels[k]))
    layer.setRenderer(QgsCategorizedSymbolRenderer(field, categories))


def filled(layer, colour, width):
    """A polygon symbol filled and outlined (width in mm) in one colour."""
    symbol = QgsSymbol.defaultSymbol(layer.geometryType())
    symbol.setColor(QColor(colour))
    symbol.symbolLayer(0).setStrokeColor(QColor(colour))
    symbol.symbolLayer(0).setStrokeWidth(width)
    return symbol


# The project

def raster(out, name, title):
    """A raster layer of the package."""
    layer = QgsRasterLayer(str(out / name), title)
    assert layer.isValid(), name
    return layer


def vector(out, name, title):
    """A vector layer of the package (a GeoJSON file, or the one layer of a GeoPackage)."""
    layer = QgsVectorLayer(str(out / name), title, "ogr")
    assert layer.isValid(), name
    return layer


def result_layers(site, out):
    """(key, layer, shown when the project opens) for the results group, top of the layer list first."""
    radar_note = " (not measurable by radar here)" if radar_cannot_see_water(site) else ""
    study = vector(out, "study_area.geojson", "Study area")
    outline(study, INK, 0.6, dashed=True)
    bodies = vector(out, "water_bodies.gpkg", "Water bodies by form (polygons with attributes)")
    by_class(bodies, "class_code", WATER_COLORS, WATER_LABELS)
    water = raster(out, "water_class_10m.tif", "Water bodies by form (10 m raster)")
    paletted(water, WATER_COLORS, WATER_LABELS)
    cardinal = raster(out, "water_cardinal_10m.tif", "Orientation of rectilinear water")
    paletted(cardinal, {1: ACCENT}, {1: "Rectilinear, within 7.5° of N-S or E-W"})
    ponds = raster(out, "pond_density_20m.tif",
                   "Ponds and small water bodies per km² (all under 10 ha), 1 km window")
    classed(ponds, whole_classes(DENSITY_CLASSES), opacity=0.75, hide_below=0.5)  # counts, rounded as on the map
    change = raster(out, "water_change_20m.tif", f"Water, dry versus wet season{radar_note}")
    paletted(change, CHANGE_COLORS, CHANGE_LABELS)
    months = raster(out, "s1_water_months_20m.tif", f"Months with open water in {config.YEAR}{radar_note}")
    classed(months, whole_classes(MONTH_CLASSES), opacity=0.85, hide_below=0.5)  # 0 months not drawn
    grain = raster(out, "field_edge_direction_100m.tif", "Dominant field-boundary direction (100 m)")
    classed(grain, [(7.5, GRAIN_COLORS[0], "Within 7.5° of cardinal"), (15, GRAIN_COLORS[1], "7.5° to 15°"),
                    (30, GRAIN_COLORS[2], "15° to 30°"), (OPEN_END, GRAIN_COLORS[3], "30° to 45° (diagonal)")],
            opacity=0.8)
    ndvi = raster(out, "ndvi_change_20m.tif", "NDVI, wet minus dry season")
    classed(ndvi, ndvi_classes(), opacity=0.85)
    cycles = raster(out, "crop_cycles_20m.tif", f"Crop cycles on agricultural land, {config.YEAR}")
    paletted(cycles, {**CYCLE_COLORS, TOO_FEW_MONTHS: "#99ffffff"},  # white at 60%, as on the map
             {**CYCLE_LABELS, TOO_FEW_MONTHS: "Too few clear months"})
    share = raster(out, "ag_share_1km_20m.tif", "Share of land under agriculture, 1 km (%)")
    classed(share, share_classes(), opacity=0.7, hide_below=round(AG_SHARE_CLASSES[0][0] * 100) - 0.5)
    fields = raster(out, "field_class_20m.tif", f"Agricultural field systems, {config.YEAR}")
    paletted(fields, FIELD_COLORS, FIELD_LABELS)
    return [("study", study, True), ("bodies", bodies, True), ("water", water, False),
            ("cardinal", cardinal, False), ("ponds", ponds, False), ("change", change, False),
            ("months", months, False), ("grain", grain, False), ("ndvi", ndvi, False), ("cycles", cycles, False),
            ("share", share, False), ("fields", fields, True)]


def input_layers(site, out):
    """(key, layer, shown when the project opens) for the open-data group, top of the layer list first."""
    if (out / "ftw_fields_2024.gpkg").exists():
        ftw_layer = vector(out, "ftw_fields_2024.gpkg", f"Fields of The World {config.YEAR}, field polygons")
        outline(ftw_layer, FTW_COLOUR, 0.2)
    else:
        ftw_layer = raster(out, "ftw_field_outlines_10m.tif", f"Fields of The World {config.YEAR}, field outlines")
        paletted(ftw_layer, {1: FTW_COLOUR}, {1: "Field boundary"})
    worldcover = raster(out, "worldcover_2021_10m.tif", "ESA WorldCover 2021")
    paletted(worldcover, {k: c for k, (c, _) in WORLDCOVER.items()}, {k: n for k, (_, n) in WORLDCOVER.items()})
    jrc = raster(out, "jrc_occurrence_20m.tif", "JRC Global Surface Water, occurrence 1984 to 2020 (%)")
    classed(jrc, whole_classes(JRC_CLASSES), hide_below=0.5)  # the classes of the web map; 0 is never water
    transitions = raster(out, "jrc_transitions_20m.tif", "JRC Global Surface Water, transitions 1984 to 2020")
    trajectory_groups(transitions)
    layers = [("ftw", ftw_layer, False), ("worldcover", worldcover, False), ("jrc", jrc, False),
              ("transitions", transitions, False)]
    for season, composite in (("dry", "median"), ("wet", "lowest monthly median")):
        vv = raster(out, f"s1_vv_{season}_20m.tif",
                    f"Sentinel-1 VV (dB), {season} season {composite}, {season_text(site, season)}")
        grey(vv, *VV_SHOWN_DB)
        layers.append((f"vv_{season}", vv, False))
    for season in ("wet", "dry"):
        tc = raster(out, f"s2_truecolour_{season}_20m.tif",
                    f"Sentinel-2 true colour, {season} season, {season_text(site, season)}")
        renderer = QgsMultiBandColorRenderer(tc.dataProvider(), 1, 2, 3)
        renderer.setAlphaBand(alpha_band(tc))  # cloud gaps and the area outside the study area stay transparent
        tc.setRenderer(renderer)
        layers.append((f"truecolour_{season}", tc, season == "dry"))
    hillshade = raster(out, "dem_20m.tif", "Copernicus DEM, hillshade")
    renderer = QgsHillshadeRenderer(hillshade.dataProvider(), 1, 315, 45)
    renderer.setZFactor(DEM_Z_FACTOR)
    hillshade.setRenderer(renderer)
    return layers + [("hillshade", hillshade, False)]


def trajectory_groups(layer):
    """JRC transition classes in the groups of the long-term water map (maps.TRAJECTORY), one legend entry each."""
    renderer = QgsPalettedRasterRenderer(layer.dataProvider(), 1, [])
    renderer.setMultiValueClasses([QgsPalettedRasterRenderer.MultiValueClass(list(TRAJECTORY_CODES[group]),
                                                                             QColor(colour), label)
                                   for group, colour, label in TRAJECTORY])
    layer.setRenderer(renderer)


def make_project(site, out):
    """Both layer groups, styled, and the map layouts in <site>.qgz with relative paths, opening on the study area."""
    project = QgsProject.instance()
    project.clear()
    crs = QgsCoordinateReferenceSystem(f"EPSG:{site.meta['epsg']}")
    project.setCrs(crs)
    project.setTitle(f"SRTSR0576 rapid mapping: {site.meta['name']}, {config.YEAR}")
    project.setFilePathStorage(Qgis.FilePathType.Relative)
    metadata = QgsProjectMetadata()
    metadata.setTitle(project.title())
    metadata.setAuthor("Muhammad Dita Farel")
    metadata.setAbstract("Rapid mapping of water bodies, agricultural field systems and seasonal change from open "
                         "satellite data. README.txt in this folder describes the layers.")
    project.setMetadata(metadata)
    root = project.layerTreeRoot()
    keyed = {}
    for title, layers in ((RESULTS, result_layers(site, out)), (INPUTS, input_layers(site, out))):
        group = root.addGroup(title)
        for key, layer, shown in layers:
            project.addMapLayer(layer, False)
            group.addLayer(layer).setItemVisibilityChecked(shown)
            keyed[key] = layer
    extent = study_extent(site)
    project.viewSettings().setDefaultViewExtent(QgsReferencedRectangle(extent, crs))
    project.viewSettings().setPresetFullExtent(QgsReferencedRectangle(extent, crs))
    make_layouts(site, project, keyed, out)
    assert project.write(str(out / f"{site.key}.qgz")), "project not written"


def study_extent(site):
    """The study area's bounding box in the site CRS, with a margin of 2% of its longer side."""
    xs, ys = zip(*aoi_utm(site.key))
    extent = QgsRectangle(min(xs), min(ys), max(xs), max(ys))
    extent.grow(0.02 * max(extent.width(), extent.height()))
    return extent


# Map layouts

LAYOUT_DIR = "layouts"          # folder of a site package that holds the exported PDFs
LAYOUT_DPI = 200                # resolution of the raster images in the PDFs
PDF_JPEG_QUALITY = 55           # JPEG quality of the photographic images in the PDFs (the background)
MASK_LEVEL_SHARE = 0.002        # an opacity held by this share of an image's pixels is one of its true levels
PHOTO_COLOURS = 256             # an image with at least this many colours in a sample is photographic
PAGE_MM = (420.0, 297.0)        # A3, long side first
MARGIN_MM = 10.0                # page margin
TITLE_BAND_MM = 19.0            # title and subtitle, above the map
GRID_MARGIN_MM = 7.0            # graticule ticks and degree labels, outside each side of the map frame
GAP_MM = 6.0                    # space between the map and the panel blocks, and between blocks
MIN_SIDE_MM = 70.0              # a strip beside the map narrower than this takes no legend or text
MIN_COLUMN_MM = 70.0            # narrowest column of text under a wide map
TAIL_MM = 70.0                  # width of the north arrow, scale bar, CRS and credit block under a wide map
TALL_ASPECT = 0.8               # a study area narrower than this (width over height) goes on a portrait page
GRATICULE_DEG = 0.1             # graticule interval, in degrees
LEGEND_TEXT_MM = (56.0, 68.0, 84.0)  # widths legend labels may be wrapped to, narrowest first
WRAP = "~"                      # line break in legend labels
WRAP_SHARE = 0.95               # share of a text box's width that its wrapped lines may take
TEXT_SLACK_MM = 4.0             # room under a text box for a line more, where another font sets the text longer
CRS_SLACK_MM = 1.5              # room under the CRS text, above the credit line
NBSP = "\u00a0"                 # joins a number to its unit, so that a line does not break between them
FONT = "Segoe UI"               # the font of the static maps, where it can be found
FONT_FILES = ("/mnt/c/Windows/Fonts/segoeui.ttf", "/mnt/c/Windows/Fonts/segoeuib.ttf",
              "/mnt/c/Windows/Fonts/seguisb.ttf")
# QGIS sets the layouts in another font on a computer without Segoe UI: on Linux, usually DejaVu Sans,
# which sets text about a seventh wider (Segoe UI sets it taller). Text is wrapped and boxed for both.
FALLBACK_FONTS = ("DejaVu Sans",)
MONTHS = ("January", "February", "March", "April", "May", "June", "July", "August", "September", "October",
          "November", "December")
TITLE_PT, SUBTITLE_PT, HEADING_PT, TEXT_PT, NOTE_PT, CREDIT_PT = 20, 11, 10, 9, 8, 7.5
NORTH_MM = (8.0, 14.0)          # north arrow, width and height
NORTH_ARROW = "arrows/NorthArrow_02.svg"  # one of the SVG symbols that come with QGIS
SCALE_BAR_MM = 60.0             # longest scale bar
SCALE_SEGMENTS_KM = (0.5, 1, 2, 2.5, 5, 10, 20)
CREDIT = "Muhammad Dita Farel, PhD applicant, SRTSR0576 (supervisor: Assoc Prof Scott Hawken), 2026"
# when the TanDEM-X radar data of the Copernicus DEM were acquired (Copernicus DEM product handbook)
DEM_ACQUIRED = "December 2010 to January 2015"
SENTINEL_CREDIT = f"Contains modified Copernicus Sentinel data {config.YEAR}."
# Named styles of the project's layers that the layouts draw through their map themes
STYLE_BACKGROUND = "layout background"
STYLE_FORM = "layout form classes"
STYLE_ONE_COLOUR = "layout dry-season water"
STYLE_REFERENCE = "layout reference water"
STYLE_ORIENTATION = "layout orientation"
# The background is the dry-season true colour at lower contrast, lighter and partly desaturated, close to
# 36% of the image over white. Fading the pixels themselves, rather than drawing the image part
# transparent, keeps the JPEG in the PDF small.
BACKGROUND_CONTRAST = -40
BACKGROUND_BRIGHTNESS = 80
BACKGROUND_SATURATION = -45
# Water bodies get an outline so that the smallest (0.1 ha, about 0.3 mm across on the page) show, wider
# where open water is sparse; as the enlarged small bodies of the static maps, this is display only.
WATER_OUTLINE_MM = 0.3
SPARSE_WATER_OUTLINE_MM = 0.5
# The three groups of the map of rectilinear water by orientation, in the colours of maps.cardinal_map
ORIENTATION_GROUPS = {1: (ACCENT, "Rectilinear, within 7.5° of N-S or E-W"),
                      2: ("#3a5f7d", "Rectilinear, other orientation"), 3: ("#9fb8cf", "Other water bodies")}

_font = []  # the layout font family, once found


@dataclass
class Sheet:
    """One map layout: what its map draws and what its panel says."""
    number: int
    name: str           # file name of its PDF, after the number
    title: str
    layers: list        # (layer key, named style or None) drawn over the faded background
    legend_title: str
    legend: list        # (layer key, named style or None, [(label in that style or None if it has one, label)])
    analysis: str
    data: list          # lines of the data-source box

    @property
    def theme(self):
        """Name of the layout and of its map theme."""
        return f"{self.number:02d} {self.title}"

    @property
    def pdf(self):
        return f"{self.number:02d}_{self.name}.pdf"


def pct_text(v):
    """A percentage as the deck writes it (km2_text: one decimal, two under 1)."""
    return f"{km2_text(float(v))}%"


def km2(v):
    """An area in km² as the deck writes it (km2_text), the number joined to its unit."""
    return f"{km2_text(float(v))}{NBSP}km²"


def count_text(n, one, many):
    """A count with its noun, singular for one."""
    return f"{n:,} {one if n == 1 else many}"


def comparative(key):
    """The site's rows of the comparative CSV as {(metric, unit): value}; empty when the CSV has not been made."""
    path = config.OUT / "dataset" / f"waterfields_comparative_{config.YEAR}.csv"
    if not path.exists():
        return {}
    with open(path, newline="") as fh:
        return {(r["metric"], r["unit"]): float(r["value"]) for r in csv.DictReader(fh)
                if r["site_key"] == key and r["unit"] != "date" and r["value"] != ""}  # "" is not available


def months_text(site, season):
    """The calendar months a season touches, for example 'February to April 2024' (the monthly radar series)."""
    months = site.season_months(season)
    first, last = (date(config.YEAR, m, 1) for m in (months[0], months[-1]))
    return f"{first:%B} to {last:%B} {config.YEAR}"


def data_sources(site):
    """The data-source lines that the layouts use, by key."""
    dry, wet, year = season_text(site, "dry"), season_text(site, "wet"), config.YEAR
    return {
        "s2_dry": f"Copernicus Sentinel-2 L2A, dry-season median composite, 10 m, {dry}.",
        "s2_seasons": f"Copernicus Sentinel-2 L2A, median composites, 20 m: dry season {dry}; wet season {wet}.",
        "s2_monthly": f"Copernicus Sentinel-2 L2A, monthly NDVI medians, 20 m, January to December {year}.",
        "s2_dry_ndvi": f"Copernicus Sentinel-2 L2A, dry-season median NDVI, 10 m, {dry}.",
        "s1_monthly": f"Copernicus Sentinel-1 RTC, VV and VH monthly medians, 20 m, January to December {year}.",
        "s1_vv_monthly": f"Copernicus Sentinel-1 RTC, VV monthly medians, 20 m, January to December {year}.",
        "s1_seasons": (f"Copernicus Sentinel-1 RTC, VV monthly medians, 20 m: dry season {months_text(site, 'dry')}; "
                       f"wet season {months_text(site, 'wet')}."),
        "worldcover": "ESA WorldCover 2021, 10 m (cropland, where the agricultural land starts).",
        "dem": f"Copernicus DEM GLO-30, 30 m, from TanDEM-X radar acquired {DEM_ACQUIRED} (slope mask).",
        "jrc": "JRC Global Surface Water v1.3, transitions between 1984 and 2020 (Pekel et al. 2016), 30 m.",
        "farmland": (f"Agricultural land as on map 2: Copernicus Sentinel-2 L2A and Sentinel-1 RTC, January to "
                     f"December {year}; ESA WorldCover 2021; Copernicus DEM GLO-30, {DEM_ACQUIRED}."),
        "background": f"Background: Sentinel-2 L2A dry-season true colour, {dry}, faded.",
    }


def radar_threshold(water_stats):
    """The Sentinel-1 water threshold of a site and how it was set, from water_stats.json."""
    db = water_stats["thresholds"]["s1_vv_db_edge_otsu"]
    if "Otsu" in water_stats["threshold_method"]["s1_vv"]:
        return f"{db:.1f} dB, set from the data by Otsu's method"
    return f"a fixed {db:.0f} dB, as the data gave no clear split"


def radar_blind(site, drawn):
    """Analysis text for a radar map at a site where the radar cannot see the water (maps.radar_cannot_see_water)."""
    monthly = read_json(site.products / "accuracy_stats.json")["s1_vv_over_optical_water_db"]["monthly_median"]
    db = [v for v in monthly if v is not None]
    months = "every month" if len(db) == 12 else "every month with data"
    return (f"Not measurable by radar at this site. Over the open water mapped from Sentinel-2, the Sentinel-1 VV "
            f"backscatter stays between {min(db):.1f} and {max(db):.1f} dB in {months} of {config.YEAR}, well above "
            f"the level of open water, so the radar cannot tell water from land here. {drawn}")


def raster_values(out, name):
    """The distinct values of a one-band raster of the package (no data left out)."""
    with rasterio.open(out / name) as src:
        arr = src.read(1, masked=True)
    return set(np.unique(arr.compressed()).tolist())


def map_sheets(site, out):
    """The eleven layouts of a site, in order, with their legends and texts, all numbers from the products."""
    w, f, c = (site.load_json(name) for name in ("water_stats", "field_stats", "change_stats"))
    table = comparative(site.key)
    src = data_sources(site)
    blind = radar_cannot_see_water(site)
    year = config.YEAR
    tail = [src["background"], SENTINEL_CREDIT]
    study = ("study", None, [(None, "Study area")])
    outline_mm = SPARSE_WATER_OUTLINE_MM if w["dry_season_water_km2"] < SPARSE_KM2 else WATER_OUTLINE_MM
    sheets = []

    counts, areas = w["count_by_class"], w["area_km2_by_class"]
    sheets.append(Sheet(
        1, "water_bodies", "Water bodies by form", [("study", None), ("bodies", STYLE_FORM)],
        "Form class: number of bodies, mapped area",
        [("bodies", STYLE_FORM, [(WATER_LABELS[k], f"{WATER_LABELS[k]}: "
                                                   f"{count_text(counts[CLASSES[k]], 'body', 'bodies')}, "
                                                   f"{km2(areas[CLASSES[k]])}")
                                 for k in (1, 6, 2, 3, 4, 5, 7, 8) if counts.get(CLASSES[k])]), study],
        f"Open water in the dry season, mapped from Sentinel-2 and grouped by the shape of each body with fixed rules. "
        f"{w['feature_count']:,} water bodies of 0.1 ha or more cover {km2(w['dry_season_water_km2'])}, "
        f"{pct_text(w['dry_season_water_share_pct'])} of the study area, and ponds and small bodies under 10 ha number "
        f"{w['ponds_per_100km2']:.1f} per 100{NBSP}km². Names such as moat and canal describe shapes and do not "
        f"date a feature. Each body has a {outline_mm}{NBSP}mm outline so that the smallest show; the areas are as "
        f"mapped.",
        [src["s2_dry"], src["dem"]] + tail))

    area, share = f["area_km2_by_class"], f["share_of_agriculture_pct"]
    largest = FIELD_LABELS[max(FIELD_LABELS, key=lambda k: area[FIELD_LABELS[k]])]
    shown, _ = field_classes_shown(site)  # as on the static map, classes under MIN_CLASS_KM2 stay off the legend
    left_out = [FIELD_LABELS[k] for k in (1, 2, 3, 4) if k not in shown]
    sheets.append(Sheet(
        2, "field_systems", "Agricultural field systems", [("study", None), ("bodies", STYLE_ONE_COLOUR),
                                                           ("fields", None)],
        "Field system: mapped area",
        [("fields", None, [(FIELD_LABELS[k], f"{FIELD_LABELS[k]}: {km2(area[FIELD_LABELS[k]])}") for k in shown]),
         ("bodies", STYLE_ONE_COLOUR, [(None, f"Dry-season open water: {km2(w['dry_season_water_km2'])}")]),
         study],
        f"Agricultural land, sorted with fixed rules by the seasonal behaviour of each 20 m cell: flooding, "
        f"dry-season greenness, crop cycles and the radar signature of paddy. It covers "
        f"{km2(f['agricultural_km2'])}, {pct_text(f['agricultural_share_pct'])} of the study area. The "
        f"largest class is {lower_first(largest)}, {km2(area[largest])} or "
        f"{pct_text(share[largest])} of the farmland. Dry-season open water is drawn for reference."
        + (f" Left out of the legend as too small to see, under {MIN_CLASS_KM2}{NBSP}km²: "
           f"{', '.join(f'{lower_first(n)} ({km2(area[n])})' for n in left_out)}." if left_out else ""),
        [src["s2_seasons"], src["s2_monthly"], src["s1_monthly"], src["worldcover"], src["dem"]] + tail))

    water_km2 = c["water_km2"]
    by_change = {1: "Permanent water (both seasons)", 2: "Seasonal inundation (wet season only)",
                 3: "Dry-season only water"}
    if blind:
        change_text = radar_blind(site, "The few cells drawn (radar water in any wet-season month: "
                                        f"{km2(c['radar_water_change_km2']['wet_season_maximum_water'])}) "
                                        "are not a measure of seasonal change.")
    else:
        lo, hi = w["sensitivity"]["seasonal_inundation_km2_range"]
        change_text = (f"Open water seen by the Sentinel-1 radar, with one VV threshold ({radar_threshold(w)}) for "
                       f"both seasons: the dry-season median against each month of the wet season. Water in the wet "
                       f"season only covers {km2(water_km2[by_change[2]])} ({km2_text(lo)} to "
                       f"{km2(hi)} with the threshold 1 dB lower or higher), and water in both seasons "
                       f"{km2(water_km2[by_change[1]])}.")
    sheets.append(Sheet(
        3, "water_seasons_radar", "Water, dry versus wet season (radar)", [("study", None), ("change", None)],
        "Radar water class, not a measure here" if blind else "Radar water class: mapped area",
        [("change", None, [(CHANGE_LABELS[k], CHANGE_LABELS[k] if blind
                            else f"{CHANGE_LABELS[k]}: {km2(water_km2[by_change[k]])}") for k in (1, 2, 3)]),
         study],
        change_text, [src["s1_seasons"], src["dem"]] + tail))

    sheets.append(Sheet(
        4, "ndvi_change", "NDVI change, wet minus dry season", [("study", None), ("ndvi", None)],
        "NDVI, wet minus dry season",
        [("ndvi", None, [(label, label) for _, label in NDVI_CLASSES]), study],
        f"Change in the greenness of vegetation between the seasons: the Sentinel-2 NDVI of the wet-season "
        f"composite minus that of the dry-season composite. {km2(c['ndvi_greening_km2'])} turned greener "
        f"by more than {NDVI_CHANGE} and {km2(c['ndvi_browning_km2'])} browner. On agricultural land the "
        f"median change is {c['ag_ndvi_change_median']:+.2f}, and {pct_text(c['ag_ndvi_greening_share_pct'])} of the "
        f"farmland turned greener. Cells without a clear view in both seasons are left blank; "
        f"{pct_text(c['ndvi_change_coverage_pct'])} of the study area had a clear view in both. Over open water the "
        f"colours do not describe plants.",
        [src["s2_seasons"]] + tail))

    present = raster_values(out, "s1_water_months_20m.tif")
    months_legend = []
    for (lo, hi), _, label in MONTH_CLASSES:
        if any(lo <= v <= hi for v in present):
            permanent = hi == 12 and not blind
            months_legend.append((label, f"{label}: {km2(w['s1_permanent_water_km2'])}" if permanent
                                  else label))
    if blind:
        months_text_ = radar_blind(site, "The month counts drawn (radar water in 1 to 9 months: "
                                         f"{km2(w['s1_intermittent_water_km2'])}) do not show where water "
                                         "lies or for how long.")
    else:
        months_text_ = (f"The number of months of {year} in which the Sentinel-1 monthly median VV lies below the "
                        f"water threshold ({radar_threshold(w)}). Water in 10 or more months covers "
                        f"{km2(w['s1_permanent_water_km2'])}, and water in 1 to 9 months "
                        f"{km2(w['s1_intermittent_water_km2'])}, a wider measure than the wet-season-only water "
                        f"of map 3.")
    sheets.append(Sheet(
        5, "water_months", f"Months with open water, {year}", [("study", None), ("months", None)],
        "Months with radar water, not a measure here" if blind else f"Months with radar water in {year}",
        [("months", None, months_legend), study], months_text_, [src["s1_vv_monthly"], src["dem"]] + tail))

    long_term = c["long_term"]["area_km2"]
    jrc_km2 = {"permanent 1984 and 2020": table.get(("Water permanent in 1984 and 2020 (JRC 1984-2020)", "km2")),
               "seasonal 1984 and 2020": table.get(("Water seasonal in 1984 and 2020 (JRC 1984-2020)", "km2")),
               **{group: long_term[group] for group in ("new since 1984", "lost since 1984", "ephemeral")}}
    sheets.append(Sheet(
        6, "water_since_1984", "Water since 1984 (JRC transitions)", [("study", None), ("transitions", None)],
        "Water history, 1984 to 2020: area",
        [("transitions", None, [(label, label.replace(": ", ", ") if jrc_km2[group] is None
                                 else f"{label.replace(': ', ', ')}: {km2(jrc_km2[group])}")
                                for group, _, label in TRAJECTORY]), study],
        f"Long-term history of surface water from the JRC Global Surface Water transition classes, which compare "
        f"the first and the last year of the 1984 to 2020 record. Water present in both years covers "
        f"{km2(long_term['water since 1984'])}. A further {km2(long_term['new since 1984'])} is new, "
        f"water in 2020 not recorded in 1984, which does not by itself mean newly built, and "
        f"{km2(long_term['lost since 1984'])} was lost."
        + (" JRC records little water at this site, so the classes are hard to see at this scale."
           if max(long_term.values()) < SPARSE_KM2 else ""),
        [src["jrc"]] + tail))

    large = table.get(("Large ponds (1 to 10 ha)", "count"))
    sheets.append(Sheet(
        7, "pond_density", "Ponds and small water bodies per km²",
        [("study", None), ("bodies", STYLE_REFERENCE), ("ponds", None)],
        "Ponds and small water bodies per km²",
        [("ponds", None, [(label, label) for _, _, label in DENSITY_CLASSES]),
         ("bodies", STYLE_REFERENCE, [(None, f"Mapped water body: "
                                             f"{count_text(w['feature_count'], 'body', 'bodies')}, "
                                             f"{km2(w['dry_season_water_km2'])}")]), study],
        f"Density of ponds and small water bodies under 10 ha (rectilinear ponds, irregular ponds and bodies too "
        f"small to classify), counted at their centroid in a 1 km moving window. The study area holds "
        f"{w['ponds_per_100km2']:.1f} per 100{NBSP}km²"
        + (f", and {int(large):,} ponds measure 1 to 10 ha" if large is not None else "")
        + ". Mapped water bodies are drawn in grey.",
        [src["s2_dry"], src["dem"]] + tail))

    cycles, rule = f["cropping_cycles_share_pct"], f["crop_cycle_rule"]
    twice = table.get(("Farmland cropped twice or more", "%"))
    sheets.append(Sheet(
        8, "crop_cycles", f"Crop cycles in {year}", [("study", None), ("cycles", None)],
        "Crop cycles: share of farmland",
        [("cycles", None, [(CYCLE_LABELS[k], f"{CYCLE_LABELS[k]}: {pct_text(cycles[str(k)])}")
                           for k in (1, 2, 3, 0)]
          + [("Too few clear months", f"Too few clear months (under {rule['min_clear_months']})")]), study],
        f"Crop cycles of {year} on agricultural land, counted from peaks in the monthly Sentinel-2 NDVI: a peak of "
        f"{rule['peak_ndvi_min']} or more that rises {rule['prominence_min']} above the lowest values around it. The "
        f"shares are of the farmland with {rule['min_clear_months']} or more clear months."
        + (f" {pct_text(twice)} of the farmland was cropped twice or more, probably a lower bound, since monthly "
           f"optical images miss some second crops." if twice is not None else ""),
        [src["s2_monthly"], src["farmland"]] + tail))

    top = AG_SHARE_CLASSES[-1][2]
    sheets.append(Sheet(
        9, "agriculture_concentration", "Concentration of agricultural land (1 km)",
        [("study", None), ("share", None)],
        "Share of agricultural land within 1 km",
        [("share", None, [(label, f"{label}: {km2(f['concentrated_ag_km2'])}" if label == top else label)
                          for _, _, label in AG_SHARE_CLASSES]), study],
        f"Share of agricultural land within a 1 km moving window, which shows where farming is concentrated. "
        f"Agricultural land covers {km2(f['agricultural_km2'])}, {pct_text(f['agricultural_share_pct'])} of "
        f"the study area, and {km2(f['concentrated_ag_km2'])} lies where {CONCENTRATED:.0%} or more of the "
        f"surrounding window is farmland. Shares under {AG_SHARE_CLASSES[0][0]:.0%} are not drawn.",
        [src["farmland"]] + tail))

    fo = f["field_orientation"]
    ci = fo["cardinal_share_ci95_pct"]
    test = ("It passes the rotation test, in which the image is turned by 30°." if passes_rotation_test(fo) else
            "It fails the rotation test, in which the image is turned by 30°, so it is not interpreted.")
    sheets.append(Sheet(
        10, "field_grain", "Field-edge orientation (field grain)", [("study", None), ("grain", None)],
        "Dominant field-edge direction: deviation from N-S or E-W",
        [("grain", None, [(label, label) for label in ("Within 7.5° of cardinal", "7.5° to 15°", "15° to 30°",
                                                        "30° to 45° (diagonal)")]), study],
        f"The dominant direction of field edges in 100 m cells, from the structure tensor of the 10 m dry-season "
        f"NDVI over agricultural land; cells whose edges are only weakly directional are not drawn. "
        f"{pct_text(fo['cardinal_share_pct'])} of the strong field edges lie within 7.5° of north-south or east-west "
        f"(95% interval {pct_text(ci[0])} to {pct_text(ci[1])}), against {pct_text(fo['uniform_expectation_pct'])} if "
        f"directions were random and {pct_text(fo['controls']['non_agricultural_edges_pct'])} for edges outside "
        f"farmland. {test}",
        [src["s2_dry_ndvi"], src["farmland"]] + tail))

    passed = (w.get("orientation_control") or {}).get("passed")
    test = {True: "The share passes the rotation test, in which the image is turned by 30°.",
            False: "The share fails the rotation test, in which the image is turned by 30°, so it is reported but "
                   "not interpreted.",
            None: "There are too few rectilinear bodies to run the rotation test, so the share is reported but not "
                  "interpreted."}[passed]
    # bodies in each group of ORIENTATION_GROUPS, the CASE of the orientation style (layout_styles)
    rectilinear = sum(counts[CLASSES[k]] for k in RECTILINEAR)
    cardinal = sum(w["cardinal_rectilinear_by_class"][CLASSES[k]] for k in RECTILINEAR)
    in_group = {1: cardinal, 2: rectilinear - cardinal, 3: w["feature_count"] - rectilinear}
    sheets.append(Sheet(
        11, "water_orientation", "Orientation of rectilinear water bodies",
        [("study", None), ("bodies", STYLE_ORIENTATION)], "Water bodies by orientation: number of bodies",
        [("bodies", STYLE_ORIENTATION, [(label, f"{label}: {count_text(in_group[g], 'body', 'bodies')}")
                                        for g, (_, label) in ORIENTATION_GROUPS.items() if in_group[g]]), study],
        f"Water bodies with straight, perpendicular sides (straight canals, enclosing moats, rectilinear reservoirs "
        f"and rectilinear ponds), split by whether their axes lie within 7.5° of north-south and east-west. "
        f"{pct_text(w['cardinal_share_of_rectilinear_pct'])} of them do, against "
        f"{pct_text(fo['uniform_expectation_pct'])} if directions were random. {test}",
        [src["s2_dry"], src["dem"]] + tail))
    for sheet in sheets:  # no line break between a number and its unit
        sheet.legend_title = keep_units(sheet.legend_title)
        sheet.analysis = keep_units(sheet.analysis)
        sheet.data = [keep_units(line) for line in sheet.data]
        sheet.legend = [(key, style, [(original, keep_units(text)) for original, text in items])
                        for key, style, items in sheet.legend]
    return sheets


def keep_units(text):
    """A text with a no-break space between each number and the unit after it (m, km, ha, dB, mm), and
    between a day and its month (1 February)."""
    return re.sub(rf"(\d) (m|km|ha|dB|mm|{'|'.join(MONTHS)})\b", rf"\1{NBSP}\2", text)


def add_style(layer, name, restyle):
    """Add a named style to a layer: its current style changed by restyle(layer). The layer keeps its current style."""
    manager = layer.styleManager()
    current = manager.currentStyle()
    manager.addStyleFromLayer(name)
    manager.setCurrentStyle(name)
    restyle(layer)
    manager.setCurrentStyle(current)


def layout_styles(site, layers):
    """The named styles that the layouts draw: the faded background and three ways of drawing the water bodies."""
    def fade(layer):
        layer.brightnessFilter().setContrast(BACKGROUND_CONTRAST)
        layer.brightnessFilter().setBrightness(BACKGROUND_BRIGHTNESS)
        layer.hueSaturationFilter().setSaturation(BACKGROUND_SATURATION)
    add_style(layers["truecolour_dry"], STYLE_BACKGROUND, fade)

    bodies = layers["bodies"]
    sparse = site.load_json("water_stats")["dry_season_water_km2"] < SPARSE_KM2
    width = SPARSE_WATER_OUTLINE_MM if sparse else WATER_OUTLINE_MM
    add_style(bodies, STYLE_FORM, lambda layer: by_class(layer, "class_code", WATER_COLORS, WATER_LABELS, width))
    add_style(bodies, STYLE_ONE_COLOUR,
              lambda layer: layer.setRenderer(QgsSingleSymbolRenderer(filled(layer, WATER_COLORS[1], width))))
    add_style(bodies, STYLE_REFERENCE,
              lambda layer: layer.setRenderer(QgsSingleSymbolRenderer(filled(layer, WATER_FOR_REFERENCE, width))))
    rectilinear = ", ".join(str(k) for k in RECTILINEAR)
    group = (f"CASE WHEN class_code IN ({rectilinear}) AND cardinal = 1 THEN 1 "
             f"WHEN class_code IN ({rectilinear}) THEN 2 ELSE 3 END")
    add_style(bodies, STYLE_ORIENTATION, lambda layer: by_class(
        layer, group, {k: colour for k, (colour, _) in ORIENTATION_GROUPS.items()},
        {k: label for k, (_, label) in ORIENTATION_GROUPS.items()}, width))


def add_theme(project, sheet, layers):
    """A map theme showing a sheet's layers, in their named styles, over the faded background."""
    records = []
    for key, style in sheet.layers + [("truecolour_dry", STYLE_BACKGROUND)]:
        record = QgsMapThemeCollection.MapThemeLayerRecord(layers[key])
        if style:
            record.usingCurrentStyle = True
            record.currentStyle = style
        records.append(record)
    theme = QgsMapThemeCollection.MapThemeRecord()
    theme.setLayerRecords(records)
    project.mapThemeCollection().insert(sheet.theme, theme)


def layout_font():
    """The font family of the layouts: Segoe UI, as on the static maps, where its files can be read."""
    if not _font:
        for path in FONT_FILES:
            if os.path.exists(path):
                QFontDatabase.addApplicationFont(path)
        _font.append(FONT if FONT in QFontDatabase().families() else QFont().defaultFamily())
    return _font[0]


def text_format(size, bold=False, colour=INK):
    """A text format of the layout font, in points."""
    font = QFont(layout_font())
    font.setBold(bold)
    fmt = QgsTextFormat()
    fmt.setFont(font)
    fmt.setSize(size)
    fmt.setSizeUnit(Qgis.RenderUnit.Points)
    fmt.setColor(QColor(colour))
    return fmt


class TextMeasure:
    """Widths and heights of text in mm, as the layout's text renderer sets it, for wrapping and placing text.

    The project names the layout font, which a computer without it replaces by another, so each text
    is measured in the layout font and in each of FALLBACK_FONTS that this computer has (families),
    and the largest measure is kept. Text wrapped here (legend labels) then fits its width in all of
    them, and a text item, which QGIS wraps itself in whatever font it finds, gets the height that the
    tallest of them needs.
    """

    _context = None
    _families = None

    @property
    def context(self):
        """A render context in dots per mm, made on first use, once QGIS has started."""
        if self._context is None:
            self._context = QgsRenderContext()
            self._context.setScaleFactor(LAYOUT_DPI / 25.4)
        return self._context

    @property
    def families(self):
        """The layout font, then the fallback fonts that can be found."""
        if self._families is None:
            known = set(QFontDatabase().families())
            self._families = [layout_font()] + [f for f in FALLBACK_FONTS if f in known and f != layout_font()]
        return self._families

    def formats(self, fmt):
        """A text format in each of the families."""
        result = []
        for family in self.families:
            font = QFont(fmt.font())
            font.setFamily(family)
            copy = QgsTextFormat(fmt)
            copy.setFont(font)
            result.append(copy)
        return result

    def _width(self, text, fmt):
        return QgsTextRenderer.textWidth(self.context, fmt, [text]) / self.context.scaleFactor()

    def _height(self, lines, fmt):
        return (QgsTextRenderer.textHeight(self.context, fmt, lines, Qgis.TextLayoutMode.Rectangle)
                / self.context.scaleFactor())

    @staticmethod
    def _wrap(text, width, measure):
        """The lines of a text broken at spaces to fit a width; each paragraph (line of `text`) starts a line.

        Lines are kept to WRAP_SHARE of the width, as the layout sets text a little wider than measured
        here. A no-break space (as between a number and its unit) does not break.
        """
        lines = []
        for paragraph in text.split("\n"):
            line = ""
            for word in paragraph.split(" "):
                trial = f"{line} {word}" if line else word
                if line and measure(trial) > width * WRAP_SHARE:
                    lines.append(line)
                    line = word
                else:
                    line = trial
            lines.append(line)
        return lines

    def wrap(self, text, fmt, width):
        """The lines of a text broken to fit a width in every family."""
        formats = self.formats(fmt)
        return self._wrap(text, width, lambda t: max(self._width(t, f) for f in formats))

    def height(self, text, fmt, width):
        """Height of a text as a text item of a width wraps it, in the family that needs the most height."""
        return max(self._height(self._wrap(text, width, lambda t, f=f: self._width(t, f)), f)
                   for f in self.formats(fmt))


TEXT = TextMeasure()


def panel_text(sheet, block):
    """(heading, body, body size in points) of the analysis or the data-source block of a sheet."""
    if block == "analysis":
        return "What the map shows", sheet.analysis, TEXT_PT
    return "Data sources", "\n".join(sheet.data), NOTE_PT


def text_block(heading, body, body_pt, width):
    """(heading height, body height, total height), in mm, of a headed block of text at a width."""
    head = TEXT.height(heading, text_format(HEADING_PT, bold=True), width) if heading else 0.0
    body_h = TEXT.height(body, text_format(body_pt), width)
    return head, body_h, body_h + (head + 1.5 if heading else 0.0)


def tail_lines(site):
    """(item name, text, format, room under it in mm) of the CRS and the credit, from the top down."""
    return (("CRS", crs_text(site), text_format(NOTE_PT), CRS_SLACK_MM),
            ("Credit", CREDIT, text_format(CREDIT_PT, colour=MUTED), 0.5))


def tail_height(site, width):
    """Height of the block with the north arrow, scale bar, CRS and credit, at a width."""
    return NORTH_MM[1] + 3 + sum(TEXT.height(text, fmt, width) + room
                                 for _, text, fmt, room in tail_lines(site)) + 2


def crs_text(site):
    """The coordinate reference system of the maps, in words."""
    epsg = site.meta["epsg"]
    zone = f"{epsg % 100}{'N' if epsg // 100 == 326 else 'S'}"
    return (f"Coordinate reference system: WGS{NBSP}84 / UTM zone{NBSP}{zone} (EPSG:{epsg}). "
            f"Graticule in degrees of WGS{NBSP}84.")


def page_plan(site, legends, sheets):
    """The page, the map frame and the panel blocks (mm), the same for all the layouts of a site.

    The map frame keeps the shape of the study area and is made as large as the page allows while
    the panel blocks keep their room: the legend, the analysis and the data sources beside the map
    where the strip there is wide enough, otherwise in columns under it; the north arrow, scale bar,
    CRS and credit under the map or at the foot of the strip beside it. Blocks are sized for the
    longest text among the site's layouts, and legends is a list of (label width, legend width,
    legend height) for the label widths the legends were measured at. A study area narrower than
    TALL_ASPECT goes on a portrait page. Returns the page size, the map frame, the blocks
    {name: QRectF} and the width the legend labels are wrapped to.
    """
    extent = study_extent(site)
    aspect = extent.width() / extent.height()
    page_w, page_h = (PAGE_MM[1], PAGE_MM[0]) if aspect < TALL_ASPECT else PAGE_MM
    m, g, a = MARGIN_MM, GAP_MM, GRID_MARGIN_MM
    top = m + TITLE_BAND_MM
    heights = {}

    def need(block, width):
        """Height (mm) that the analysis, the data sources or the tail takes at a width."""
        key = (block, round(width, 1))
        if key not in heights:
            heights[key] = (tail_height(site, width) if block == "tail" else
                            max(text_block(*panel_text(sheet, block), width)[2] for sheet in sheets))
        return heights[key]

    map_w = min(page_w - 2 * m - 2 * a, (page_h - top - m - 2 * a) * aspect)
    while map_w > 100:
        map_h = map_w / aspect
        right = m + 2 * a + map_w             # right edge of the map with its graticule labels
        bottom = top + 2 * a + map_h
        side = QRectF(right + g, top + a, page_w - m - right - g, page_h - m - top - a)  # from the frame's top
        foot = QRectF(m, bottom + g, page_w - 2 * m, page_h - m - bottom - g)
        plan = beside(side, foot, right, need, legends) or under(side, foot, need, legends)
        if plan:
            return (page_w, page_h), QRectF(m + a, top + a, map_w, map_h), *plan
        map_w -= 2
    raise RuntimeError(f"{site.key}: the layout panel does not fit on the page")


def beside(side, foot, map_right, need, legends):
    """Legend, analysis and data stacked in the strip beside the map, the tail directly under the map or at the
    foot of the strip; (blocks, legend label width) or None if they do not fit."""
    fitting = [legend for legend in legends if legend[1] <= side.width()]
    if side.width() < MIN_SIDE_MM or not fitting:
        return None
    text_mm, legend_w, legend_h = fitting[-1]
    blocks, y, g = {}, side.y(), GAP_MM
    for name in ("legend", "analysis", "data"):
        h = legend_h if name == "legend" else need(name, side.width())
        blocks[name] = QRectF(side.x(), y, side.width(), h)
        y += h + g
    under_map = QRectF(foot.x(), foot.y(), min(map_right - foot.x(), 2 * TAIL_MM), foot.height())
    if under_map.height() >= need("tail", under_map.width()):  # directly under the map
        h = need("tail", under_map.width())
        blocks["tail"] = QRectF(under_map.x(), under_map.y(), under_map.width(), h)
        fits = y - g <= side.bottom()
    else:
        h = need("tail", side.width())
        blocks["tail"] = QRectF(side.x(), side.bottom() - h, side.width(), h)
        fits = y <= side.bottom() - h
    # the labels may take the whole width of the strip, which only makes the legend shorter
    return (blocks, text_mm + side.width() - legend_w) if fits else None


def under(side, foot, need, legends):
    """Legend, analysis and data in columns under the map, and the tail beside the map at the foot of the strip
    there if it fits, else in a last column; (blocks, legend label width) or None if they do not fit.

    Of the legend widths, the one that leaves the shortest columns is used.
    """
    g, blocks = GAP_MM, {}
    columns = ["legend", "analysis", "data"]
    tail_h = need("tail", side.width()) if side.width() >= TAIL_MM else None
    if tail_h is not None and side.height() - foot.height() - g >= tail_h:
        blocks["tail"] = QRectF(side.x(), foot.y() - g - tail_h, side.width(), tail_h)
    else:
        columns.append("tail")
    best = None
    for text_mm, legend_w, legend_h in legends:
        text_w = (foot.width() - legend_w - (TAIL_MM if "tail" in columns else 0) - g * (len(columns) - 1)) / 2
        if text_w < MIN_COLUMN_MM:
            continue
        widths = {"legend": legend_w, "analysis": text_w, "data": text_w, "tail": TAIL_MM}
        heights = {name: legend_h if name == "legend" else need(name, widths[name]) for name in columns}
        tallest = max(heights.values())
        if tallest <= foot.height() and (best is None or tallest < best[0]):
            best = (tallest, text_mm, widths, heights)
    if best is None:
        return None
    _, text_mm, widths, heights = best
    x = foot.x()
    for name in columns:
        y = foot.bottom() - heights[name] if name == "tail" else foot.y()
        blocks[name] = QRectF(x, y, widths[name], heights[name])
        x += widths[name] + g
    return blocks, text_mm


def legend_labels(layer, style):
    """The labels of a layer's legend entries in one of its named styles (None: its current style), in order."""
    manager = layer.styleManager()
    current = manager.currentStyle()
    if style:
        manager.setCurrentStyle(style)
    nodes = layer.legend().createLayerTreeModelLegendNodes(QgsLayerTreeLayer(layer))
    labels = [node.data(Qt.ItemDataRole.DisplayRole) for node in nodes]
    manager.setCurrentStyle(current)
    return labels


def legend_label(text, width):
    """A legend label wrapped to a width (mm) with the legend's line-break character."""
    return WRAP.join(TEXT.wrap(text, text_format(TEXT_PT), width))


def add_legend(layout, sheet, layers, text_width):
    """The sheet's legend: chosen entries of its layers with the sheet's labels, in the sheet's order.

    Labels are wrapped to text_width (mm).
    """
    legend = QgsLayoutItemLegend(layout)
    layout.addLayoutItem(legend)
    legend.setId("Legend")
    legend.setAutoUpdateModel(False)
    legend.setTitle(WRAP.join(TEXT.wrap(sheet.legend_title, text_format(HEADING_PT, bold=True), text_width)))
    legend.setWrapString(WRAP)
    legend.rstyle(QgsLegendStyle.Title).setTextFormat(text_format(HEADING_PT, bold=True))
    legend.rstyle(QgsLegendStyle.Title).setMargin(QgsLegendStyle.Bottom, 2.0)
    for style in (QgsLegendStyle.Group, QgsLegendStyle.Subgroup, QgsLegendStyle.SymbolLabel):
        legend.rstyle(style).setTextFormat(text_format(TEXT_PT))
    legend.rstyle(QgsLegendStyle.Subgroup).setMargin(QgsLegendStyle.Top, 1.6)
    legend.rstyle(QgsLegendStyle.Symbol).setMargin(QgsLegendStyle.Top, 1.6)
    legend.setSymbolWidth(6)
    legend.setSymbolHeight(4)
    legend.setDrawRasterStroke(True)  # a thin grey edge, so that pale swatches show on white
    legend.setRasterStrokeColor(QColor(SWATCH_OUTLINE))
    legend.setRasterStrokeWidth(0.15)
    legend.setBoxSpace(0)
    legend.setBackgroundEnabled(False)
    root = legend.model().rootGroup()
    root.clear()
    for key, style, items in sheet.legend:
        layer = layers[key]
        labels = legend_labels(layer, style)
        node = root.addLayer(layer)
        order, texts = [], []
        for original, text in items:
            index = 0 if original is None else labels.index(original) if original in labels else None
            if index is None:  # a class that this site does not have
                continue
            order.append(index)
            texts.append(text)
            QgsMapLayerLegendUtils.setLegendNodeUserLabel(node, index, legend_label(text, text_width))
        QgsMapLayerLegendUtils.setLegendNodeOrder(node, order)
        if isinstance(layer, QgsVectorLayer) and len(order) == 1:
            node.setCustomProperty("legend/title-label", legend_label(texts[0], text_width))  # in the layer's row
        else:
            QgsLegendRenderer.setNodeLegendStyle(node, QgsLegendStyle.Hidden)
        legend.model().refreshLayerLegend(node)
    return legend


LEGEND_TEXT_STYLES = (QgsLegendStyle.Title, QgsLegendStyle.Group, QgsLegendStyle.Subgroup, QgsLegendStyle.SymbolLabel)


def legend_size(layout, legend):
    """Width and height (mm) that a legend takes in the family that makes it largest (TextMeasure.families)."""
    kept = {style: QgsTextFormat(legend.rstyle(style).textFormat()) for style in LEGEND_TEXT_STYLES}
    width = height = 0.0
    for i in range(len(TEXT.families)):
        for style in LEGEND_TEXT_STYLES:
            legend.rstyle(style).setTextFormat(TEXT.formats(kept[style])[i])
        size = QgsLegendRenderer(legend.model(), legend.legendSettings()).minimumSize(
            QgsLayoutUtils.createRenderContextForLayout(layout, None))
        width, height = max(width, size.width()), max(height, size.height())
    for style in LEGEND_TEXT_STYLES:
        legend.rstyle(style).setTextFormat(kept[style])
    return width, height


def add_label(layout, name, text, fmt, rect):
    """A text item without margins."""
    label = QgsLayoutItemLabel(layout)
    layout.addLayoutItem(label)
    label.setId(name)
    label.setText(text)
    label.setTextFormat(fmt)
    label.setMarginX(0)
    label.setMarginY(0)
    label.attemptSetSceneRect(rect)
    return label


def add_text_block(layout, name, heading, body, body_pt, rect):
    """A heading and a text under it, from the top of rect.

    The text items get the text unbroken, so that QGIS wraps it to their width in whatever font a
    computer has; the boxes are sized from the wrap measured here in the layout font and the fallback
    fonts (TextMeasure), which allows a little more room (WRAP_SHARE), with TEXT_SLACK_MM to spare
    under each.
    """
    head_h, body_h, _ = text_block(heading, body, body_pt, rect.width())
    y = rect.y()
    if heading:
        add_label(layout, f"{name} heading", heading, text_format(HEADING_PT, bold=True),
                  QRectF(rect.x(), y, rect.width(), head_h + 0.5))
        y += head_h + 1.5
    add_label(layout, name, body, text_format(body_pt), QRectF(rect.x(), y, rect.width(), body_h + TEXT_SLACK_MM))


def add_map(layout, site, sheet, rect):
    """The map frame, drawing the sheet's theme over the study area, with a graticule."""
    item = QgsLayoutItemMap(layout)
    layout.addLayoutItem(item)
    item.setId("Map")
    item.attemptSetSceneRect(rect)
    item.setCrs(QgsCoordinateReferenceSystem(f"EPSG:{site.meta['epsg']}"))
    item.zoomToExtent(study_extent(site))
    item.setFollowVisibilityPreset(True)
    item.setFollowVisibilityPresetName(sheet.theme)
    item.setFrameEnabled(True)
    item.setFrameStrokeColor(QColor(INK))
    item.setFrameStrokeWidth(QgsLayoutMeasurement(0.25))
    grid = QgsLayoutItemMapGrid("Graticule (WGS 84 degrees)", item)
    grid.setCrs(QgsCoordinateReferenceSystem("EPSG:4326"))
    grid.setIntervalX(GRATICULE_DEG)
    grid.setIntervalY(GRATICULE_DEG)
    grid.setStyle(QgsLayoutItemMapGrid.Solid)
    grid.setLineSymbol(QgsLineSymbol.createSimple({"color": "29,29,31,80", "width": "0.12"}))
    grid.setFrameStyle(QgsLayoutItemMapGrid.ExteriorTicks)
    grid.setFrameWidth(1.6)
    grid.setFramePenSize(0.25)
    grid.setFramePenColor(QColor(INK))
    grid.setAnnotationEnabled(True)
    grid.setAnnotationFormat(QgsLayoutItemMapGrid.DecimalWithSuffix)
    grid.setAnnotationPrecision(1)
    grid.setAnnotationTextFormat(text_format(NOTE_PT, colour=MUTED))
    grid.setAnnotationFrameDistance(0.8)
    for side, shown, direction in ((QgsLayoutItemMapGrid.Left, QgsLayoutItemMapGrid.LatitudeOnly,
                                    QgsLayoutItemMapGrid.Vertical),
                                   (QgsLayoutItemMapGrid.Right, QgsLayoutItemMapGrid.LatitudeOnly,
                                    QgsLayoutItemMapGrid.Vertical),
                                   (QgsLayoutItemMapGrid.Top, QgsLayoutItemMapGrid.LongitudeOnly,
                                    QgsLayoutItemMapGrid.Horizontal),
                                   (QgsLayoutItemMapGrid.Bottom, QgsLayoutItemMapGrid.LongitudeOnly,
                                    QgsLayoutItemMapGrid.Horizontal)):
        grid.setAnnotationPosition(QgsLayoutItemMapGrid.OutsideMapFrame, side)
        grid.setAnnotationDisplay(shown, side)
        grid.setAnnotationDirection(direction, side)
    item.grids().addGrid(grid)
    return item


def add_tail(layout, site, item, rect):
    """North arrow and scale bar side by side, then the CRS and the credit line, from the top of rect."""
    x, y, width = rect.x(), rect.y(), rect.width()
    arrow = QgsLayoutItemPicture(layout)
    layout.addLayoutItem(arrow)
    arrow.setId("North arrow")
    arrow.setPicturePath(QgsSymbolLayerUtils.svgSymbolNameToPath(NORTH_ARROW, QgsProject.instance().pathResolver()),
                         Qgis.PictureFormat.SVG)
    arrow.setSvgFillColor(QColor(INK))
    arrow.setSvgStrokeColor(QColor(INK))
    arrow.setLinkedMap(item)
    arrow.setNorthMode(QgsLayoutItemPicture.TrueNorth)
    arrow.attemptSetSceneRect(QRectF(x, y, *NORTH_MM))

    bar = QgsLayoutItemScaleBar(layout)
    layout.addLayoutItem(bar)
    bar.setId("Scale bar")
    bar.setLinkedMap(item)
    bar.setStyle("Single Box")
    bar.setUnits(Qgis.DistanceUnit.Kilometers)
    bar.setUnitLabel("km")
    bar.setSegmentSizeMode(Qgis.ScaleBarSegmentSizeMode.Fixed)
    mm_per_km = item.sizeWithUnits().width() / (item.extent().width() / 1000)
    longest = min(SCALE_BAR_MM, width - NORTH_MM[0] - 12)
    bar.setUnitsPerSegment(max([s for s in SCALE_SEGMENTS_KM if 2 * s * mm_per_km <= longest] or
                               [SCALE_SEGMENTS_KM[0]]))
    bar.setNumberOfSegments(2)
    bar.setNumberOfSegmentsLeft(0)
    bar.setHeight(1.8)
    bar.setLabelBarSpace(1.2)
    bar.setBoxContentSpace(0)
    bar.setTextFormat(text_format(NOTE_PT))
    bar.attemptMove(QgsLayoutPoint(x + NORTH_MM[0] + 6, y + NORTH_MM[1] - 7))
    bar.resizeToMinimumWidth()

    y += NORTH_MM[1] + 3
    for name, text, fmt, room in tail_lines(site):  # QGIS wraps each text to the width, in the font it has
        h = TEXT.height(text, fmt, width) + room
        add_label(layout, name, text, fmt, QRectF(x, y, width, h))
        y += h + 2


def fill_layout(site, sheet, layout, legend, page, map_rect, blocks, count):
    """Page size, title, map, legend, texts, north arrow, scale bar, CRS and credit of one layout."""
    layout.pageCollection().page(0).setPageSize(QgsLayoutSize(*page))
    layout.renderContext().setDpi(LAYOUT_DPI)
    layout.setCustomProperty("waterfields/pdf", sheet.pdf)
    # the PDF and SVG export settings that the layout designer offers: text as text, rasters in one piece
    layout.setCustomProperty("pdfTextFormat", int(Qgis.TextRenderFormat.AlwaysText))
    layout.setCustomProperty("svgTextFormat", int(Qgis.TextRenderFormat.AlwaysText))
    layout.setCustomProperty("pdfDisableRasterTiles", True)
    m, width = MARGIN_MM, page[0] - 2 * MARGIN_MM
    add_label(layout, "Title", sheet.title, text_format(TITLE_PT, bold=True), QRectF(m, m - 1, width, 9))
    add_label(layout, "Subtitle", f"{site.meta['name']} study area, {config.YEAR}. SRTSR0576 rapid mapping, "
                                  f"map {sheet.number} of {count}.",
              text_format(SUBTITLE_PT, colour=MUTED), QRectF(m, m + 9, width, 6))
    item = add_map(layout, site, sheet, map_rect)
    legend.setLinkedMap(item)  # the legend draws the layers in the styles of the map's theme
    legend.attemptMove(QgsLayoutPoint(blocks["legend"].x(), blocks["legend"].y()))
    legend.setResizeToContents(True)
    legend.adjustBoxSize()
    for name in ("analysis", "data"):
        add_text_block(layout, "Analysis" if name == "analysis" else "Data sources", *panel_text(sheet, name),
                       blocks[name])
    add_tail(layout, site, item, blocks["tail"])


def make_layouts(site, project, layers, out):
    """One print layout per sheet, each drawing its own map theme, laid out alike for all the sheets of a site."""
    layout_styles(site, layers)
    sheets = map_sheets(site, out)
    built = []
    for sheet in sheets:
        add_theme(project, sheet, layers)
        layout = QgsPrintLayout(project)
        layout.initializeDefaults()
        layout.setName(sheet.theme)
        project.layoutManager().addLayout(layout)
        built.append((sheet, layout))
    legends = []  # the largest legend of the site at each label width
    for text_mm in LEGEND_TEXT_MM:
        sizes = []
        for sheet, layout in built:
            legend = add_legend(layout, sheet, layers, text_mm)
            sizes.append(legend_size(layout, legend))
            layout.removeLayoutItem(legend)
        legends.append((text_mm, max(w for w, _ in sizes), max(h for _, h in sizes)))
    page, map_rect, blocks, text_mm = page_plan(site, legends, sheets)
    # a label width that page_plan widened may still set a legend wider or taller than its block in one
    # of the fonts; it then comes down until every legend fits (at the width measured above, they do)
    box, measured = blocks["legend"], max(t for t in LEGEND_TEXT_MM if t <= text_mm + 0.01)
    while True:
        made = [(sheet, layout, add_legend(layout, sheet, layers, text_mm)) for sheet, layout in built]
        sizes = [legend_size(layout, legend) for _, layout, legend in made]
        if all(w <= box.width() + 0.01 and h <= box.height() + 0.01 for w, h in sizes) or text_mm <= measured:
            break
        for _, layout, legend in made:
            layout.removeLayoutItem(legend)
        text_mm = max(measured, text_mm - 1)
    for sheet, layout, legend in made:
        fill_layout(site, sheet, layout, legend, page, map_rect, blocks, len(sheets))
    if hasattr(project, "labelingEngineSettings"):  # text exported as text, also when exported from QGIS
        settings = project.labelingEngineSettings()
        settings.setDefaultTextRenderFormat(Qgis.TextRenderFormat.AlwaysText)
        project.setLabelingEngineSettings(settings)


def export_layouts(key, out):
    """Each layout of the open project to layouts/NN_<name>.pdf, and all of them to <site>_map_atlas.pdf.

    The images go in at LAYOUT_DPI: the class maps lossless and the background as JPEG where
    PyMuPDF is importable (see recompress), all as JPEG (QGIS's lossy default) where it is not.
    Returns the paths of the PDFs written.
    """
    folder = out / LAYOUT_DIR
    folder.mkdir(exist_ok=True)
    layouts = sorted(QgsProject.instance().layoutManager().printLayouts(), key=lambda layout: layout.name())
    settings = QgsLayoutExporter.PdfExportSettings()
    settings.dpi = LAYOUT_DPI
    settings.rasterizeWholeImage = False
    settings.forceVectorOutput = False
    settings.appendGeoreference = True
    settings.exportMetadata = True
    settings.simplifyGeometries = True
    settings.textRenderFormat = Qgis.TextRenderFormat.AlwaysText
    lossless = pdf_tools() is not None  # lossless images, of which recompress turns the background into JPEG
    paths = []
    for layout in layouts:
        context = layout.renderContext()
        # each raster drawn as one image, so that no seams show between tiles
        context.setFlag(QgsLayoutRenderContext.FlagDisableTiledRasterLayerRenders, True)
        context.setFlag(QgsLayoutRenderContext.FlagLosslessImageRendering, lossless)
        settings.flags = context.flags()
        path = folder / layout.customProperty("waterfields/pdf")
        result = QgsLayoutExporter(layout).exportToPdf(str(path), settings)
        context.setFlag(QgsLayoutRenderContext.FlagLosslessImageRendering, False)  # QGIS's default for later exports
        assert result == QgsLayoutExporter.Success, (path, result)
        if lossless:
            recompress(path)
        paths.append(path)
    atlas = join_pdfs(paths, [layout.name() for layout in layouts], folder / f"{key}_map_atlas.pdf")
    return paths + ([atlas] if atlas else [])


def pdf_tools():
    """PyMuPDF (fitz) when it can be imported, else None."""
    try:
        import fitz
    except ImportError:
        return None
    return fitz


def recompress(path):
    """Tidy the images of an exported PDF with PyMuPDF: transparent pixels coloured, the background as JPEG.

    The layouts are exported with lossless images where PyMuPDF is at hand: class maps of a few flat
    colours stay exact and compress well that way. An image with PHOTO_COLOURS colours or more in a
    sample of every fourth pixel in both directions, which is the faded background, is stored as JPEG
    at PDF_JPEG_QUALITY instead, many times smaller; one of 256 colours or fewer gets a palette. Before
    that, each soft mask is snapped to its true opacities (see snap_opacity), and becomes a 1-bit
    stencil where it has no other opacity than full and none; the pixels a mask hides are coloured
    (see bleed).
    """
    fitz = pdf_tools()
    doc = fitz.open(path)
    seen = set()
    for page in doc:
        for xref, mask, *_ in page.get_images(full=True):
            if xref in seen:
                continue
            seen.add(xref)
            pix = fitz.Pixmap(doc, xref)
            if pix.n != 3 or pix.alpha:
                continue
            rgb = np.frombuffer(pix.samples, np.uint8).reshape(pix.height, pix.width, 3).copy()
            photo = len(np.unique(packed(rgb[::4, ::4]))) >= PHOTO_COLOURS
            if mask:
                alpha = fitz.Pixmap(doc, mask)
                if (alpha.width, alpha.height, alpha.n) == (pix.width, pix.height, 1):
                    opacity = np.frombuffer(alpha.samples, np.uint8)
                    soft = (doc.xref_get_key(mask, "ImageMask")[1] != "true"
                            and doc.xref_get_key(mask, "BitsPerComponent")[1] == "8")
                    if soft and mask not in seen:  # an 8-bit soft mask, not a 1-bit stencil
                        opacity = snap_opacity(opacity)
                        if set(np.unique(opacity).tolist()) <= {0, 255}:
                            stencil(doc, xref, mask, opacity.reshape(pix.height, pix.width))
                        else:
                            doc.update_stream(mask, opacity.tobytes(), compress=True)
                            doc.xref_set_key(mask, "Filter", "/FlateDecode")
                            doc.xref_set_key(mask, "DecodeParms", "null")
                    seen.add(mask)
                    bleed(rgb, opacity.reshape(pix.height, pix.width) == 0)
            colours, index = np.unique(packed(rgb).ravel(), return_inverse=True)
            if photo:
                image = fitz.Pixmap(fitz.csRGB, pix.width, pix.height, rgb.tobytes(), False)
                doc.update_stream(xref, image.tobytes("jpeg", jpg_quality=PDF_JPEG_QUALITY), compress=False)
                doc.xref_set_key(xref, "Filter", "/DCTDecode")
            elif len(colours) <= 256:  # a palette and one byte per pixel, a third of the size before compression
                palette = np.stack([colours >> 16, colours >> 8, colours], axis=1).astype(np.uint8)
                doc.update_stream(xref, index.astype(np.uint8).tobytes(), compress=True)
                doc.xref_set_key(xref, "ColorSpace",
                                 f"[/Indexed /DeviceRGB {len(colours) - 1} <{palette.tobytes().hex()}>]")
            else:
                doc.update_stream(xref, rgb.tobytes(), compress=True)
            doc.xref_set_key(xref, "DecodeParms", "null")
    tmp = path.with_suffix(".tmp")
    doc.save(tmp, garbage=4, deflate=True)
    doc.close()
    tmp.replace(path)


def stencil(doc, image, mask, opacity):
    """Replace the soft mask of an image, wholly opaque or wholly clear at each pixel, by a 1-bit stencil mask.

    The stencil is written as Qt writes its own: one bit per pixel, 1 where the image shows (Decode
    [1 0]), each row padded to a whole byte. It takes an eighth of the soft mask before compression.
    """
    doc.update_stream(mask, np.packbits(opacity > 0, axis=1).tobytes(), compress=True)
    for key, value in (("ColorSpace", "null"), ("BitsPerComponent", "null"), ("ImageMask", "true"),
                       ("Decode", "[1 0]"), ("Filter", "/FlateDecode"), ("DecodeParms", "null")):
        doc.xref_set_key(mask, key, value)
    doc.xref_set_key(image, "SMask", "null")
    doc.xref_set_key(image, "Mask", f"{mask} 0 R")


def snap_opacity(alpha):
    """An image mask with each value moved to the nearest of its true levels (0 and the opacities of its layers).

    Qt smooths the edges of a scaled image, which leaves every value from 0 to 255 along them; snapped
    to its few true levels, a mask compresses to about half. The edges move by less than a pixel.
    """
    values, counts = np.unique(alpha, return_counts=True)
    levels = np.union1d(values[counts >= MASK_LEVEL_SHARE * alpha.size], [0]).astype(np.int16)
    return levels[np.abs(alpha[:, None].astype(np.int16) - levels[None, :]).argmin(axis=1)].astype(np.uint8)


def packed(rgb):
    """The colours of an RGB image as single integers, 0xRRGGBB."""
    rgb = rgb.astype(np.uint32)
    return (rgb[..., 0] << 16) | (rgb[..., 1] << 8) | rgb[..., 2]


def bleed(rgb, hidden, reach=2):
    """Colour the hidden pixels of an image (in place): the nearest visible colour within reach pixels, white beyond.

    Qt leaves them black. A viewer that smooths an image mixes them into the edge of every patch
    beside them, which then shows a dark rim; coloured like their neighbours they leave no trace.
    """
    if not hidden.any() or hidden.all():
        return
    dist, (iy, ix) = ndimage.distance_transform_edt(hidden, return_indices=True)
    near = hidden & (dist <= reach)
    rgb[near] = rgb[iy[near], ix[near]]
    rgb[hidden & ~near] = 255


def join_pdfs(paths, titles, atlas):
    """One PDF of all the layouts, with a bookmark for each, by PyMuPDF, pypdf or PyPDF2; None when none is
    importable."""
    fitz = pdf_tools()
    if fitz:
        doc = fitz.open()
        for path in paths:
            with fitz.open(path) as part:
                doc.insert_pdf(part)
        doc.set_toc([[1, title, page] for page, title in enumerate(titles, start=1)])
        doc.save(atlas, garbage=4, deflate=True)  # garbage=4 keeps one copy of identical images
        doc.close()
        return atlas
    try:
        from pypdf import PdfWriter as Merger
    except ImportError:
        try:
            from PyPDF2 import PdfMerger as Merger
        except ImportError:
            return None
    merger = Merger()
    for path, title in zip(paths, titles):
        merger.append(str(path), title)  # the second argument is the bookmark in pypdf and in every PyPDF2
    if hasattr(merger, "compress_identical_objects"):  # pypdf only
        merger.compress_identical_objects()
    merger.write(str(atlas))
    merger.close()
    return atlas


# README.txt and the zip

def readme(site, out, ftw_note):
    """The README.txt of a site folder."""
    name, epsg, year = site.meta["name"], site.meta["epsg"], config.YEAR
    qgis_version = Qgis.version().split("-")[0]
    zone = f"{epsg % 100}{'N' if epsg // 100 == 326 else 'S'}"
    boundary = site.meta["aoi_source"]
    if boundary.startswith("SRT AOI.kml: "):
        boundary = f'the polygon "{boundary.removeprefix("SRT AOI.kml: ")}" in SRT AOI.kml'
    else:
        boundary = f"a {boundary}"  # the Angkor square
    notes = ftw_note
    if radar_cannot_see_water(site):
        notes = ("\n- At this site the Sentinel-1 radar cannot separate open water from land: backscatter over the\n"
                 "  water mapped from Sentinel-2 stays close to -11 dB in every month. The two radar water layers\n"
                 "  (dry versus wet season, and months with open water) are therefore not a measure here, as\n"
                 "  layouts 3 and 5 say." + ftw_note)
    extent = study_extent(site)
    page = "portrait" if extent.width() / extent.height() < TALL_ASPECT else "landscape"
    sheets = map_sheets(site, out)
    layouts = "\n".join(f"   {LAYOUT_DIR}/{sheet.pdf:<35} {sheet.title}" for sheet in sheets)
    return f"""SRTSR0576 rapid mapping of fields and water bodies: QGIS project for {name}, {year}

This folder holds a QGIS project of the rapid mapping of water bodies, agricultural field systems and
seasonal change, shown together with the open data that the mapping used. To use it, keep the files
together in this folder and open {site.key}.qgz in QGIS. The project was saved with QGIS {qgis_version}. It
refers to its layers by relative paths, so the folder can be copied or moved as a whole.

Study area: {name}. Its boundary is {boundary}.
Everything outside the study area is no data.
Coordinate reference system: WGS 84 / UTM zone {zone} (EPSG:{epsg}).

LAYER GROUPS

The project opens on the water bodies and the agricultural field systems, drawn over the dry-season
Sentinel-2 true colour. The other layers can be switched on in the Layers panel.

1. Interpretation (results): our interpretation, in the colours and with the legends of the maps in the deck.
   - Study area: the outline of the study area.
   - Water bodies by form: one polygon per water body (water_bodies.gpkg), with its form class and shape
     measures in the attribute table, and the same classes as a 10 m raster (water_class_10m.tif). The
     table is also given as water_features.csv.
   - Orientation of rectilinear water: the rectilinear water bodies that lie within 7.5 degrees of
     north-south or east-west.
   - Ponds and small water bodies per km²: rectilinear ponds, irregular ponds and small water bodies, all
     under 10 ha, counted at their centroid in a 1 km moving window.
   - Water, dry versus wet season, and months with open water in {year}: water seen by the Sentinel-1 radar.
   - Dominant field-boundary direction: how far field edges deviate from the cardinal directions, in
     100 m cells.
   - NDVI, wet minus dry season: the change in the Sentinel-2 vegetation index between the two seasons.
   - Crop cycles on agricultural land: the number of crop cycles in {year}, from monthly Sentinel-2 NDVI.
   - Share of land under agriculture: the share of agricultural land within 1 km, which shows where
     agriculture is concentrated.
   - Agricultural field systems: flood-recession fields, irrigated or multi-cropped land, rainfed rice-field
     landscape and open agricultural land.

2. Open data (inputs): the datasets named in the task, clipped to the study area.
   - Fields of The World {year}: the reference field map.
   - ESA WorldCover 2021: land cover at 10 m.
   - JRC Global Surface Water, occurrence: the share of observations between 1984 and 2020 in which each
     pixel held water.
   - JRC Global Surface Water, transitions: the change in each pixel's water between the first and the last
     year of the record, drawn in the groups of the long-term water map (layout 6).
   - Sentinel-1 VV backscatter in dB: the median of the dry season, and the lowest monthly median of the wet
     season, which is the state that the wet-season water map is drawn from.
   - Sentinel-2 true colour: median composites of the dry season ({season_text(site, 'dry')}) and
     of the wet season ({season_text(site, 'wet')}).
   - Copernicus DEM: a hillshade of the 30 m elevation model, kept as elevation in metres on a 20 m grid.

MAP LAYOUTS

The project holds {len(sheets)} print layouts, one per analysis, on A3 {page} pages (Project > Layout
Manager). Each has a title; the map, over the faded dry-season Sentinel-2 true colour, with the study-area
outline and a graticule in degrees; a legend that gives the area or count of each class where the
products give one; a north arrow; a scale bar; a short analysis with its key numbers; the data sources
with their dates; the coordinate reference system; and a credit line. The text is set in Segoe UI; on a
computer without it, QGIS uses another font, for which the text boxes leave room. The map of each layout
follows the map theme of the same name, which can also be shown in the main window (Layers panel, Manage
Map Themes); the faded background and the other ways of drawing the water bodies are named styles of those
layers (Layer Properties, Style). The layouts are exported as PDF, with the raster images at {LAYOUT_DPI} dpi,
to the {LAYOUT_DIR} folder, and all of them together to {LAYOUT_DIR}/{site.key}_map_atlas.pdf:

{layouts}

NOTES ON READING THE LAYERS

- Form classes such as "Enclosing moat (candidate)" and "Straight canal (candidate)" describe the shape of
  a water body. They are morphological candidates and do not date a feature.
- Areas are mapped areas, not error-adjusted estimates.
- The value layers are drawn in the classes of their legends, as on the maps and the web GIS. As on the
  maps, 1 km windows without a pond and agricultural shares under 20% are not drawn, and these pixels
  read as no data.
- NDVI change, radar backscatter and agricultural share are stored as small integers with a scale factor,
  which QGIS applies, so pixel values read in NDVI units, in dB and in per cent. NDVI change is stored in
  steps of 0.01 and agricultural share in whole per cent, whereas the statistics and maps were computed at
  full precision. Each value is rounded to the nearest step unless that would move it into another class,
  in which case it takes the step on its own side of the class break, so every pixel keeps the class it
  has on the maps and lies within one step of the computed value. Pixels without a clear wet-season view
  are left blank.
- The true colour and radar images are JPEG-compressed for viewing, so their pixel values are approximate.{notes}

STATISTICS

All statistics, including the water and agricultural areas, the seasonal changes and the agreement with
reference maps, are reported in the deck and in the comparative dataset (waterfields_comparative_{year}.csv,
which is included in this folder). This folder holds the maps behind them.

DATA SOURCES AND LICENCES

- Sentinel-2 L2A and Sentinel-1 RTC: contains modified Copernicus Sentinel data {year}, accessed through
  Microsoft Planetary Computer.
- JRC Global Surface Water v1.3, 1984 to 2020 (Pekel et al. 2016): EC JRC and Google, free to use with
  attribution.
- ESA WorldCover 2021 v200 (Zanaga et al. 2022): CC BY 4.0.
- Copernicus DEM GLO-30: produced using Copernicus WorldDEM-30, (c) DLR e.V. 2010-2014 and (c) Airbus Defence
  and Space GmbH 2014-2018, provided under COPERNICUS by the European Union and ESA.
- Fields of The World {year} (Robinson et al. 2026): CC BY 4.0.
- Study area: the SRT AOI.kml supplied with the task, or for Angkor a 44.72 km square.
- Derived layers: Muhammad Dita Farel, SRTSR0576 task assignment. They keep the terms of the data they were made from,
  and the code that produced them is available on request.
"""


def ftw_outline_note(key, count):
    """README note for a package that holds the Fields of The World outlines instead of the polygons."""
    return (f"\n- Fields of The World: the {count:,} field polygons of this study area would take the data layers\n"
            f"  over {MAX_ZIP_MIB} MiB, so their outlines are given instead as a 10 m raster\n"
            f"  (ftw_field_outlines_10m.tif). The polygons, with their areas and confidence, can be streamed\n"
            f"  from Source Cooperative as GeoParquet, from the one file in the folder of the country code:\n"
            f"  {ftw_source(key)}")


def finish(site, out, ftw_note=""):
    """The project with its layouts, README.txt and a fresh copy of the comparative CSV."""
    make_project(site, out)
    (out / "README.txt").write_text(readme(site, out, ftw_note), encoding="utf-8")
    table = config.OUT / "dataset" / f"waterfields_comparative_{config.YEAR}.csv"
    if table.exists():
        shutil.copy(table, out / table.name)


def zip_folder(key, out):
    """Zip a site folder, with its layouts folder, into out/gis_zip/, inside a folder named like the zip.

    Returns the path and the size in MiB.
    """
    ZIP_DIR.mkdir(parents=True, exist_ok=True)
    path = ZIP_DIR / f"{key}_waterfields_qgis.zip"
    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED, compresslevel=9) as zf:
        for f in sorted(out.rglob("*")):
            if f.is_file() and not f.name.endswith(".aux.xml"):  # statistics a GIS may have cached beside a raster
                zf.write(f, f"{path.stem}/{f.relative_to(out).as_posix()}")
    return path, path.stat().st_size / 2**20


def build(key):
    """The QGIS folder and zip of one site."""
    site = Site(key)
    out = config.OUT / "gis" / key
    shutil.rmtree(out, ignore_errors=True)  # rebuilt from scratch, so no file of an earlier version stays
    out.mkdir(parents=True)
    write_inputs(site, out)
    write_results(site, out)
    fields = ftw.read_fields(key, site.g10.crs, shapely.Polygon(aoi_utm(key)))
    write_ftw_polygons(site, out, fields)
    finish(site, out)
    path, mib = zip_folder(key, out)  # the layers alone, before the PDFs of the layouts are added
    if mib > MAX_ZIP_MIB:  # the polygons do not fit, so their outlines go in instead
        (out / "ftw_fields_2024.gpkg").unlink()
        write_ftw_outlines(site, out, fields)
        finish(site, out, ftw_outline_note(key, len(fields[0])))
        path, mib = zip_folder(key, out)
    if mib > MAX_ZIP_MIB:
        sys.exit(f"{path} is {mib:.1f} MiB without the layouts, over the {MAX_ZIP_MIB} MiB limit")
    pdfs = export_layouts(key, out)
    path, mib = zip_folder(key, out)
    form = "polygons" if (out / "ftw_fields_2024.gpkg").exists() else "outlines"
    with zipfile.ZipFile(path) as zf:  # what each part takes in the zip
        stored = {info.filename.split("/", 1)[1]: info.compress_size / 2**20 for info in zf.infolist()}
    layouts = sum(v for name, v in stored.items() if name.startswith(f"{LAYOUT_DIR}/"))
    atlas = stored.get(f"{LAYOUT_DIR}/{key}_map_atlas.pdf", 0.0)
    print(key, "->", out / f"{key}.qgz", "|", path, f"{mib:.1f} MiB, of which the {len(pdfs)} PDFs of the layouts "
          f"{layouts:.1f} MiB (the atlas {atlas:.1f} MiB); Fields of The World as {form}", flush=True)
    if mib > MAX_ZIP_MIB:
        print(f"  over {MAX_ZIP_MIB} MiB: without the atlas the zip would take {mib - atlas:.1f} MiB, and with the "
              f"{LAYOUT_DIR} folder in a zip of its own {mib - layouts:.1f} MiB", flush=True)


if __name__ == "__main__":
    os.environ["GDAL_PAM_ENABLED"] = "NO"  # no .aux.xml files beside the rasters
    QgsApplication.setPrefixPath("/usr", True)
    app = QgsApplication([], False)
    app.initQgis()
    for k in sys.argv[1:] or list(config.SITES):
        build(k)
    app.exitQgis()
