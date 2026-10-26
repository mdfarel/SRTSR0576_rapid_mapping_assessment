"""Static maps and charts (PNG) for the deck, one set per site, plus the cross-site figures.

Inputs: the products in data/<site>/products/ (GeoTIFFs and *_stats.json), and the cached
Sentinel-2 dry-season composite and Copernicus DEM for the basemap.

Outputs, in out/figures/<site>/, for each map <name>:
    <name>.png           the map with its legend and notes in a side panel
    <name>_map.png       the map alone, for the deck, which draws the legend itself
    <name>_legend.json   that legend: title, items (colour, label), note and study-area source
plus the charts field_water_relation.png, phenology.png and rose_*.png. In out/figures/: the
compare_*.png cross-site charts (compare_water_sources.png once pipeline.sources has run at every site).

Every map is drawn on the 20 m site grid at one image pixel per grid cell, with its longest side
10 inches. The maps carry the study-area outline, a scale bar and a north arrow, sized for the
space the deck gives the map so that they read at about 10 pt on the slide, and no place names:
this is a rapid interpretation of the landscape, not a study of named sites. The area outside the
study polygon is faded. On the two water-body maps, bodies under 1 ha are drawn enlarged so that
they stay visible at slide size (see small_bodies_grown). Nothing here changes a product or a
statistic.

    python -m pipeline.maps [site ...]

The cross-site charts are drawn by `python -m pipeline.dataset` (make dataset).
"""
import json
import sys
import textwrap
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib import font_manager, patheffects
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Patch
from scipy import ndimage

import config
from pipeline.change import GROUPS, NDVI_CHANGE
from pipeline.common import aoi_utm
from pipeline.layers import Site
from pipeline.sources import MIN_GROUP_KM2
from pipeline.water import CLASSES, RECTILINEAR

for path in ("/mnt/c/Windows/Fonts/segoeui.ttf", "/mnt/c/Windows/Fonts/segoeuib.ttf",
             "/mnt/c/Windows/Fonts/seguisb.ttf", "/mnt/c/Windows/Fonts/segoeuii.ttf"):
    if Path(path).exists():
        font_manager.fontManager.addfont(path)
if any(f.name == "Segoe UI" for f in font_manager.fontManager.ttflist):
    plt.rcParams["font.family"] = "Segoe UI"

# Colours. The class colours are shared with the web map (tiles.py) and the QGIS project.
# Each colour keeps one meaning across the maps: navy is water present all year, light blue is
# seasonal inundation (wet season only), and the two purples are the small non-rectilinear bodies.
# Classes drawn over the light basemap are dark or saturated enough to show at slide size.
INK = "#1d1d1f"
MUTED = "#6e6e73"
ACCENT = "#c2410c"
LINE = "#d2d2d7"

SEASONAL = "#86b6e6"
WATER_COLORS = {
    1: "#3b7fc4", 2: "#0b7a75", 3: ACCENT, 4: "#14a39c", 5: "#5cc8c2",
    6: "#1f4e89", 7: "#6a3d9a", 8: "#a75cc0",
}
WATER_LABELS = {
    1: "River or stream", 2: "Straight canal (candidate)", 3: "Enclosing moat (candidate)",
    4: "Rectilinear reservoir or moat segment, 10 ha+", 5: "Rectilinear pond", 6: "Lake or irregular reservoir",
    7: "Irregular pond", 8: "Small water body, shape not resolved",
}
FIELD_COLORS = {1: "#178c96", 2: "#2e8540", 3: "#b5c64b", 4: "#e8dccb"}
FIELD_LABELS = {1: "Flood-recession fields", 2: "Irrigated or multi-cropped",
                3: "Rice-field landscape (rainfed)", 4: "Open agricultural land"}
CHANGE_COLORS = {1: "#1f4e89", 2: SEASONAL, 3: "#e07b39"}
CHANGE_LABELS = {1: "Water in both seasons", 2: "Seasonal inundation (wet only)", 3: "Water in dry season only"}
CYCLE_COLORS = {0: "#a3876a", 1: "#b5c64b", 2: "#2e8540", 3: "#0f4d2a"}
CYCLE_LABELS = {0: "No clear crop cycle", 1: "One crop cycle", 2: "Two crop cycles", 3: "Three crop cycles"}
GRAIN_COLORS = (ACCENT, "#e8a27c", "#c9c3bc", "#8e9aaf")  # within 7.5°, 7.5 to 15°, 15 to 30°, 30 to 45°
# The groups of the long-term water map. pipeline.change counts water present in both 1984 and 2020 as one
# group; the map splits it into water permanent in both years (JRC class 1) and water seasonal in at least
# one of them (4 seasonal, 7 seasonal to permanent, 8 permanent to seasonal), so that permanent and
# seasonal water are not merged. New, lost and some-years-only (JRC ephemeral) water are change.GROUPS.
JRC_BOTH_YEARS = {"permanent 1984 and 2020": (1,), "seasonal 1984 and 2020": (4, 7, 8)}
assert sorted(c for codes in JRC_BOTH_YEARS.values() for c in codes) == sorted(GROUPS["water since 1984"])
TRAJECTORY_CODES = {**JRC_BOTH_YEARS, **{g: GROUPS[g] for g in ("new since 1984", "lost since 1984", "ephemeral")}}
TRAJECTORY = [  # JRC transition group, colour, legend label
    ("permanent 1984 and 2020", "#1f4e89", "Permanent, 1984 and 2020"),
    ("seasonal 1984 and 2020", "#4a8fd4", "Seasonal, 1984 and 2020"),
    ("new since 1984", ACCENT, "New: not recorded in 1984"),
    ("lost since 1984", "#e8a27c", "Lost: not recorded in 2020"),
    ("ephemeral", "#b9c6d2", "Some years only, neither 1984 nor 2020"),
]

# Mapped water on the pond-density map, in slate grey so it is not read as the darkest density.
WATER_FOR_REFERENCE = "#6b7785"
# Bars of the cross-site charts that are not a map class
NEUTRAL = "#8e9aaf"
# A typical basemap tone. Legend swatches for colours drawn part-transparent over the basemap are
# mixed with it, so that each swatch matches what the map shows.
BASEMAP_TONE = "#cfcdc9"
# Legend swatches lighter than this (luminance, 0 to 1) get a thin grey outline so they show on white.
PALE_SWATCH = 0.85
SWATCH_OUTLINE = "#b0b0b0"
# Open agricultural land is drawn in a pale sand on the maps; its line on the NDVI chart is darker.
OPEN_AG_LINE = "#a99577"

# Map layout
CELL_M = 20                           # one map pixel is one cell of the 20 m site grid
CELL_KM2 = CELL_M * CELL_M / 1e6
LONG_SIDE_IN = 10.0                   # longest side of every map image, in inches
PANEL_IN = 3.9                        # width of the legend panel beside the map
SCALE_BAR_KM = 10
MAP_TEXT_PT = 10                      # scale bar and north arrow, in points as seen on the slide
# The box (width, height in inches) the deck fits a map into, for tall maps (width under 0.8 of the
# height), wide maps (1.8 times the height or more) and the others.
TALL_MAP, WIDE_MAP = 0.8, 1.8
DECK_BOX_IN = {                       # (tall, other, wide)
    "full": ((6.4, 5.75), (6.4, 5.75), (7.6, 3.3)),       # one map on a slide (map_slide)
    "half": ((3.67, 5.75), (4.57, 3.75), (6.6, 2.8)),     # two maps on a slide (pair_slide)
    "quarter": ((2.56, 3.2), (3.7, 2.78), (3.7, 2.78)),   # appendix, three maps and a chart (grid_slide)
}
# The per-site charts are drawn at the size the deck places them (grid_slide), so their text reads
# at its stated size: about 4.2 inches wide, in a grid cell 2.7 inches high, or in the taller column
# beside tall maps.
CHART_W_IN = 4.2
CHART_H_IN = 2.7
CHART_H_TALL_IN = 4.4
CHART_TEXT_PT = 11
NOTE_CHARS = 50                       # longest line of a note under a chart
# Display thresholds (they change what a figure shows, never a statistic)
SPARSE_KM2 = 2.0          # a layer smaller than this cannot be seen at slide size, so the map says so
MIN_CLASS_KM2 = 0.5       # field classes smaller than this are left out of legends and charts
CONCENTRATED = 0.8        # share of agriculture in 1 km counted as concentrated (as in fields.py)
GRAIN_MIN_COHERENCE = 0.15  # field-grain cells less directional than this are not coloured
RANDOM_CARDINAL_PCT = 2 * 7.5 / 90 * 100  # share of directions within 7.5° of cardinal if random
RADAR_BLIND_DB = -15      # radar cannot see a site's water if VV over it stays above this (dB)
# Water bodies under SMALL_BODY_HA cover only a few map pixels and vanish once the deck shrinks a
# map, so the water-body maps draw them grown by SMALL_BODY_GROW_PX map pixels (20 m each) on every
# side. Display only: the products, statistics and web layers keep the mapped outlines.
SMALL_BODY_HA = 1.0
SMALL_BODY_GROW_PX = 2
# Where the open water in all covers under SPARSE_KM2 (Prambanan), even grown bodies are hard to find on
# the slide, so they are grown by this many map pixels instead. Display only, as above.
SPARSE_BODY_GROW_PX = 4
SMALL_BODY_NOTE = "Bodies under 1 ha drawn enlarged to show\nat this scale; areas are as mapped."

# Value layers are drawn in classes, so that every pixel takes the exact colour of one legend row.
# The web map (tiles.py, with its legends in webgis/app.js) and the QGIS project use the same classes
# and colours. Each class is ((lowest, highest) value, both included, colour, legend label).
MONTH_CLASSES = [((1, 2), "#c6dbef", "Water 1 to 2 months"), ((3, 5), "#6baed6", "3 to 5 months"),
                 ((6, 9), "#2171b5", "6 to 9 months"), ((10, 12), "#08306b", "10 to 12 months (permanent)")]
# Pond-class bodies counted in a 1 km² window are whole counts, rounded before they are classed.
DENSITY_CLASSES = [((1, 1), "#c7e9b4", "1"), ((2, 3), "#7fcdbb", "2 to 3"), ((4, 7), "#1d91c0", "4 to 7"),
                   ((8, np.inf), "#081d58", "8 or more")]
# JRC occurrence: per cent of the 1984 to 2020 observations with water (0 is never water, not drawn).
JRC_CLASSES = [((1, 24), "#c6dbef", "Under 25% of observations"), ((25, 49), "#6baed6", "25 to 49%"),
               ((50, 74), "#2171b5", "50 to 74%"), ((75, 100), "#08306b", "75% or more")]
# NDVI change, wet minus dry season, in legend order (see ndvi_change_class): (colour, legend label).
NDVI_STRONG = 0.3  # a change of this much or more either way is strongly greener or browner
NDVI_CLASSES = [("#01665e", "Strongly greener, +0.3 or more"), ("#80cdc1", "Greener, +0.15 to +0.3"),
                ("#f5f5f5", "Little change, -0.15 to +0.15"), ("#dfc27d", "Browner, -0.15 to -0.3"),
                ("#8c510a", "Strongly browner, -0.3 or less")]
NO_CLEAR_VIEW = "#c8c8cc"
# Share of agricultural land within 1 km, a continuous value: (lowest share of the class, colour, legend
# label), each class running up to the lowest share of the next; shares under 20% are not drawn.
AG_SHARE_CLASSES = [(0.2, "#fed98e", "20 to 40% agricultural"), (0.4, "#fe9929", "40 to 60%"),
                    (0.6, "#d95f0e", "60 to 80%"),
                    (CONCENTRATED, "#993404", f"{CONCENTRATED:.0%} or more (concentrated)")]


def figdir(key):
    """Output folder for one site's figures."""
    d = config.OUT / "figures" / key
    d.mkdir(parents=True, exist_ok=True)
    return d


def read_json(path):
    """Contents of a JSON file, or None when it does not exist."""
    return json.loads(path.read_text()) if path.exists() else None


# Basemap and colouring

def hillshade(dem, res, az=315, alt=45):
    """Hillshade (0 to 1) of a DEM, lit from the north-west."""
    gy, gx = np.gradient(dem.astype(np.float32), res)
    slope = np.arctan(np.hypot(gx, gy) * 2.0)  # vertical exaggeration 2
    aspect = np.arctan2(-gx, gy)
    a, z = np.radians(az), np.radians(90 - alt)
    hs = np.cos(z) * np.cos(slope) + np.sin(z) * np.sin(slope) * np.cos(a - aspect)
    return np.clip(hs, 0, 1)


def basemap(site, true_colour_weight=0.45):
    """Light grey hillshade softly tinted with the dry-season true colour, on the 20 m grid."""
    dem = site.static10("dem")[::2, ::2]
    hs = hillshade(np.where(dem < -1000, 0, dem), CELL_M)
    grey = 0.84 + 0.22 * (hs - 0.7)
    rgb = np.dstack([grey] * 3)
    d = site.s2_composite("dry", 20)
    tc = np.dstack([d["B04"], d["B03"], d["B02"]])
    tc = np.clip(tc / 0.18, 0, 1) ** 0.8
    lum = tc.mean(axis=2, keepdims=True)
    tc = 0.6 * tc + 0.4 * lum  # partly desaturate so the thematic colours stay dominant
    ok = np.isfinite(tc).all(axis=2, keepdims=True)
    rgb = np.where(ok, (1 - true_colour_weight) * rgb + true_colour_weight * np.nan_to_num(tc), rgb)
    return np.clip(rgb, 0, 1)


def paint(rgb, mask, colour, opacity=1.0):
    """Draw one colour over the image (in place) where mask is true."""
    rgb[mask] = (1 - opacity) * rgb[mask] + opacity * np.array(to_rgb(colour))


def paint_classes(rgb, values, classes, opacity):
    """Draw each ((lowest, highest), colour, label) class of a value layer (in place); NaN takes no class."""
    for (lo, hi), colour, _ in classes:
        paint(rgb, (values >= lo) & (values <= hi), colour, opacity)


def ndvi_change_class(dn):
    """The NDVI_CLASSES index of each pixel of an NDVI change (wet minus dry), or -1 where it is NaN.

    Greener and browner are changes beyond NDVI_CHANGE (more than +0.15, less than -0.15), as
    pipeline/change.py counts them; a change of NDVI_STRONG or more either way is strongly so.
    """
    d = np.nan_to_num(dn, nan=0.0)
    cls = np.select([d >= NDVI_STRONG, d > NDVI_CHANGE, d >= -NDVI_CHANGE, d > -NDVI_STRONG], [0, 1, 2, 3], 4)
    return np.where(np.isfinite(dn), cls, -1).astype(np.int8)


def over_basemap(colour, opacity):
    """The colour shown when `colour` is drawn at `opacity` over the basemap (for legend swatches)."""
    return to_hex(opacity * np.array(to_rgb(colour)) + (1 - opacity) * np.array(to_rgb(BASEMAP_TONE)))


def classes_to_20m(cls10, order):
    """Collapse a 10 m class map to 20 m, letting any present class win in the given priority order."""
    h, w = cls10.shape
    blocks = cls10[: h // 2 * 2, : w // 2 * 2].reshape(h // 2, 2, w // 2, 2)
    out = np.zeros((h // 2, w // 2), np.uint8)
    for k in reversed(order):
        out[(blocks == k).any(axis=(1, 3))] = k
    return out


def small_bodies_grown(cls10, order, grow=SMALL_BODY_GROW_PX):
    """Display only: the 20 m classes of water bodies under SMALL_BODY_HA, grown by `grow` pixels.

    Bodies are the 8-connected groups of water pixels, as in pipeline/water.py. A water map draws
    this layer under the mapped classes, so it adds a rim around each small body and leaves the
    larger bodies as mapped; no product or statistic reads it.
    """
    lab, _ = ndimage.label(cls10 > 0, np.ones((3, 3)))
    small = np.bincount(lab.ravel()) < SMALL_BODY_HA * 100  # 10 m pixels of 0.01 ha
    small[0] = False  # label 0 is the background
    cls = classes_to_20m(np.where(small[lab], cls10, 0).astype(np.uint8), order)
    grown = np.zeros_like(cls)
    square = np.ones((3, 3), bool)
    for k in reversed(order):  # where grown bodies meet, the class first in order wins
        grown[ndimage.binary_dilation(cls == k, square, iterations=grow)] = k
    return grown


def display_growth(site):
    """Map pixels by which small water bodies are grown for display: more where open water is sparse."""
    sparse = site.load_json("water_stats")["dry_season_water_km2"] < SPARSE_KM2
    return SPARSE_BODY_GROW_PX if sparse else SMALL_BODY_GROW_PX


def coarse_to_20m(grid, factor, shape, fill):
    """Repeat each cell of a coarser grid factor x factor times, cut or padded with `fill` to the map shape."""
    big = np.kron(grid, np.ones((factor, factor)))[: shape[0], : shape[1]]
    out = np.full(shape, fill, dtype=float)
    out[: big.shape[0], : big.shape[1]] = big
    return out


# Map frame: outlines, names, scale bar, north arrow, legend

def to_pixels(site, xy):
    """Pixel (column, row) positions on the 20 m grid of UTM (x, y) points."""
    t = site.g20.transform
    return [((x - t.c) / CELL_M, (t.f - y) / CELL_M) for x, y in xy]


def map_text_size(h, w, shown):
    """Font size on the 10-inch image that reads at MAP_TEXT_PT once the deck has shrunk the map."""
    if w / h < TALL_MAP:
        box_w, box_h = DECK_BOX_IN[shown][0]
    elif w / h >= WIDE_MAP:
        box_w, box_h = DECK_BOX_IN[shown][2]
    else:
        box_w, box_h = DECK_BOX_IN[shown][1]
    image_w, image_h = LONG_SIDE_IN * w / max(h, w), LONG_SIDE_IN * h / max(h, w)
    shrink = min(box_w / image_w, box_h / image_h)
    return MAP_TEXT_PT / shrink


def halo(fs, line_width=0):
    """White outline around text and lines (of the given width), so they read over any map colour."""
    return [patheffects.withStroke(linewidth=line_width + fs / 4, foreground="white")]


def draw_boundaries(ax, site, fs):
    """Study-area outline, dashed."""
    ring = to_pixels(site, aoi_utm(site.key))
    ax.plot(*zip(*ring), color=INK, lw=fs / 14, ls=(0, (4, 3)))


def draw_scale_bar(ax, h, w, em, fs):
    """A SCALE_BAR_KM bar with its length written above it, in the lower-right corner."""
    length = SCALE_BAR_KM * 1000 / CELL_M
    y = h - em
    x1 = w - em
    x0 = x1 - length
    ax.plot([x0, x1], [y, y], color=INK, lw=fs / 4, solid_capstyle="butt", path_effects=halo(fs, fs / 4), zorder=5)
    ax.text((x0 + x1) / 2, y - 0.4 * em, f"{SCALE_BAR_KM} km", ha="center", va="bottom", fontsize=fs,
            color=INK, path_effects=halo(fs), zorder=5)


def draw_north_arrow(ax, w, em, fs):
    """An N with an arrow under it pointing up, in the upper-right corner."""
    x = w - 1.2 * em
    ax.text(x, 0.6 * em, "N", ha="center", va="top", fontsize=fs, weight="bold", color=INK,
            path_effects=halo(fs), zorder=5)
    ax.annotate("", xy=(x, 2.0 * em), xytext=(x, 4.0 * em), zorder=5,
                arrowprops=dict(arrowstyle="-|>", color=INK, lw=fs / 8, mutation_scale=fs * 1.2,
                                shrinkA=0, shrinkB=0, path_effects=halo(fs, fs / 8)))


def swatch_edge(colour):
    """Outline colour of a legend swatch: thin grey for near-white swatches, none for the others."""
    luminance = np.dot(to_rgb(colour), [0.2126, 0.7152, 0.0722])
    return SWATCH_OUTLINE if luminance > PALE_SWATCH else "none"


def draw_legend_panel(fig, left, legend, title, foot):
    """Legend and notes in a panel to the right of the map (for the composed figure only)."""
    pan = fig.add_axes([left, 0, 1 - left, 1])
    pan.set_axis_off()
    pan.set_xlim(0, 1)
    pan.set_ylim(0, 1)
    if legend:
        handles = [Patch(facecolor=c, edgecolor=swatch_edge(c), linewidth=0.6, label=label) for c, label in legend]
        leg = pan.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.06, 0.98), frameon=False,
                         fontsize=11, title=title, title_fontsize=12, labelspacing=0.6, handlelength=1.3,
                         handleheight=1.0, alignment="left")
        leg.get_title().set_color(INK)
        leg.get_title().set_fontweight("bold")
    pan.text(0.08, 0.03, "\n".join(foot), fontsize=9, color=MUTED, va="bottom", ha="left", wrap=True)


def frame(site, rgb, path, legend, title, note=None, shown="full", map_note=None):
    """Save a map at one pixel per 20 m cell: alone with its legend as JSON (for the deck), and with a legend panel.

    shown is how much of a slide the deck gives the map ("full", "half" or "quarter"), which sets
    the text size. map_note, for a layer too small to see, is added to the legend note rather than
    written on the map, which carries no text but its scale bar and north arrow.
    """
    if map_note:
        map_note = map_note.replace("\n", " ")
        note = "\n".join(part for part in (note, map_note if map_note.endswith(".") else map_note + ".") if part)
    h, w = rgb.shape[:2]
    dpi = max(h, w) / LONG_SIDE_IN
    map_w, map_h = w / dpi, h / dpi
    fs = map_text_size(h, w, shown)
    em = fs * dpi / 72  # one text height, in map pixels

    rgb = rgb.copy()
    outside = ~site.mask(20)
    rgb[outside] = 0.4 * rgb[outside] + 0.6  # fade the area outside the study polygon

    fig = plt.figure(figsize=(map_w, map_h), dpi=dpi)
    ax = fig.add_axes([0, 0, 1, 1])
    ax.imshow(rgb, interpolation="none")
    ax.set_xlim(-0.5, w - 0.5)
    ax.set_ylim(h - 0.5, -0.5)
    ax.set_axis_off()
    draw_boundaries(ax, site, fs)
    draw_north_arrow(ax, w, em, fs)
    draw_scale_bar(ax, h, w, em, fs)

    stem = Path(path).with_suffix("")
    fig.savefig(f"{stem}_map.png", dpi=dpi, facecolor="white")
    spec = {"title": title, "items": legend, "note": note, "aoi": site.meta["aoi_source"]}
    Path(f"{stem}_legend.json").write_text(json.dumps(spec))

    # The same figure, widened, with the legend panel on the right.
    fig.set_size_inches(map_w + PANEL_IN, map_h)
    ax.set_position([0, 0, map_w / (map_w + PANEL_IN), 1])
    foot = ([note] if note else []) + [f"Study area: {site.meta['aoi_source']}"]
    draw_legend_panel(fig, map_w / (map_w + PANEL_IN), legend, title, foot)
    fig.savefig(path, dpi=dpi, facecolor="white")
    plt.close(fig)


# One function per map

def water_map(site, out):
    """Water bodies by form class, over seasonal (wet-season-only) water where the radar can see water."""
    cls10 = site.load("water_class_10m")
    order = [3, 2, 4, 5, 1, 6, 7, 8]  # which class wins when four 10 m cells merge: thin moats and canals first
    cls = classes_to_20m(cls10, order)
    grown = small_bodies_grown(cls10, order, display_growth(site))
    rgb = basemap(site)
    blind = radar_cannot_see_water(site)  # then the radar's seasonal water is not shown, as on the deck
    if not blind:
        seasonal = site.load("water_change_20m") == 2
        paint(rgb, seasonal & (cls == 0), SEASONAL)
    for layer in (grown, cls):  # enlarged small bodies first, then every body as mapped
        for k in reversed(order):
            paint(rgb, layer == k, WATER_COLORS[k])
    stats = site.load_json("water_stats")
    counts = stats.get("count_by_class", {})
    legend = []
    for k in [1, 6, 2, 3, 4, 5, 7, 8]:
        if (cls10 == k).any():
            legend.append((WATER_COLORS[k], f"{WATER_LABELS[k]}: {counts.get(CLASSES[k], 0):,}"))
    if blind:
        seasonal_note = "Seasonal water: not measurable by radar here."
    else:
        legend.append((SEASONAL, CHANGE_LABELS[2]))
        seasonal_note = "Seasonal water: Sentinel-1, wet season 2024."
    note = f"Open water: Sentinel-2 10 m, dry season 2024.\n{seasonal_note}\n{SMALL_BODY_NOTE}"
    if stats["dry_season_water_km2"] < SPARSE_KM2:
        note += f"\nLittle open water: {km2_text(stats['dry_season_water_km2'])} km² in all."
    frame(site, rgb, out / "water.png", legend, "Water bodies by form", note=note)


def cardinal_map(site, out):
    """Rectilinear water bodies, split by whether they lie within 7.5° of the cardinal directions."""
    cls10 = site.load("water_class_10m")
    card10 = site.load("water_cardinal_10m")
    rect = np.isin(cls10, RECTILINEAR)
    # 2 rectilinear and cardinal, 1 rectilinear in another direction, 3 any other water body
    groups10 = np.where(rect, np.where(card10 == 1, 2, 1), np.where(cls10 > 0, 3, 0)).astype(np.uint8)
    order = [2, 1, 3]
    groups = classes_to_20m(groups10, order)
    grown = small_bodies_grown(groups10, order, display_growth(site))
    rgb = basemap(site, 0.2)
    colours = {2: ACCENT, 1: "#3a5f7d", 3: "#9fb8cf"}
    for layer in (grown, groups):  # enlarged small bodies first, then every body as mapped
        for k in reversed(order):
            paint(rgb, layer == k, colours[k])
    legend = [(ACCENT, "Rectilinear, within 7.5° of N-S or E-W"),
              ("#3a5f7d", "Rectilinear, other orientation"), ("#9fb8cf", "Other water bodies")]
    frame(site, rgb, out / "water_cardinal.png", legend, "Orientation of rectilinear water",
          note="Rectilinear: outline edges follow one pair of\nperpendicular directions (coherence 0.40 or more).\n"
               f"{SMALL_BODY_NOTE}")


def density_map(site, out):
    """Pond-class bodies per km² in a 1 km window, with the mapped water bodies for reference."""
    dens = np.round(site.load("pond_density_20m"))  # whole counts in the 1 km² window
    rgb = basemap(site, 0.2)
    opacity = 0.75
    paint_classes(rgb, dens, DENSITY_CLASSES, opacity)
    water = classes_to_20m((site.load("water_class_10m") > 0).astype(np.uint8), [1])
    paint(rgb, water == 1, WATER_FOR_REFERENCE)
    legend = [(over_basemap(colour, opacity), label) for _, colour, label in DENSITY_CLASSES]
    legend.append((WATER_FOR_REFERENCE, "Mapped water body"))
    frame(site, rgb, out / "water_density.png", legend, "Ponds per km² (under 10 ha)", shown="quarter",
          note="Ponds: rectilinear and irregular ponds and\nsmall water bodies, all under 10 ha, counted\n"
               "at their centroid in a 1 km moving window.")


def months_map(site, out):
    """Number of months in 2024 with open water on Sentinel-1."""
    months = site.load("s1_water_months_20m")  # 255 is no data
    rgb = basemap(site, 0.2)
    opacity = 0.85
    ok = (months >= 1) & (months <= 12)
    paint_classes(rgb, months, MONTH_CLASSES, opacity)
    legend = [(over_basemap(colour, opacity), label) for (lo, hi), colour, label in MONTH_CLASSES
              if ((months >= lo) & (months <= hi)).any()]  # only the classes that occur
    km2 = ok.sum() * CELL_KM2
    map_note = None
    if radar_cannot_see_water(site):  # the deck's wording for such a site, with the radar figure named
        map_note = f"Not measurable by radar\n(radar estimate {km2_text(km2)} km²\nin any month of 2024, not reliable)"
    elif km2 < SPARSE_KM2:
        map_note = f"Radar detects almost\nno open water here\n({km2:.2f} km²)"
    frame(site, rgb, out / "water_months.png", legend, "Months with open water in 2024", shown="quarter",
          note=radar_threshold_note(site.load_json("water_stats")), map_note=map_note)


def radar_threshold_note(water_stats):
    """How the Sentinel-1 water threshold was set at this site, from water_stats.json threshold_method."""
    if "Otsu" in water_stats["threshold_method"]["s1_vv"]:
        return "Sentinel-1 VV monthly medians below the\nOtsu water threshold."
    db = water_stats["thresholds"]["s1_vv_db_edge_otsu"]  # the threshold applied, here the fixed fallback
    return f"Sentinel-1 VV monthly medians below a fixed\n{db:.0f} dB threshold (edges not bimodal)."


def field_classes_shown(site):
    """Field classes large enough to show, and a note on any left out.

    A class whose area rounds to 0.0 km² in field_stats.json is called practically absent, not small.
    """
    area = site.load_json("field_stats")["area_km2_by_class"]
    shown = [k for k in (1, 2, 3, 4) if area[FIELD_LABELS[k]] >= MIN_CLASS_KM2]
    absent = [FIELD_LABELS[k] for k in (1, 2, 3, 4) if area[FIELD_LABELS[k]] == 0]
    small = [f"{FIELD_LABELS[k]} ({area[FIELD_LABELS[k]]:.1f} km²)" for k in (1, 2, 3, 4)
             if k not in shown and FIELD_LABELS[k] not in absent]
    notes = []
    if small:
        notes.append(f"Left out, under {MIN_CLASS_KM2} km²: {', '.join(lower_first(n) for n in small)}.")
    if absent:
        notes.append(f"Practically none mapped: {', '.join(f'{lower_first(n)} (0.0 km²)' for n in absent)}.")
    return shown, " ".join(notes) or None


def lower_first(label):
    """A legend label with a lower-case first letter, for use after a colon in a note."""
    return label[:1].lower() + label[1:]


def km2_text(v):
    """An area as the deck writes it: one decimal, or two without a trailing zero under 1 km² (0.89, 0.4)."""
    if 0 < abs(v) < 1:
        return f"{v:.2f}".rstrip("0")
    return f"{v:,.1f}"


def area_label(site, k):
    """A field class name with its area from field_stats.json, for the charts."""
    area = site.load_json("field_stats")["area_km2_by_class"][FIELD_LABELS[k]]
    return f"{FIELD_LABELS[k]}, {area:,.1f} km²"


def fields_map(site, out):
    """Agricultural field-system classes, with dry-season water."""
    cls = site.load("field_class_20m")
    rgb = basemap(site, 0.25)
    for k in (4, 3, 2, 1):
        paint(rgb, cls == k, FIELD_COLORS[k])
    water = classes_to_20m((site.load("water_dry_10m") > 0).astype(np.uint8), [1])
    paint(rgb, water == 1, WATER_COLORS[1])
    shown, left_out = field_classes_shown(site)
    legend = [(FIELD_COLORS[k], FIELD_LABELS[k]) for k in shown] + [(WATER_COLORS[1], "Dry-season water"),
                                                                     (BASEMAP_TONE, "Other land")]
    note = "Sentinel-1 paddy cycle, Sentinel-2 dry-season\nNDVI and monthly crop cycles, ESA WorldCover."
    if left_out:
        note += "\n" + left_out
    frame(site, rgb, out / "fields.png", legend, "Agricultural field systems, 2024", note=note)


def cycles_map(site, out):
    """Number of crop cycles in 2024 on agricultural land."""
    cyc = site.load("crop_cycles_20m")
    ag = site.load("field_class_20m") > 0
    rgb = basemap(site, 0.25)
    for k in (0, 1, 2, 3):
        paint(rgb, ag & (cyc == k), CYCLE_COLORS[k])
    paint(rgb, ag & (cyc == 255), "#ffffff", 0.6)  # 255: too few clear months to count cycles
    # Classes covering less than MIN_CLASS_KM2 cannot be seen at slide size: they stay off the legend,
    # and the note names those that occur at all.
    classes = [(k, CYCLE_COLORS[k], CYCLE_LABELS[k]) for k in (1, 2, 3, 0)]
    classes.append((255, over_basemap("#ffffff", 0.6), "Too few clear months"))
    legend, left_out = [], []
    for k, colour, label in classes:
        km2 = (ag & (cyc == k)).sum() * CELL_KM2
        if km2 >= MIN_CLASS_KM2:
            legend.append((colour, label))
        elif km2 > 0:
            left_out.append(label)
    note = "Peaks in monthly Sentinel-2 NDVI (0.45 or more,\nrising 0.15 above neighbouring troughs)."
    if left_out:
        note += f"\nLeft out, under {MIN_CLASS_KM2} km²: {', '.join(lower_first(n) for n in left_out)}."
    frame(site, rgb, out / "crop_cycles.png", legend, "Crop cycles on agricultural land, 2024", shown="quarter",
          note=note)


def intensity_map(site, out):
    """Share of agricultural land within a 1 km window."""
    share = site.load("ag_intensity_20m")
    rgb = basemap(site, 0.2)
    opacity = 0.7
    # 0 below the lowest class, then 1 to 4 for the classes; NaN is drawn as no class
    cls = np.digitize(np.nan_to_num(share, nan=-1.0), [lo for lo, _, _ in AG_SHARE_CLASSES])
    for k, (_, colour, _) in enumerate(AG_SHARE_CLASSES, start=1):
        paint(rgb, cls == k, colour, opacity)
    legend = [(over_basemap(colour, opacity), label) for _, colour, label in AG_SHARE_CLASSES]
    frame(site, rgb, out / "ag_intensity.png", legend, "Share of land under agriculture, 1 km", shown="quarter")


def change_maps(site, out):
    """Two maps of wet against dry season: water (Sentinel-1) and vegetation (Sentinel-2 NDVI)."""
    wc = site.load("water_change_20m")
    rgb = basemap(site, 0.2)
    for k in (3, 2, 1):
        paint(rgb, wc == k, CHANGE_COLORS[k])
    legend = [(CHANGE_COLORS[k], CHANGE_LABELS[k]) for k in (1, 2, 3) if (wc == k).any()]
    km2 = np.isin(wc, (1, 2, 3)).sum() * CELL_KM2
    map_note = None
    if radar_cannot_see_water(site):  # the deck's wording for such a site, with the radar figure named
        map_note = f"Not measurable by radar\n(radar estimate {km2_text(km2)} km², not reliable)"
    elif km2 < SPARSE_KM2:
        map_note = f"Radar detects almost\nno water here\n({km2:.2f} km²)"
    frame(site, rgb, out / "change_water.png", legend, "Water, dry versus wet season", shown="half",
          note="Sentinel-1 VV, monthly medians.", map_note=map_note)

    # "Greener" or "browner" in the deck counts both shades on its side (a change of more than 0.15),
    # so the labels give each shade's range.
    ndvi = ndvi_change_class(site.load("ndvi_change_20m"))
    rgb = basemap(site, 0.0)
    opacity = 0.85
    for k, (colour, _) in enumerate(NDVI_CLASSES):
        paint(rgb, ndvi == k, colour, opacity)
    paint(rgb, ndvi < 0, NO_CLEAR_VIEW, 0.4)
    legend = [(over_basemap(colour, opacity), label) for colour, label in NDVI_CLASSES]
    legend.append((over_basemap(NO_CLEAR_VIEW, 0.4), "No clear wet-season view"))
    note = "Sentinel-2 NDVI difference of seasonal medians.\nGreener or browner: a change of more than 0.15."
    if (wc == 1).any():  # lakes and rivers take a colour too, but NDVI over open water is not about plants
        note += "\nNot interpreted over open water (lakes, rivers)."
    frame(site, rgb, out / "change_ndvi.png", legend, "NDVI, wet minus dry season", shown="half", note=note)


def grain_map(site, out):
    """Dominant direction of field boundaries in 100 m cells, as its deviation from the cardinal directions."""
    dev = site.load("field_grain_deviation_100m")
    coh = site.load("field_grain_coherence_100m")
    rgb = basemap(site, 0.15)
    shape = rgb.shape[:2]
    deviation = coarse_to_20m(np.nan_to_num(dev, nan=-1), 5, shape, -1)
    coherence = coarse_to_20m(np.nan_to_num(coh), 5, shape, 0)
    has = (deviation >= 0) & (coherence >= GRAIN_MIN_COHERENCE)
    palette = np.array([to_rgb(c) for c in GRAIN_COLORS])
    colours = palette[np.digitize(deviation, [7.5, 15, 30])]
    rgb[has] = 0.2 * rgb[has] + 0.8 * colours[has]
    labels = ["Within 7.5° of cardinal", "7.5° to 15°", "15° to 30°", "30° to 45° (diagonal)"]
    legend = [(over_basemap(c, 0.8), label) for c, label in zip(GRAIN_COLORS, labels)]
    frame(site, rgb, out / "field_grain.png", legend, "Dominant field-boundary direction",
          note="Structure tensor of Sentinel-2 10 m NDVI\nover cropland, 1 km window.")


def trajectory_map(site, out):
    """Long-term water change, 1984 to 2020, from the JRC Global Surface Water transitions.

    Where every JRC group covers less than SPARSE_KM2 the map looks empty at slide size, so it says
    so, with the areas from change_stats.json.
    """
    if not (site.static / "jrc_transitions.tif").exists():
        return
    tr = site.static10("jrc_transitions")
    rgb = basemap(site, 0.25)
    for group, colour, _ in reversed(TRAJECTORY):
        m = classes_to_20m(np.isin(tr, TRAJECTORY_CODES[group]).astype(np.uint8), [1]) == 1
        paint(rgb, m, colour)
    legend = [(colour, label) for _, colour, label in TRAJECTORY]
    area = site.load_json("change_stats")["long_term"]["area_km2"]
    map_note = None
    if max(area.values()) < SPARSE_KM2:
        map_note = (f"JRC records little water here:\n{km2_text(area['water since 1984'])} km² in 1984 and 2020, "
                    f"{km2_text(area['new since 1984'])} km² new,\n{km2_text(area['lost since 1984'])} km² lost, "
                    f"{km2_text(area['ephemeral'])} km² in some years only")
    frame(site, rgb, out / "trajectory.png", legend, "Water, 1984 to 2020", shown="quarter",
          note="JRC Global Surface Water transitions\n(Pekel et al. 2016), 30 m.", map_note=map_note)


# Per-site charts, drawn at the size the deck places them (CHART_W_IN wide) with 11 pt text.

def style_axes(ax):
    """Plain axes: no top or right spine, light grey left and bottom spines."""
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(LINE)


def share_in_band(dist, mask, lo, hi):
    """Percentage of the masked pixels whose distance lies in [lo, hi)."""
    return 100 * float(((dist >= lo) & (dist < hi) & mask).sum()) / max(int(mask.sum()), 1)


def chart_size(site):
    """(width, height) in inches of a per-site chart: taller in the column beside tall maps."""
    tall = site.g20.width / site.g20.height < TALL_MAP
    return CHART_W_IN, CHART_H_TALL_IN if tall else CHART_H_IN


def relation_chart(site, out):
    """Share of each field class within distance bands of perennial water, against all land."""
    dist = site.load("distance_to_perennial_water_20m")
    fields = site.load("field_class_20m")
    inside = site.mask(20)
    shown, left_out = field_classes_shown(site)
    groups = [("All land", inside & (dist > 0))] + [(FIELD_LABELS[k], inside & (fields == k)) for k in shown]
    bands = [(0, 500, "under 0.5 km", "#1f4e89"), (500, 1000, "0.5 to 1 km", "#5b8fc7"),
             (1000, 2000, "1 to 2 km", "#a9cbea"), (2000, 1e9, "over 2 km", "#c5cfd8")]
    fig, ax = plt.subplots(figsize=chart_size(site), dpi=200, layout="constrained")
    y = np.arange(len(groups))
    left = np.zeros(len(groups))
    for lo, hi, label, colour in bands:
        vals = np.array([share_in_band(dist, mask, lo, hi) for _, mask in groups])
        ax.barh(y, vals, left=left, color=colour, height=0.65, label=label, edgecolor="white", linewidth=0.5)
        left += vals
    ax.set_yticks(y, [name for name, _ in groups], fontsize=CHART_TEXT_PT, color=INK)
    ax.invert_yaxis()
    ax.set_xlim(0, 100)
    ax.set_xticks([0, 50, 100], ["0", "50%", "100%"])
    ax.tick_params(axis="x", labelsize=CHART_TEXT_PT, colors=MUTED)
    ax.tick_params(axis="y", length=0)
    style_axes(ax)
    fig.legend(fontsize=CHART_TEXT_PT, frameon=False, ncol=2, loc="outside upper left", handlelength=1,
               columnspacing=1.0, title="Distance to perennial water", title_fontsize=CHART_TEXT_PT,
               alignment="left")
    if left_out:  # under the whole figure, as the class names take the left of the axes
        fig.supxlabel(textwrap.fill(left_out, NOTE_CHARS), x=0, ha="left", fontsize=CHART_TEXT_PT, color=MUTED)
    fig.savefig(out / "field_water_relation.png", dpi=200, transparent=True)
    plt.close(fig)


def phenology_chart(site, out):
    """Mean monthly NDVI of each field class in 2024, with the dry and wet seasons shaded."""
    curves = site.load_json("field_stats")["ndvi_monthly_mean_by_class"]
    shown, left_out = field_classes_shown(site)
    line_colours = {**FIELD_COLORS, 4: OPEN_AG_LINE}
    fig, ax = plt.subplots(figsize=chart_size(site), dpi=200, layout="constrained")
    months = np.arange(1, 13)
    for k in shown:
        vals = curves.get(FIELD_LABELS[k])
        if vals and all(v is not None for v in vals):
            ax.plot(months, vals, color=line_colours[k], lw=2.2, label=area_label(site, k))
    # Season labels sit above the plot, where no line can cross them.
    for season, (start, end) in site.meta["seasons"].items():
        m0, m1 = int(start[5:7]), int(end[5:7])
        ax.axvspan(m0 - 0.5, m1 + 0.5, color="#1f4e89" if season == "wet" else "#e07b39", alpha=0.07, lw=0)
        ax.text((m0 + m1) / 2, 1.02, f"{season} season", transform=ax.get_xaxis_transform(), ha="center",
                va="bottom", fontsize=CHART_TEXT_PT, color=MUTED, clip_on=False)
    ax.set_xticks(months, ["J", "F", "M", "A", "M", "J", "J", "A", "S", "O", "N", "D"], fontsize=CHART_TEXT_PT,
                  color=INK)
    ax.set_xlim(0.5, 12.5)
    ax.set_ylabel("Mean NDVI", fontsize=CHART_TEXT_PT, color=MUTED)
    ax.yaxis.set_major_locator(plt.MaxNLocator(4))
    ax.tick_params(axis="y", labelsize=CHART_TEXT_PT, colors=MUTED)
    style_axes(ax)
    fig.legend(fontsize=CHART_TEXT_PT, frameon=False, loc="outside lower left", ncol=1, handlelength=1.5)
    if left_out:
        ax.set_xlabel(textwrap.fill(left_out, NOTE_CHARS), loc="left", fontsize=CHART_TEXT_PT, color=MUTED)
    fig.savefig(out / "phenology.png", dpi=200, transparent=True)
    plt.close(fig)


def rose(hist, path):
    """Half-rose of boundary azimuths (0 to 180°, 5° bins), mirrored to a full circle.

    Orange petals are bins lying wholly within 7.5° of north-south or east-west.
    """
    h = np.asarray(hist)
    centres = np.arange(0, 360, 5) + 2.5
    fig = plt.figure(figsize=(4, 4), dpi=200)
    ax = fig.add_subplot(projection="polar")
    ax.set_theta_zero_location("N")
    ax.set_theta_direction(-1)
    colours = [ACCENT if min(d % 90, 90 - d % 90) < 7.5 else "#9aa5b1" for d in centres]
    ax.bar(np.radians(centres), np.concatenate([h, h]), width=np.radians(5), color=colours,
           edgecolor="white", linewidth=0.4)
    ax.set_yticklabels([])
    ax.set_xticks(np.radians([0, 90, 180, 270]))
    ax.set_xticklabels(["N", "E", "S", "W"], fontsize=22, color=INK)
    ax.grid(color=LINE, lw=0.6)
    ax.spines["polar"].set_visible(False)
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


def run(key):
    """All maps and charts of one site."""
    site = Site(key)
    out = figdir(key)
    water_map(site, out)
    cardinal_map(site, out)
    density_map(site, out)
    months_map(site, out)
    fields_map(site, out)
    cycles_map(site, out)
    intensity_map(site, out)
    change_maps(site, out)
    grain_map(site, out)
    phenology_chart(site, out)
    relation_chart(site, out)
    trajectory_map(site, out)
    rose(site.load_json("field_stats")["field_orientation"]["histogram_5deg"], out / "rose_fields.png")
    ftw = read_json(site.products / "ftw_stats.json")
    if ftw:
        rose(ftw["orientation"]["histogram_5deg"], out / "rose_ftw.png")
    print(key, "figures in", out, flush=True)


# Cross-site figures

def radar_cannot_see_water(site):
    """True when Sentinel-1 backscatter over optically mapped water stays well above the open-water range.

    At Prambanan VV over the mapped water is about -11 dB every month, so radar water there is not a measure.
    """
    acc = read_json(site.products / "accuracy_stats.json") or {}
    monthly = (acc.get("s1_vv_over_optical_water_db") or {}).get("monthly_median")
    if not monthly:
        return False
    return np.nanmedian([np.nan if v is None else v for v in monthly]) > RADAR_BLIND_DB


def site_rows():
    """Per-site statistics used by the cross-site charts, in config.SITES order."""
    short = {"angkor": "Angkor", "srei_santhor": "Srei Santhor", "sukhothai": "Sukhothai", "prambanan": "Prambanan"}
    rows = []
    for key in config.SITES:
        site = Site(key)
        fields = site.load_json("field_stats")
        rotation = fields["field_orientation"]["controls"]["rotation_test"]
        rows.append({
            "name": short[key],
            "water": site.load_json("water_stats"),
            "fields": fields,
            "ftw": read_json(site.products / "ftw_stats.json"),
            "sources": read_json(site.products / "water_sources.json"),
            "radar_blind": radar_cannot_see_water(site),
            # field edges stay near the image axes after a 30° rotation: the signal is not in the ground
            "rotation_fails": rotation["after_30deg_rotation_image_axes_pct"] >= rotation["window_share_pct"],
        })
    return rows


def bar_panel(ax, names, values, title, unit, colour="#3a5f7d", ref=None, digits=1, notes=None, ci=None,
              control=None):
    """Horizontal bars, one per site, with the value written after each bar.

    notes replaces a site's bar and value with a short text, ref draws a dashed reference line, ci
    draws 95% interval lines (with end ticks) and control draws a circle for a control value. A note
    starts just right of the reference line, so the line never runs through it.
    """
    notes = notes or {}
    ci = ci or [None] * len(values)
    control = control or [None] * len(values)
    y = np.arange(len(values))
    vals = [0 if (v is None or i in notes) else v for i, v in enumerate(values)]
    ax.barh(y, vals, color=colour, height=0.6)
    for i in range(len(values)):
        if i in notes:
            ax.annotate(notes[i], (ref or 0, i), xytext=(6, 0), textcoords="offset points", va="center",
                        fontsize=12, color=INK)
            continue
        if ci[i]:
            ax.plot(ci[i], [i, i], color=INK, lw=1.2, marker="|", ms=8, mew=1.2)
        if control[i] is not None:
            ax.plot([control[i]], [i], marker="o", ms=7, mfc="white", mec=INK, mew=1.2, ls="none")
        end = max(vals[i], ci[i][1] if ci[i] else 0, control[i] if control[i] is not None else 0)
        shown = digits + 1 if values[i] is not None and digits == 1 and 0 < abs(values[i]) < 1 else digits
        text = "n/a" if values[i] is None else f"{values[i]:,.{shown}f}{unit}"  # 0.08%, as the deck writes it
        ax.annotate(text, (end, i), xytext=(8, 0), textcoords="offset points", va="center", fontsize=12, color=INK)
    if ref is not None:
        ax.axvline(ref, color=MUTED, lw=1, ls=(0, (3, 2)), zorder=0)
    ax.set_yticks(y, names, fontsize=12, color=INK)
    ax.invert_yaxis()
    ax.set_title(title, loc="left", fontsize=13, color=INK, pad=8)
    for side in ("top", "right", "bottom"):
        ax.spines[side].set_visible(False)
    ax.spines["left"].set_color(LINE)
    ax.tick_params(axis="x", bottom=False, labelbottom=False)
    ax.tick_params(axis="y", length=0)
    ax.set_xlim(0, max(max(vals), ref or 0) * 1.4 or 1)


def core_chart(rows, path):
    """Six headline measures side by side for the four sites."""
    names = [r["name"] for r in rows]
    water = [r["water"] for r in rows]
    fig, axs = plt.subplots(3, 2, figsize=(11, 8.0), dpi=200)
    # Class colours only where a panel means what that colour means on the maps (dry-season water,
    # seasonal inundation); every other panel is neutral, so no bar reads as a map class.
    bar_panel(axs[0, 0], names, [w["dry_season_water_share_pct"] for w in water],
              "Open water, dry season (% of area)", "%", WATER_COLORS[1])
    bar_panel(axs[0, 1], names, [w["seasonal_inundation_km2"] / w["study_area_km2"] * 100 for w in water],
              "Seasonal inundation (% of area)", "%", SEASONAL,
              notes={i: "not measurable by radar" for i, r in enumerate(rows) if r["radar_blind"]})
    bar_panel(axs[1, 0], names, [r["fields"]["agricultural_share_pct"] for r in rows],
              "Agricultural land (% of area)", "%", NEUTRAL)
    bar_panel(axs[1, 1], names, [w["ponds_per_100km2"] for w in water], "Ponds per 100 km²", "", NEUTRAL)
    cycles = [r["fields"]["cropping_cycles_share_pct"] for r in rows]
    bar_panel(axs[2, 0], names, [c["2"] + c["3"] for c in cycles], "Agricultural land cropped twice or more (%)",
              "%", NEUTRAL)
    bar_panel(axs[2, 1], names, [r["ftw"]["median_field_ha"] if r["ftw"] else None for r in rows],
              "Median field size, Fields of The World (ha)", " ha", NEUTRAL, digits=2)
    fig.tight_layout(h_pad=2.5, w_pad=3)
    fig.savefig(path, dpi=200, transparent=True)
    plt.close(fig)


def water_sources_chart(rows, path):
    """Irrigated or multi-cropped share of farmland by the water source within reach, per site.

    One bar per water-source group of pipeline.sources, the dashed line across each site is the share
    over all its farmland, and the thin vertical line on the two compared groups is the range over the
    settings of the sweep in which that group holds 1 km² or more of farmland (all 27, or fewer where a
    setting leaves the group smaller). A group with under 1 km² of farmland is marked "none".
    """
    names = [r["name"] for r in rows]
    groups = [("reservoir reach only", WATER_COLORS[4]), ("natural water only", WATER_COLORS[1]),
              ("both", NEUTRAL), ("neither", "#d6d1ca")]
    fig, ax = plt.subplots(figsize=(8.5, 3.6), dpi=200)
    x = np.arange(len(rows))
    width = 0.2
    for i, (group, colour) in enumerate(groups):
        offset = (i - 1.5) * width
        vals = []
        for r in rows:
            g = r["sources"]["groups"][group]
            vals.append(g["irrigated_share_pct"] if g["farmland_km2"] >= MIN_GROUP_KM2 else np.nan)
        ax.bar(x + offset, vals, width, color=colour, label=group)
        for j, (v, r) in enumerate(zip(vals, rows)):
            if np.isnan(v):  # no farmland in this group (no rectilinear reservoir of 10 ha or more)
                # set upright, so the word fits within its own bar's width and clear of its neighbours
                ax.annotate("none", (j + offset, 0), xytext=(0, 4), textcoords="offset points", ha="center",
                            va="bottom", rotation=90, fontsize=11, color=INK)
                continue
            sweep = ((r["sources"].get("sensitivity") or {}).get("irrigated_share_pct_range") or {}).get(group)
            if sweep:
                ax.plot([j + offset] * 2, sweep, color=INK, lw=1.1, marker="_", ms=6, mew=1.1)
    for j, r in enumerate(rows):
        ax.plot([j - 0.45, j + 0.45], [r["sources"]["irrigated_share_all_farmland_pct"]] * 2, color=INK, lw=1,
                ls=(0, (3, 2)), label="all farmland" if j == 0 else None)
    ax.set_xticks(x, names, fontsize=12, color=INK)
    ax.set_ylabel("Irrigated or multi-cropped\n(% of group's farmland)", fontsize=11, color=MUTED)
    ax.tick_params(axis="y", labelsize=11, colors=MUTED)
    ax.tick_params(axis="x", length=0)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    for side in ("left", "bottom"):
        ax.spines[side].set_color(LINE)
    ax.legend(fontsize=11, frameon=False, ncol=5, loc="upper center", bbox_to_anchor=(0.5, 1.2), handlelength=1)
    sweeps = [(r["sources"].get("sensitivity") or {}) for r in rows]
    settings = max((s.get("settings_total") or 0) for s in sweeps)
    used = sorted({n for s in sweeps for n in (s.get("settings_with_group") or {}).values() if n})
    counts = f"{used[0]} or {used[-1]}" if len(used) > 1 else f"{used[0]}" if used else "0"
    ax.text(0, -0.16, "Thin line: range over the settings of distance and drop in which the group holds 1 km² or "
            f"more of farmland\n({counts} of the {settings} settings). Dashed: all farmland of the site.",
            transform=ax.transAxes, fontsize=11, color=MUTED, va="top", linespacing=1.3)
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


def orientation_chart(rows, path):
    """Share of rectilinear water, field edges and Fields of The World fields within 7.5° of cardinal."""
    names = [r["name"] for r in rows]
    fig, axs = plt.subplots(1, 3, figsize=(13, 3.0), dpi=200)
    passed = [(r["water"].get("orientation_control") or {}).get("passed") for r in rows]
    bar_panel(axs[0], names, [r["water"]["cardinal_share_of_rectilinear_pct"] for r in rows],
              "Rectilinear water within 7.5° of N-S/E-W", "%", ACCENT, ref=RANDOM_CARDINAL_PCT,
              notes={i: ("fails the rotation test" if p is False else "too few bodies to test")
                     for i, p in enumerate(passed) if p is not True})
    fo = [r["fields"]["field_orientation"] for r in rows]
    bar_panel(axs[1], names, [f["cardinal_share_pct"] for f in fo],
              "Strong field edges within 7.5° (S-2 NDVI)", "%", ACCENT, ref=RANDOM_CARDINAL_PCT,
              ci=[f.get("cardinal_share_ci95_pct") for f in fo],
              control=[f["controls"]["non_agricultural_edges_pct"] for f in fo],
              notes={i: "fails the rotation test" for i, r in enumerate(rows) if r["rotation_fails"]})
    key = f"Line: 95% interval. Circle: edges outside farmland. Dashed: {RANDOM_CARDINAL_PCT:.1f}% if random."
    axs[1].text(0, -0.04, key, transform=axs[1].transAxes, fontsize=11, color=MUTED, va="top")
    bar_panel(axs[2], names, [r["ftw"]["orientation"]["cardinal_share_pct"] if r["ftw"] else None for r in rows],
              "Fields of The World within 7.5°", "%", ACCENT, ref=RANDOM_CARDINAL_PCT)
    fig.tight_layout(w_pad=3)
    fig.savefig(path, dpi=200, bbox_inches="tight", transparent=True)
    plt.close(fig)


def compare(path_dir=None):
    """Cross-site bar charts from the per-site stats."""
    out = path_dir or config.OUT / "figures"
    out.mkdir(parents=True, exist_ok=True)
    rows = site_rows()
    core_chart(rows, out / "compare_core.png")
    if all(r["sources"] for r in rows):  # pipeline.sources has run at every site
        water_sources_chart(rows, out / "compare_water_sources.png")
    orientation_chart(rows, out / "compare_orientation.png")


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
