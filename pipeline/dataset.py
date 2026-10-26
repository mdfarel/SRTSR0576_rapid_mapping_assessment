"""Cross-site comparative dataset: one CSV of every headline metric with its range and caveats,
plus the comparison charts. Run after every site is processed.

Inputs: the season windows in config.py, the statistics JSON files of every site (water, fields,
change, water sources, Fields of The World and agreement with the reference maps) and
water_features.csv.
Output: out/dataset/waterfields_comparative_2024.csv, one row per site and metric; the charts are
drawn by pipeline.maps.

The note column flags the values the deck does not read, by the same tests as the deck: radar water
where Sentinel-1 cannot see the site's water (pipeline.maps.radar_cannot_see_water), and the
cardinal share of rectilinear water or of field edges where the 30-degree rotation test fails or
has too few bodies (orientation_control in water_stats.json, rotation_test in field_stats.json).

    python -m pipeline.dataset
"""
import csv
import json

import numpy as np
import rasterio
from pyproj import Transformer

import config
from pipeline.layers import Site
from pipeline.maps import JRC_BOTH_YEARS, compare, radar_cannot_see_water
from pipeline.sources import MIN_GROUP_KM2

HEADER = ["site_key", "site", "metric", "unit", "value", "range_across_thresholds", "note"]
# Notes on values that are not read, worded as on the deck. ROTATION_NOTE is keyed on "passed" in
# orientation_control of water_stats.json; a passed test (True) needs no note.
RADAR_BLIND_NOTE = "not measurable: radar does not see water here"
ROTATION_NOTE = {False: "fails the rotation test, not interpreted", None: "too few bodies to test, not interpreted"}
# Where the paddy radar signature covers under this share of the mapped farmland, the farmland map misses
# the small sawah plots (Prambanan), so its area is a lower bound (the rule the deck uses too).
PADDY_BLIND_SHARE = 0.01
LOWER_BOUND_NOTE = ("lower bound: 10 m land cover labels many small sawah and home-garden plots as tree cover or "
                    "built-up land, and the paddy radar signature fails here")
# Reporting only, no class or total changes: water bodies whose centroid lies above this height (Copernicus
# DEM) are counted apart, because on volcano slopes (Prambanan) dark surfaces there may be shadow or
# volcanic deposit rather than water.
HIGH_GROUND_M = 500
# The water-source groups of pipeline.sources, as named in the CSV, and the rule behind each.
SOURCE_GROUPS = {"reservoir reach only": "reservoir reach only", "natural water only": "natural water only",
                 "both": "both reservoir reach and natural water", "neither": "neither"}
SOURCE_NOTES = {
    "reservoir reach only": "within 5 km of, and at least 1 m below, a rectilinear reservoir or moat segment of 10 ha or more",
    "natural water only": "within 1 km of a river, stream, lake or irregular reservoir, or flooded 1 month or more in 2024",
    "both": "",
    "neither": "rain, groundwater or works too small to see",
}
HIGH_GROUND_NOTE = ("on high ground; dark surfaces there may be shadow or volcanic deposit rather than water; "
                    "included in the water totals")


def load(site, name):
    """A site's products/<name>.json, or {} when that step has not run."""
    p = site.products / f"{name}.json"
    return json.loads(p.read_text()) if p.exists() else {}


def get(d, *path):
    """Nested value d[k1][k2]..., or None where any level is missing."""
    for k in path:
        d = d.get(k) if isinstance(d, dict) else None
    return d


def range_text(r):
    """A [low, high] pair as "low to high"; blank when there is no range."""
    return "" if not r else f"{r[0]} to {r[1]}"


def f1(a, section, name):
    """F1 score of one comparison in accuracy_stats.json."""
    return get(a, section, name, "f1")


def add_note(note, flag):
    """A note with a flag added after it ("; " between them); the note alone when flag is empty."""
    return "; ".join(part for part in (note, flag) if part)


def count_large_ponds(site):
    """Rectilinear and irregular ponds of 1 to 10 ha, or None before pipeline.water has run."""
    fpath = site.products / "water_features.csv"
    if not fpath.exists():
        return None
    n = 0
    with open(fpath) as fh:
        for r in csv.DictReader(fh):
            if 1 <= float(r["area_ha"]) < 10 and r["class"] in ("Rectilinear pond", "Irregular pond"):
                n += 1
    return n


def high_ground_water(site):
    """Count and area (km2) of the water bodies whose centroid lies above HIGH_GROUND_M, or (None, None)."""
    fpath, dem = site.products / "water_features.csv", site.static / "dem.tif"
    if not (fpath.exists() and dem.exists()):
        return None, None
    with open(fpath) as fh:
        bodies = list(csv.DictReader(fh))
    if not bodies:
        return 0, 0.0
    with rasterio.open(dem) as ds:
        to_grid = Transformer.from_crs(4326, ds.crs, always_xy=True)
        xs, ys = to_grid.transform([float(b["lon"]) for b in bodies], [float(b["lat"]) for b in bodies])
        rows, cols = rasterio.transform.rowcol(ds.transform, xs, ys)
        height = ds.read(1)
    high = [b for b, r, c in zip(bodies, rows, cols)
            if 0 <= r < height.shape[0] and 0 <= c < height.shape[1] and height[r, c] > HIGH_GROUND_M]
    return len(high), round(sum(float(b["area_ha"]) for b in high) / 100, 2)


def site_rows(site, w):
    """The study area, and the dry and wet season windows of config.py (start and end, inclusive)."""
    rows = [("Study area", "km2", w.get("study_area_km2"), "", site.meta["aoi_source"])]
    for season in ("dry", "wet"):
        start, end = site.meta["seasons"][season]
        rows.append((f"{season.capitalize()} season window, start", "date", start, "",
                     "Sentinel-2 seasonal composites use these dates; the monthly radar and NDVI series use "
                     "every calendar month the window touches"))
        rows.append((f"{season.capitalize()} season window, end", "date", end, "", "inclusive"))
    return rows


def water_rows(site, w, c, blind):
    """Open water in both seasons, and the water bodies; blind flags the radar rows as not measurable.

    Expansion and contraction are read on the radar basis alone (Sentinel-1 dry-season water against
    the wet-season maximum), never against the optical dry-season water.
    """
    ws = w.get("sensitivity", {})
    wk = c.get("water_km2", {})
    radar = RADAR_BLIND_NOTE if blind else ""
    return [
        ("Open water, dry season", "km2", w.get("dry_season_water_km2"), range_text(ws.get("dry_season_water_km2_range")),
         f"Sentinel-2 NDWI {get(w, 'threshold_method', 'ndwi') or ''}"),
        ("Open water, dry season", "% of area", w.get("dry_season_water_share_pct"), "", ""),
        ("Radar water, dry season", "km2", get(c, "radar_water_change_km2", "dry_season_water"), "",
         add_note("Sentinel-1 VV, the radar basis for expansion and contraction", radar)),
        ("Radar water in any wet-season month", "km2", w.get("s1_wet_max_water_km2"),
         range_text(ws.get("s1_wet_max_water_km2_range")),
         add_note(f"Sentinel-1 VV {get(w, 'threshold_method', 's1_vv') or ''}; water in at least one month of the "
                  "wet season, not an extent reached at one time", radar)),
        ("Wet-season expansion (radar water in the wet season only)", "km2", w.get("seasonal_inundation_km2"),
         range_text(ws.get("seasonal_inundation_km2_range")), add_note("expansion of water", radar)),
        ("Seasonal inundation (radar water 1 to 9 months of 2024)", "km2", w.get("s1_intermittent_water_km2"), "",
         add_note("Sentinel-1 VV monthly medians", radar)),
        ("Water in both seasons", "km2", wk.get("Permanent water (both seasons)"), "", add_note("Sentinel-1", radar)),
        ("Water in the dry season only", "km2", wk.get("Dry-season only water"), "",
         add_note("contraction in the wet season", radar)),
        ("Permanent water (10 or more months)", "km2", w.get("s1_permanent_water_km2"), "", radar),
        ("Water bodies", "count", w.get("feature_count"), "", "0.1 ha minimum"),
        ("Water bodies", "per 100 km2", w.get("features_per_100km2"), "", ""),
        ("Ponds and small water bodies", "per 100 km2", w.get("ponds_per_100km2"), "",
         "rectilinear ponds, irregular ponds and small water bodies, all under 10 ha"),
        ("Large ponds (1 to 10 ha)", "count", count_large_ponds(site), "", "rectilinear and irregular ponds"),
        *high_ground_rows(site),
    ]


def high_ground_rows(site):
    """Water bodies above HIGH_GROUND_M: their count and area, with a caveat where there are any."""
    n, km2 = high_ground_water(site)
    note = HIGH_GROUND_NOTE if n else ""
    return [(f"Water bodies above {HIGH_GROUND_M} m", "count", n, "", note),
            (f"Water bodies above {HIGH_GROUND_M} m", "km2", km2, "", note)]


def water_class_rows(w):
    """Count and area of each water-body class, and the share of rectilinear bodies that are cardinal.

    The cardinal share is flagged as not interpreted unless it passes the rotation test of pipeline/water.py.
    """
    rows = []
    for name, n in (w.get("count_by_class") or {}).items():
        rows.append((f"Water class: {name}", "count", n, "", ""))
        rows.append((f"Water class: {name}", "km2", get(w, "area_km2_by_class", name), "", ""))
    passed = get(w, "orientation_control", "passed")
    rows.append(("Rectilinear water bodies within 7.5 deg of cardinal", "%", w.get("cardinal_share_of_rectilinear_pct"), "",
                 add_note("16.7% if random", ROTATION_NOTE.get(passed, ""))))
    return rows


def jrc_both_years(site):
    """Area (km2) of the water present in both 1984 and 2020 that was permanent in both years and that was
    seasonal in at least one of them (pipeline.maps.JRC_BOTH_YEARS), or {} without the JRC transitions."""
    if not (site.static / "jrc_transitions.tif").exists():
        return {}
    tr, inside = site.static10("jrc_transitions"), site.mask(10)
    return {name: round(float((np.isin(tr, codes) & inside).sum() * site.g10.pixel_area_km2), 2)
            for name, codes in JRC_BOTH_YEARS.items()}


def long_term_rows(site, c):
    """Long-term water history from the JRC transitions, with the water of both years split into permanent
    and seasonal water."""
    lt = (c.get("long_term") or {}).get("area_km2", {})
    split = jrc_both_years(site)
    return [
        ("Water since 1984 (JRC 1984-2020)", "km2", lt.get("water since 1984"), "",
         "JRC transitions; water in both 1984 and 2020"),
        ("Water permanent in 1984 and 2020 (JRC 1984-2020)", "km2", split.get("permanent 1984 and 2020"), "",
         "JRC class permanent"),
        ("Water seasonal in 1984 and 2020 (JRC 1984-2020)", "km2", split.get("seasonal 1984 and 2020"), "",
         "JRC classes seasonal, seasonal to permanent and permanent to seasonal"),
        ("Water new since 1984 (JRC 1984-2020)", "km2", lt.get("new since 1984"), "",
         "JRC transitions; water in 2020 not recorded in 1984, not newly built"),
        ("Water lost since 1984 (JRC 1984-2020)", "km2", lt.get("lost since 1984"), "", "JRC transitions"),
        ("Ephemeral water (JRC 1984-2020)", "km2", lt.get("ephemeral"), "", "JRC transitions"),
    ]


def farmland_rows(f):
    """Agricultural land, its classes and how often it is cropped."""
    fs = f.get("sensitivity", {})
    cls = f.get("area_km2_by_class", {})
    cyc = f.get("cropping_cycles_share_pct", {})
    cropped_twice = (cyc.get("2") or 0) + (cyc.get("3") or 0) if cyc else None  # two or three cycles a year
    paddy, ag = f.get("paddy_signature_km2"), f.get("agricultural_km2")
    lower = LOWER_BOUND_NOTE if paddy is not None and ag and paddy < PADDY_BLIND_SHARE * ag else ""
    return [
        ("Agricultural land", "km2", ag, "", lower),
        ("Agricultural land", "% of area", f.get("agricultural_share_pct"), "", lower),
        ("Irrigated or multi-cropped", "km2", cls.get("Irrigated or multi-cropped"), "", ""),
        ("Rainfed rice-field landscape", "km2", cls.get("Rice-field landscape (rainfed)"), "", ""),
        ("Paddy radar signature, any class", "km2", f.get("paddy_signature_km2"), range_text(fs.get("paddy_signature_km2_range")),
         "highly threshold-dependent"),
        ("Flood-recession fields", "km2", cls.get("Flood-recession fields"), "", ""),
        ("Open agricultural land", "km2", cls.get("Open agricultural land"), "", ""),
        ("Farmland cropped twice or more", "%", cropped_twice, "",
         "probably a lower bound: monthly optical NDVI misses some second crops"),
        ("Concentrated agriculture (80% or more of a 1 km window)", "km2", f.get("concentrated_ag_km2"), "", ""),
    ]


def seasonal_change_rows(c):
    """Vegetation change between the dry and the wet season, over the whole area and on farmland alone."""
    farmland = "on agricultural land, where both seasons are clear"
    share_note = "share of all agricultural land; a lower bound where cloud leaves gaps"
    return [
        ("Wet-season greening (NDVI change more than +0.15)", "km2", c.get("ndvi_greening_km2"), "",
         "where both seasons are clear"),
        ("Wet-season browning (NDVI change less than -0.15)", "km2", c.get("ndvi_browning_km2"), "",
         "where both seasons are clear"),
        ("Farmland greening in the wet season (NDVI change more than +0.15)", "km2", c.get("ag_ndvi_greening_km2"), "",
         farmland),
        ("Farmland greening in the wet season (NDVI change more than +0.15)", "% of farmland",
         c.get("ag_ndvi_greening_share_pct"), "", share_note),
        ("Farmland browning in the wet season (NDVI change less than -0.15)", "km2", c.get("ag_ndvi_browning_km2"), "",
         farmland),
        ("Farmland browning in the wet season (NDVI change less than -0.15)", "% of farmland",
         c.get("ag_ndvi_browning_share_pct"), "", share_note),
        ("Clear view in both seasons, farmland", "% of farmland", c.get("ag_ndvi_change_coverage_pct"), "",
         "dry and wet season both clear"),
        ("Median NDVI change on farmland, wet minus dry", "NDVI", c.get("ag_ndvi_change_median"), "", ""),
        ("Farmland with radar VH rise over 3 dB", "km2", c.get("ag_vh_rise_over_3db_km2"), "", "cloud-free companion"),
        ("Clear view in both seasons", "% of area", c.get("ndvi_change_coverage_pct"), "", "dry and wet season both clear"),
    ]


def passes_rotation_test(fo):
    """True when the field-edge orientation passes the 30-degree image rotation test of pipeline/fields.py."""
    rot = (fo.get("controls") or {}).get("rotation_test") or {}
    return bool(rot) and rot["after_30deg_rotation_image_axes_pct"] < rot["window_share_pct"]


def orientation_rows(f):
    """Cardinal share of strong field edges (weighted by edge strength), and the same outside farmland."""
    fo = f.get("field_orientation", {})
    note = "range is a 95% block-bootstrap interval; 16.7% if random"
    if not passes_rotation_test(fo):
        note += "; fails the rotation test, not interpreted"
    return [
        ("Strong field edges within 7.5 deg of cardinal, weighted by edge strength", "%", fo.get("cardinal_share_pct"),
         range_text(fo.get("cardinal_share_ci95_pct")), note),
        ("Same, edges outside farmland", "%", get(fo, "controls", "non_agricultural_edges_pct"),
         range_text(get(fo, "controls", "non_agricultural_ci95_pct")), "control"),
    ]


def ftw_rows(t):
    """Fields of The World: field count and size, overlap with the farmland map and field orientation."""
    return [
        ("Fields of The World fields", "count", t.get("fields"), "", "2024, inside the polygon"),
        ("Fields of The World fields", "per km2", t.get("fields_per_km2"), "", ""),
        ("Fields of The World median field size", "ha", t.get("median_field_ha"), "",
         "near the resolution limit: 0.1 ha is 10 pixels at 10 m"),
        ("Mapped farmland covered by Fields of The World fields", "%", t.get("share_of_our_agriculture_covered_pct"),
         "", ""),
        ("Fields of The World field area inside the mapped farmland", "%",
         t.get("share_of_ftw_inside_our_agriculture_pct"), "", ""),
        ("Fields of The World fields within 7.5 deg of cardinal", "%", get(t, "orientation", "cardinal_share_pct"), "",
         f"descriptive only; {get(t, 'orientation', 'fields_used')} fields of 2 ha or more used; 16.7% if random"),
    ]


def group_value(src, name, field):
    """A value of one water-source group; blank when the group holds almost no farmland (under 1 km2)."""
    g = get(src, "groups", name) or {}
    return g.get(field) if (g.get("farmland_km2") or 0) >= MIN_GROUP_KM2 else None


def source_rows(src):
    """Irrigated share of farmland by the water source within reach (pipeline.sources), and the sweep.

    Reservoir reach: within 5 km of, and at least 1 m below, a rectilinear reservoir or moat segment of
    10 ha or more. Natural water: within 1 km of a river, stream, lake or irregular reservoir, or flooded
    for a month or more in 2024. The ranges are over the settings of the sweep in which the group holds
    1 km2 or more of farmland (all 27, or fewer).
    """
    if not src:
        return []
    sweep = src.get("sensitivity", {})
    n = sweep.get("settings_total")
    no_reservoir = "no farmland within reach of a rectilinear reservoir of 10 ha or more"
    rows = []
    for name, label in SOURCE_GROUPS.items():
        g = get(src, "groups", name) or {}
        needs_reservoir = name in ("reservoir reach only", "both")
        empty = no_reservoir if needs_reservoir and (g.get("farmland_km2") or 0) < MIN_GROUP_KM2 else ""
        swept = get(sweep, "irrigated_vs_all_range", name) is not None  # compared in the sweep, with a range
        n_group = get(sweep, "settings_with_group", name)  # settings in which the group holds 1 km2 or more
        rows += [
            (f"Farmland, {label}", "km2", g.get("farmland_km2"), "", add_note(SOURCE_NOTES[name], empty)),
            (f"Farmland, {label}", "% of farmland", g.get("farmland_share_pct"), "", empty),
            (f"Irrigated or multi-cropped share, {label}", "%", group_value(src, name, "irrigated_share_pct"),
             range_text(get(sweep, "irrigated_share_pct_range", name)) if swept else "",
             add_note(f"2024 pattern only; range over the {n_group} of {n} settings of the sweep in which the group "
                      "holds 1 km2 or more of farmland" if swept else "2024 pattern only",
                      empty)),
            (f"Irrigated or multi-cropped share, {label}, ratio to all farmland", "ratio",
             group_value(src, name, "irrigated_vs_all"),
             range_text(get(sweep, "irrigated_vs_all_range", name)) if swept else "", empty),
        ]
    compared = sweep.get("settings_compared")
    rows += [
        ("Irrigated or multi-cropped share, all farmland", "%", src.get("irrigated_share_all_farmland_pct"), "",
         "the reference for the water-source groups"),
        ("Water-source sweep settings compared", "count", compared, "",
         f"of {n} settings; both groups hold 1 km2 or more of farmland" if compared else no_reservoir),
        ("Water-source sweep settings with reservoir reach above natural water", "count",
         sweep.get("settings_with_reservoir_above_natural") if compared else None, "",
         f"of {compared} settings compared; a 2024 pattern, not what past reservoirs did" if compared
         else no_reservoir),
    ]
    return rows


def agreement_rows(a, blind):
    """Agreement (F1) with the reference maps of pipeline.accuracy; blind flags the radar water rows."""
    radar = RADAR_BLIND_NOTE if blind else ""
    return [
        ("F1 dry water vs JRC 1984-2020", "0-1", f1(a, "dry_season_water", "JRC occurrence 50% or more (1984-2020, 30 m)"), "",
         "agreement between maps, not accuracy"),
        ("F1 dry water vs ESA WorldCover 2021", "0-1", f1(a, "dry_season_water", "ESA WorldCover 2021 water"), "", ""),
        ("F1 wet-season expansion vs JRC seasonal water (2020)", "0-1",
         f1(a, "seasonal_inundation", "JRC seasonal water (1 to 11 months in 2020)"), "",
         add_note("2024 against 2020", radar)),
        ("F1 agriculture vs ESA WorldCover 2021 cropland", "0-1",
         f1(a, "agricultural_domain", "ESA WorldCover 2021 cropland"), "",
         "not independent: the agricultural domain starts from WorldCover cropland"),
        ("F1 agriculture vs Fields of The World", "0-1", f1(a, "agricultural_domain", "Fields of The World 2024 fields"), "", ""),
    ]


def rows_for(key):
    """Every metric row of one site, in the order of the CSV."""
    s = Site(key)
    w, f, c = load(s, "water_stats"), load(s, "field_stats"), load(s, "change_stats")
    t, a, src = load(s, "ftw_stats"), load(s, "accuracy_stats"), load(s, "water_sources")
    blind = radar_cannot_see_water(s)  # Sentinel-1 cannot see this site's water (in 2024, Prambanan only)
    rows = (site_rows(s, w) + water_rows(s, w, c, blind) + water_class_rows(w) + long_term_rows(s, c)
            + farmland_rows(f) + seasonal_change_rows(c) + orientation_rows(f) + ftw_rows(t)
            + source_rows(src) + agreement_rows(a, blind))
    return [(key, s.meta["name"], *r) for r in rows]  # each row starts with the site key and name


def run():
    """Write the comparative CSV for every site, then draw the comparison charts."""
    out = config.OUT / "dataset"
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"waterfields_comparative_{config.YEAR}.csv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(HEADER)
        for key in config.SITES:
            w.writerows(rows_for(key))
    compare()
    print("comparative dataset:", path, flush=True)


if __name__ == "__main__":
    run()
