"""Agreement of every result with the reference datasets named in the task brief.

None of these references is ground truth: each is another map with its own
errors, year and resolution, so the numbers are agreement between maps, reported
as precision (share of our mapped class that the reference also shows),
recall (share of the reference class that we map), F1 and IoU. A true accuracy
assessment needs the human-labelled stratified sample (pipeline.sample,
Olofsson et al. 2014); see pipeline.sample_accuracy.

References, all among the datasets of the brief
- JRC Global Surface Water v1.3 (Pekel et al. 2016): occurrence and seasonality, 1984-2020, 30 m.
- ESA WorldCover 2021 (Zanaga et al. 2022), 10 m: the water class, and the cropland class. WorldCover
  cropland is where our agricultural domain starts (pipeline/fields.py), so that comparison shows how far
  the domain departs from WorldCover and is not an independent check.
- Fields of The World 2024 field polygons (Robinson et al. 2026), which the brief suggests.
- Sentinel-2 and Sentinel-1 themselves, as cross-sensor checks of the radar water.

Inputs: the products of pipeline.water, pipeline.fields, pipeline.phenology, pipeline.change and
pipeline.ftw, the Sentinel-1 monthly stack, the wet-season Sentinel-2 composite and the static layers.
Output: products/accuracy_stats.json.

    python -m pipeline.accuracy [site ...]
"""
import sys

import numpy as np

import config
from pipeline.layers import Site, block_mean, ndi

MIN_SUPPORT = 50  # pixels on each side; below this a comparison is reported as n/a, not as 0
KM2_PER_10M_PIXEL = 1e-4

# class codes of the reference maps
WORLDCOVER_WATER = 80
WORLDCOVER_CROP = 40
# class codes of our own maps
CYCLES_NOT_JUDGED = 255   # crop_cycles_20m nodata (pipeline/phenology.py)
SEASONAL_INUNDATION = 2   # water_change_20m: water in the wet season only (pipeline/change.py)

RADAR_CELL_WATER_SHARE = 0.75  # 20 m cells at least 75% dry-season water at 10 m
PERMANENT_MONTHS = 10          # water in at least 10 months of 2024 counts as permanent
CROP_CELL_SHARE = 0.5          # a 20 m cell is WorldCover cropland when at least half of its 10 m pixels are
WORLDCOVER_CROP_NOTE = "not independent: the agricultural domain starts from WorldCover cropland"


def scores(ours, ref):
    """Precision, recall, F1 and IoU of our mask against a reference mask on the same grid."""
    if int(ours.sum()) < MIN_SUPPORT or int(ref.sum()) < MIN_SUPPORT:
        return {"precision": None, "recall": None, "f1": None, "iou": None,
                "note": f"too little support (ours {int(ours.sum())} px, reference {int(ref.sum())} px)"}
    tp = int((ours & ref).sum())
    fp = int((ours & ~ref).sum())
    fn = int((~ours & ref).sum())
    p = tp / max(tp + fp, 1)
    r = tp / max(tp + fn, 1)
    return {"precision": round(p, 3), "recall": round(r, 3), "f1": round(2 * p * r / max(p + r, 1e-9), 3),
            "iou": round(tp / max(tp + fp + fn, 1), 3)}


def pct(part, whole):
    """Percentage of the cells of `whole` that are in `part`, to one decimal."""
    return round(100 * float(part.sum()) / max(int(whole.sum()), 1), 1)


# --- the checks; each adds its own sections to `out` ---

def dry_season_water(site, water, jrc, out):
    """1. Dry-season open water (10 m) against JRC and WorldCover."""
    inside10 = site.mask(10)
    wc = site.static10("worldcover")
    refs = {
        "JRC occurrence 50% or more (1984-2020, 30 m)": (jrc >= 50) & (jrc <= 100) & inside10,  # above 100 is no data
        "ESA WorldCover 2021 water": (wc == WORLDCOVER_WATER) & inside10,
    }
    section = {}
    for name, ref in refs.items():
        section[name] = scores(water, ref) | {"reference_km2": round(float(ref.sum() * KM2_PER_10M_PIXEL), 2)}
    section["ours_km2"] = round(float(water.sum() * KM2_PER_10M_PIXEL), 2)
    out["dry_season_water"] = section


def wet_season_water(site, jrc, seas, out):
    """2. Wet-season water and seasonal inundation (20 m) against JRC and Sentinel-2."""
    inside20 = site.mask(20)
    wet = (site.load("s1_water_wet_20m") == 1) & inside20
    jrc20 = jrc[::2, ::2]
    out["wet_season_water"] = {
        "JRC ever water (occurrence above 0, 1984-2020)": scores(wet, (jrc20 > 0) & (jrc20 <= 100) & inside20),
    }
    inund = (site.load("water_change_20m") == SEASONAL_INUNDATION) & inside20
    out["seasonal_inundation"] = {
        "JRC seasonal water (1 to 11 months in 2020)": scores(inund, (seas >= 1) & (seas <= 11) & inside20),
    }
    # cross-sensor: Sentinel-2 wet-season MNDWI where the wet composite has a clear view
    w20 = site.s2_composite("wet", 20)
    clear = (site.s2_count("wet", 20) >= 1) & inside20
    s2wet = (ndi(w20["B03"], w20["B11"]) > 0) & clear
    out["wet_season_water"]["Sentinel-2 wet-season MNDWI above 0 (clear pixels only)"] = scores(wet & clear, s2wet)


def radar_over_optical_water(site, water, out):
    """Can the radar see the water the optical map finds? Median monthly VV over mostly-water 20 m cells.

    Open water is usually below about -18 dB.
    """
    vv, _ = site.s1_monthly()
    w10 = water.astype(np.float32)
    cells = (block_mean(w10, 2) >= RADAR_CELL_WATER_SHARE)[: vv.shape[1], : vv.shape[2]]
    if cells.any():
        with np.errstate(all="ignore"):
            med = np.nanmedian(vv[:, cells], axis=1)
        monthly = [None if not np.isfinite(v) else round(float(v), 1) for v in med]  # None for months without data
        out["s1_vv_over_optical_water_db"] = {"cells": int(cells.sum()), "monthly_median": monthly}


def permanent_water(site, seas, out):
    """Water in 10 or more months of 2024 (Sentinel-1) against JRC permanent water (12 months in 2020)."""
    inside20 = site.mask(20)
    months = site.load("s1_water_months_20m")
    perm = (months >= PERMANENT_MONTHS) & (months <= 12)
    out["permanent_water_vs_jrc_permanent"] = scores(perm & inside20, (seas == 12) & inside20)


def agricultural_domain(site, out):
    """3. Our agricultural domain (20 m) against WorldCover cropland and Fields of The World. Returns field classes."""
    inside20 = site.mask(20)
    fields = site.load("field_class_20m")
    ag = (fields > 0) & inside20
    wc_crop = (block_mean((site.static10("worldcover") == WORLDCOVER_CROP).astype(np.float32), 2) >= CROP_CELL_SHARE)
    out["agricultural_domain"] = {
        "ESA WorldCover 2021 cropland": scores(ag, wc_crop & inside20) | {"note": WORLDCOVER_CROP_NOTE}}
    if (site.products / "ftw_fields_20m.tif").exists():
        ftw = (site.load("ftw_fields_20m") == 1) & inside20
        out["agricultural_domain"]["Fields of The World 2024 fields"] = scores(ag, ftw)
    return fields


def crop_consistency(site, fields, seas, out):
    """4. Internal cross-sensor consistency of the crop classes (radar classes against optical crop cycles)."""
    cyc = site.load("crop_cycles_20m")
    judged = cyc != CYCLES_NOT_JUDGED
    recession = fields == 1
    irr = fields == 2
    rice = fields == 3
    out["crop_consistency"] = {
        "rainfed rice (radar) with at least one optical crop cycle": pct((cyc >= 1) & judged & rice, judged & rice),
        "irrigated or multi-cropped with two or more optical cycles": pct((cyc >= 2) & judged & irr, judged & irr),
        "flood-recession fields inside JRC seasonal water": pct(recession & (seas >= 1) & (seas <= 11), recession),
        "note": "Consistency between independent sensors, not accuracy; rice, irrigation and crop cycles need the labelled sample.",
    }


def orientation(site, out):
    """5. Field orientation: the controls already computed in pipeline/fields.py, plus Fields of The World."""
    fs = site.load_json("field_stats")["field_orientation"]
    ftw_cardinal = None
    if (site.products / "ftw_stats.json").exists():
        ftw_cardinal = site.load_json("ftw_stats")["orientation"]["cardinal_share_pct"]
    out["orientation"] = {"structure_tensor_cardinal_pct": fs["cardinal_share_pct"], "controls": fs["controls"],
                          "ftw_cardinal_pct": ftw_cardinal,
                          "uniform_expectation_pct": fs["uniform_expectation_pct"]}


def run(key):
    """Every agreement check for one site, written to accuracy_stats.json."""
    site = Site(key)
    out = {"site": site.meta["name"], "note": "Agreement with reference maps, not accuracy against ground truth."}
    water = site.load("water_dry_10m").astype(bool) & site.mask(10)  # our dry-season open water
    jrc = site.static10("jrc_occurrence")
    seas = site.static10("jrc_seasonality")[::2, ::2]  # months with water in 2020, at 20 m

    # each step adds its sections to `out` in the order they appear in accuracy_stats.json
    dry_season_water(site, water, jrc, out)
    wet_season_water(site, jrc, seas, out)
    radar_over_optical_water(site, water, out)
    permanent_water(site, seas, out)
    fields = agricultural_domain(site, out)
    crop_consistency(site, fields, seas, out)
    orientation(site, out)
    site.save_json("accuracy_stats", out)

    water_f1 = [v["f1"] for k, v in out["dry_season_water"].items() if k != "ours_km2"]
    ag_f1 = [v["f1"] for v in out["agricultural_domain"].values()]
    print(key, "water F1 vs JRC/WorldCover:", water_f1, "| ag F1 vs WorldCover/FTW:", ag_f1, flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
