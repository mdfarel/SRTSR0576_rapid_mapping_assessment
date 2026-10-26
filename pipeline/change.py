"""Wet-season versus dry-season change.

Water (Sentinel-1, cloud-free): permanent water present in both seasons, water
present only in the wet season (seasonal inundation), and water present only in
the dry season (storage filled for dry-season irrigation, or reservoirs that
spill into vegetation in the wet season). Expansion and contraction are both
read on this one radar basis, with the same VV threshold in both seasons: the
dry-season median against any wet-season month. They are reported together
under radar_water_change_km2 in change_stats.json; the optical dry-season water
of pipeline.water is a different measure and is not compared with them.

Vegetation: Sentinel-2 NDVI wet minus dry, where both seasons have at least one
clear observation, over the whole study area and over agricultural land alone,
since cultivation drives most of the seasonal change there. Wet-season cloud
leaves gaps, so the radar cross-polarised backscatter (VH, sensitive to
vegetation volume) is given as a gap-free companion.

Inputs: the seasonal water maps from pipeline.water, the field classes from
pipeline.fields, the Sentinel-2 seasonal composites, the Sentinel-1 monthly stack
and, when present, the JRC Global Surface Water transitions.
Outputs (data/<site>/products): water_change_20m, ndvi_change_20m and
vh_change_20m (GeoTIFFs) and change_stats.json.

    python -m pipeline.change [site ...]
"""
import sys

import numpy as np

import config
from pipeline.layers import Site, ndi
from pipeline.water import CLASSES as WATER_CLASSES

WATER = {1: "Permanent water (both seasons)", 2: "Seasonal inundation (wet season only)",
         3: "Dry-season only water"}
NDVI_CHANGE = 0.15   # wet minus dry NDVI beyond which a pixel counts as greening or browning
VH_RISE_DB = 3       # wet minus dry VH rise reported over agricultural land
MIN_CLASS_PX = 100   # 1 ha at 10 m: a water-body class covering less in total gets no long-term breakdown
RADAR_BASIS = "Sentinel-1 VV, one threshold for both seasons: the dry-season median against any wet-season month"

# JRC Global Surface Water transition classes, first year (1984) to last year (2020), Pekel et al. 2016:
# 1 permanent, 2 new permanent, 3 lost permanent, 4 seasonal, 5 new seasonal, 6 lost seasonal,
# 7 seasonal to permanent, 8 permanent to seasonal, 9 ephemeral permanent, 10 ephemeral seasonal
GROUPS = {"water since 1984": (1, 4, 7, 8), "new since 1984": (2, 5), "lost since 1984": (3, 6),
          "ephemeral": (9, 10)}


def trajectories(site):
    """Long-term water history (JRC, 1984 to 2020) of the study area and of each mapped water-body class."""
    path = site.static / "jrc_transitions.tif"
    if not path.exists():
        return None
    tr = site.static10("jrc_transitions")
    inside = site.mask(10)
    px = site.g10.pixel_area_km2
    out = {"area_km2": {n: round(float((np.isin(tr, codes) & inside).sum() * px), 2) for n, codes in GROUPS.items()}}
    cls = site.load("water_class_10m")
    by = {}
    for k, name in WATER_CLASSES.items():
        m = (cls == k) & inside
        tot = int(m.sum())
        if tot < MIN_CLASS_PX:
            continue
        by[name] = {g: round(100 * float((np.isin(tr, codes) & m).sum()) / tot, 1) for g, codes in GROUPS.items()}
        # the rest of the class has no JRC water at all
        by[name]["no water recorded 1984-2020"] = round(100 - sum(by[name].values()), 1)
    out["share_of_mapped_class_pct"] = by
    return out


def radar_water_change(water_change, px):
    """Expansion and contraction of water on the Sentinel-1 basis, in km2, from the water change map."""
    both, wet_only, dry_only = ((water_change == k).sum() for k in (1, 2, 3))
    return {
        "dry_season_water": round(float((both + dry_only) * px), 2),
        "wet_season_maximum_water": round(float((both + wet_only) * px), 2),
        "water_in_both_seasons": round(float(both * px), 2),
        "expansion_wet_season_only": round(float(wet_only * px), 2),
        "contraction_dry_season_only": round(float(dry_only * px), 2),
        "net_change_wet_minus_dry": round(float((wet_only - dry_only) * px), 2),
    }


def farmland_ndvi_change(dndvi, ag, px):
    """Greening and browning on agricultural land, in km2 and as a share of all agricultural land."""
    n_ag = max(int(ag.sum()), 1)
    greening = (dndvi > NDVI_CHANGE) & ag
    browning = (dndvi < -NDVI_CHANGE) & ag
    return {
        "ag_ndvi_change_coverage_pct": round(100 * float((np.isfinite(dndvi) & ag).sum()) / n_ag, 1),
        "ag_ndvi_greening_km2": round(float(greening.sum() * px), 1),
        "ag_ndvi_browning_km2": round(float(browning.sum() * px), 1),
        "ag_ndvi_greening_share_pct": round(100 * float(greening.sum()) / n_ag, 1),
        "ag_ndvi_browning_share_pct": round(100 * float(browning.sum()) / n_ag, 1),
    }


def run(key):
    """Seasonal change maps and change_stats.json for one site."""
    site = Site(key)
    inside = site.mask(20)

    # water: 1 both seasons, 2 wet season only, 3 dry season only (see WATER)
    dry_w = (site.load("s1_water_dry_20m") == 1) & inside
    wet_w = (site.load("s1_water_wet_20m") == 1) & inside
    water_change = np.zeros(dry_w.shape, np.uint8)
    water_change[dry_w & wet_w] = 1
    water_change[wet_w & ~dry_w] = 2
    water_change[dry_w & ~wet_w] = 3
    site.save("water_change_20m", water_change, site.g20, nodata=0, dtype="uint8")

    # NDVI wet minus dry, only where both seasons have a clear view
    dry, wet = site.s2_composite("dry", 20), site.s2_composite("wet", 20)
    ok = (site.s2_count("dry", 20) >= 1) & (site.s2_count("wet", 20) >= 1) & inside
    dndvi = np.where(ok, ndi(wet["B08"], wet["B04"]) - ndi(dry["B08"], dry["B04"]), np.nan).astype(np.float32)
    site.save("ndvi_change_20m", dndvi, site.g20, nodata=np.nan, dtype="float32")

    # radar VH wet minus dry (seasonal medians), left blank over radar water in either season
    _, vh = site.s1_monthly()
    dry_idx = [m - 1 for m in site.season_months("dry")]  # months 1-12 to array index 0-11
    wet_idx = [m - 1 for m in site.season_months("wet")]
    with np.errstate(all="ignore"):
        vh_d = np.nanmedian(vh[dry_idx], 0)
        vh_w = np.nanmedian(vh[wet_idx], 0)
    dvh = (vh_w - vh_d).astype(np.float32)
    dvh[(water_change > 0) | ~inside] = np.nan
    site.save("vh_change_20m", dvh, site.g20, nodata=np.nan, dtype="float32")

    fields = site.load("field_class_20m")
    ag = fields > 0
    px = site.g20.pixel_area_km2
    valid = np.isfinite(dndvi)
    stats = {
        "site": site.meta["name"],
        "seasons": site.meta["seasons"],
        "water_km2": {n: round(float((water_change == k).sum() * px), 2) for k, n in WATER.items()},
        "radar_water_change_km2": radar_water_change(water_change, px),
        "radar_water_change_basis": RADAR_BASIS,
        "ndvi_change_coverage_pct": round(float(valid.sum() / inside.sum() * 100), 1),
        "ndvi_greening_km2": round(float((dndvi > NDVI_CHANGE).sum() * px), 1),
        "ndvi_browning_km2": round(float((dndvi < -NDVI_CHANGE).sum() * px), 1),
        "ag_ndvi_change_median": round(float(np.nanmedian(dndvi[ag & valid])), 3) if (ag & valid).any() else None,
        **farmland_ndvi_change(dndvi, ag, px),
        "ag_vh_change_median_db": round(float(np.nanmedian(dvh[ag])), 2) if ag.any() else None,
        "ag_vh_rise_over_3db_km2": round(float(((dvh > VH_RISE_DB) & ag).sum() * px), 1),
        "s2_wet_clear_obs_median": int(np.median(site.s2_count("wet", 20)[inside])),
    }
    stats["long_term"] = trajectories(site)
    site.save_json("change_stats", stats)
    print(key, stats["water_km2"], "NDVI coverage", stats["ndvi_change_coverage_pct"])


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
