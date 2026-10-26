"""Which water sources sit next to today's irrigated and flood-recession land?

A present-day pattern for framing the long debate on whether large reservoirs
fed rice fields (van Liere 1980), not a test of what past reservoirs did:
modern canals, pumps and wells are invisible here, and 2024 is one year.

Every agricultural pixel (20 m) is placed in one of four groups:
- reservoir reach: within 5 km of a rectilinear reservoir or moat segment of
  10 ha or more (pipeline.water) and at least 1 m below that body's surface on
  the Copernicus DEM, so that gravity could carry water there;
- natural water: within 1 km of a river or stream, or of a lake or irregular
  reservoir, or flooded for at least one month in 2024 (Sentinel-1);
- both;
- neither (rain, groundwater or works too small to see).
For each group the share of its farmland that is irrigated or multi-cropped,
and the share that is flood-recession, is reported against the share over all
farmland (pipeline.fields). The three distances are then varied together over
27 settings (reach 3, 5 or 8 km; natural water 0.5, 1 or 2 km; drop 0, 1 or
3 m) to show how far the comparison depends on them. Nothing here changes a
water or field class.

Inputs: field_class_20m.tif, water_class_10m.tif and s1_water_months_20m.tif
from the products, and the Copernicus DEM.
Output: data/<site>/products/water_sources.json.

    python -m pipeline.sources [site ...]
"""
import sys

import numpy as np
from scipy import ndimage

import config
from pipeline.fetch_static import DEM_NODATA
from pipeline.fields import CLASSES as FIELD_CLASSES
from pipeline.layers import Site
from pipeline.water import CLASSES as WATER_CLASSES

REACH_KM = 5.0    # how far a reservoir could plausibly feed fields by gravity
NATURAL_KM = 1.0  # how near a river or lake counts as "next to" it
DROP_M = 1.0      # how far below the reservoir's surface a field must lie
# the sensitivity sweep: every combination of these, 27 settings in all
REACH_SWEEP_KM = (3.0, 5.0, 8.0)
NATURAL_SWEEP_KM = (0.5, 1.0, 2.0)
DROP_SWEEP_M = (0.0, 1.0, 3.0)
MIN_GROUP_KM2 = 1.0  # a group with less farmland than this is not compared in the sweep


def codes(classes, *labels):
    """Class codes of the given labels in a CLASSES dictionary (so the codes follow pipeline.water and fields)."""
    found = tuple(k for k, v in classes.items() if v in labels)
    assert len(found) == len(labels), f"class labels not found: {labels}"
    return found


RESERVOIR_CLASSES = codes(WATER_CLASSES, "Rectilinear reservoir or moat segment (10 ha or more)")
NATURAL_CLASSES = codes(WATER_CLASSES, "River or stream", "Lake or irregular reservoir")
RECESSION = codes(FIELD_CLASSES, "Flood-recession fields")[0]
IRRIGATED = codes(FIELD_CLASSES, "Irrigated or multi-cropped")[0]
GROUP_NAMES = ("reservoir reach only", "natural water only", "both", "neither")


def classes_any_20m(cls10, wanted):
    """True for each 20 m cell where any of its four 10 m pixels has one of the given water classes."""
    h, w = cls10.shape
    b = np.isin(cls10[: h // 2 * 2, : w // 2 * 2], wanted).reshape(h // 2, 2, w // 2, 2)
    return b.any(axis=(1, 3))


def distances(site):
    """Farmland, flooding and the distance (and height) of every 20 m cell relative to reservoirs and natural water."""
    inside = site.mask(20)
    fields = site.load("field_class_20m")
    ag = (fields > 0) & inside
    cls10 = site.load("water_class_10m")
    reservoir = classes_any_20m(cls10, RESERVOIR_CLASSES) & inside
    natural = classes_any_20m(cls10, NATURAL_CLASSES) & inside
    months = site.load("s1_water_months_20m")
    flooded = (months >= 1) & (months <= 12) & inside  # 255 is no data
    dem = site.static10("dem")[::2, ::2].astype(np.float32)
    dem[dem <= DEM_NODATA + 1] = np.nan  # no height: never counted as below a reservoir
    km_per_px = site.g20.res / 1000.0

    if reservoir.any():
        # distance to the nearest reservoir cell, and how far the field lies below that cell
        dist, (ri, ci) = ndimage.distance_transform_edt(~reservoir, return_indices=True)
        with np.errstate(invalid="ignore"):
            dist_res, below = dist * km_per_px, np.nan_to_num(dem[ri, ci] - dem, nan=-np.inf)
    else:
        dist_res, below = np.full(ag.shape, np.inf), np.zeros(ag.shape, np.float32)
    if natural.any():
        dist_nat = ndimage.distance_transform_edt(~natural) * km_per_px
    else:
        dist_nat = np.full(ag.shape, np.inf)
    return {"fields": fields, "ag": ag, "flooded": flooded, "dist_res": dist_res, "below": below,
            "dist_nat": dist_nat, "px": site.g20.pixel_area_km2, "total": int(ag.sum())}


def share(fields, cls, where, n):
    """Share of the n farmland cells in `where` that have field class cls."""
    return float(((fields == cls) & where).sum()) / max(n, 1)


def groups(d, irr_all, reach_km, natural_km, drop_m):
    """Farmland and its irrigated and recession shares in each of the four water-source groups."""
    reach = (d["dist_res"] <= reach_km) & (d["below"] >= drop_m)
    nat = (d["dist_nat"] <= natural_km) | d["flooded"]
    masks = dict(zip(GROUP_NAMES, (reach & ~nat, nat & ~reach, reach & nat, ~reach & ~nat)))
    out = {}
    for name, m in masks.items():
        g = d["ag"] & m
        n = int(g.sum())
        irr = share(d["fields"], IRRIGATED, g, n)
        rec = share(d["fields"], RECESSION, g, n)
        out[name] = {"farmland_km2": round(n * d["px"], 1), "farmland_share_pct": round(100 * n / max(d["total"], 1), 1),
                     "irrigated_share_pct": round(100 * irr, 1),
                     "irrigated_vs_all": round(irr / irr_all, 2) if irr_all > 0 and n else None,
                     "recession_share_pct": round(100 * rec, 1)}
    return out


def sensitivity(d, irr_all):
    """Irrigated share, and its ratio to the all-farmland share, over every setting of the sweep."""
    compared = ("reservoir reach only", "natural water only")
    ratios = {name: [] for name in compared}
    shares = {name: [] for name in compared}
    paired = []
    for rk in REACH_SWEEP_KM:
        for nk in NATURAL_SWEEP_KM:
            for dm in DROP_SWEEP_M:
                g = groups(d, irr_all, rk, nk, dm)
                ok = {}  # ratios for the groups that hold at least 1 km2 of farmland under this setting
                for name in compared:
                    if g[name]["irrigated_vs_all"] is not None and g[name]["farmland_km2"] >= MIN_GROUP_KM2:
                        ok[name] = g[name]["irrigated_vs_all"]
                        shares[name].append(g[name]["irrigated_share_pct"])
                for name, v in ok.items():
                    ratios[name].append(v)
                if len(ok) == 2:
                    paired.append(ok["reservoir reach only"] > ok["natural water only"])
    n_settings = len(REACH_SWEEP_KM) * len(NATURAL_SWEEP_KM) * len(DROP_SWEEP_M)
    return {"settings": "reach 3, 5, 8 km; natural water 0.5, 1, 2 km; drop 0, 1, 3 m",
            "settings_total": n_settings,
            "irrigated_vs_all_range": {k: ([min(v), max(v)] if v else None) for k, v in ratios.items()},
            "irrigated_share_pct_range": {k: ([min(v), max(v)] if v else None) for k, v in shares.items()},
            "settings_with_group": {k: len(v) for k, v in ratios.items()},
            "settings_compared": len(paired),
            "settings_with_reservoir_above_natural": sum(paired)}


def run(key):
    """Water-source groups of the farmland, with the sensitivity sweep, for one site."""
    site = Site(key)
    d = distances(site)
    irr_all = share(d["fields"], IRRIGATED, d["ag"], d["total"])
    rec_all = share(d["fields"], RECESSION, d["ag"], d["total"])
    out = {"rule": {"reach_km": REACH_KM, "below_reservoir_m": DROP_M, "natural_km": NATURAL_KM,
                    "reservoirs": "rectilinear reservoirs or moat segments of 10 ha or more (dry season 2024)",
                    "natural": "rivers, lakes and irregular reservoirs within 1 km, or Sentinel-1 flooding in 2024",
                    "min_group_km2_in_sweep": MIN_GROUP_KM2},
           "agricultural_km2": round(d["total"] * d["px"], 1),
           "irrigated_share_all_farmland_pct": round(100 * irr_all, 1),
           "recession_share_all_farmland_pct": round(100 * rec_all, 1),
           "groups": groups(d, irr_all, REACH_KM, NATURAL_KM, DROP_M),
           "sensitivity": sensitivity(d, irr_all),
           "caution": "a 2024 pattern of where irrigated land sits, not evidence of what past reservoirs did"}
    site.save_json("water_sources", out)
    print(key, {k: (v["farmland_share_pct"], v["irrigated_share_pct"]) for k, v in out["groups"].items()},
          "all irrigated", out["irrigated_share_all_farmland_pct"], flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
