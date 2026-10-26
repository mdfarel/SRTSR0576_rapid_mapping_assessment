"""Accuracy and error-adjusted areas from the human-labelled stratified sample.

Not part of the rapid assessment: no point has been labelled yet, and this module
documents how the labels would be scored. The labelling, and the choice of the
imagery it is done on, belong to the next stage of the work, not to this
assessment. The sample is stratified by the mapped classes (pipeline/sample.py)
and is to be labelled in two passes. Pass 1 gives land cover from imagery alone;
it is independent of the maps and is the basis of the accuracy figures. Pass 2
records wet-season flooding and crop cycles from the same Sentinel-1 and
Sentinel-2 series the maps were made from; it is reported per stratum as a
consistency check, not as accuracy.

Pass 1 is scored on what imagery can show: open water, farmland (paddy or other
cropland), and other land. Each mapped stratum belongs to one of those map
classes; seasonal inundation cells are mapped as neither dry-season water nor
farmland, so they count as other land. A land-cover label given after the
labeller saw that point's time series (pass 1 time later than pass 2) is not
independent and is left out; labels given while the Sentinel-2 layer (the maps'
own input) was open are counted and also scored separately. Estimates are given
only when every stratum has at least two judged points. Every figure is a
ratio of two stratified means with its variance from Stehman (2014, International
Journal of Remote Sensing 35:4923-4939), written for strata that differ from the
map classes; for strata equal to map classes it gives the estimators of Olofsson
et al. (2014, RSE 148:42-57). Intervals are 95%. Unsure labels are left out.

Inputs: a labels CSV (one row per point: id, site, cover_code, and the optional
pass1_at, pass2_at, s2_viewed, wet_season_flooding and crop_cycles fields read
below), and validation_sample.csv and validation_design.json from pipeline.sample.
Output: products/sample_accuracy.json for each site with labels.

    python -m pipeline.sample_accuracy waterfields_validation_labels.csv [site ...]
"""
import csv
import sys
from collections import defaultdict

import numpy as np

import config
from pipeline.layers import Site

STRATA = ["open water, dry season", "seasonal inundation", "flood-recession fields", "irrigated or multi-cropped",
          "rainfed rice landscape", "open agricultural land", "other land"]
MAP_CLASS = {"open water, dry season": "water", "seasonal inundation": "other", "flood-recession fields": "farmland",
             "irrigated or multi-cropped": "farmland", "rainfed rice landscape": "farmland",
             "open agricultural land": "farmland", "other land": "other"}
# land-cover codes of the labels CSV and the class each counts as: 1 open water, 2 paddy fields (bunded),
# 3 other cropland or fallow, 4 tree crops, orchard or plantation, 5 settlement, built-up or road,
# 6 forest, shrub or grassland, 7 bare ground, sand or wetland vegetation, 8 unsure (left out)
COVER = {1: "water", 2: "farmland", 3: "farmland", 4: "other", 5: "other", 6: "other", 7: "other"}  # 8 = unsure
PADDY = 2
REF_CLASSES = ["water", "farmland", "other"]
Z95 = 1.96  # normal quantile for a 95% interval


def ratio(y, x, h, W):
    """Stratified ratio estimator R = sum_h W_h mean(y_h) / sum_h W_h mean(x_h) and its 95% half-interval.

    y and x hold one value per judged point, h its stratum, and W the area share of each stratum.
    """
    Y = X = 0.0
    parts = []
    for s, w in W.items():
        m = h == s
        if m.sum() < 2:
            continue
        Y += w * y[m].mean()
        X += w * x[m].mean()
        parts.append((w, y[m], x[m]))
    if X <= 0:
        return None, None
    R = Y / X
    terms = []  # variance contribution of each stratum (Stehman 2014)
    for w, ys, xs in parts:
        terms.append(w * w * (np.var(ys, ddof=1) + R * R * np.var(xs, ddof=1) - 2 * R * np.cov(ys, xs, ddof=1)[0, 1])
                     / len(ys))
    v = sum(terms) / X ** 2
    return float(R), Z95 * float(np.sqrt(max(v, 0)))  # unrounded; rounded only on output


def rounded(pair, k=3):
    """An (estimate, half-interval) pair rounded to k decimals, or [None, None]."""
    return [None, None] if pair[0] is None else [round(pair[0], k), round(pair[1], k)]


def read_strata(site):
    """Map stratum of every sample point, by point id."""
    stratum = {}
    with open(site.products / "validation_sample.csv") as fh:
        for r in csv.DictReader(fh):
            stratum[r["id"]] = r["map_stratum"]
    return stratum


def labelled_late(r):
    """True when the land-cover label (pass 1) was given after the labeller saw the time series (pass 2)."""
    return bool(r.get("pass2_at") and r.get("pass1_at") and r["pass1_at"] > r["pass2_at"])


def s2_viewed(r):
    """True when the Sentinel-2 layer was open while the land-cover label was given."""
    return str(r.get("s2_viewed", "0")) in ("1", "true", "True")


def assess(h, ref, mapped, paddy, W, A):
    """Overall, user's and producer's accuracy and error-adjusted areas (km2), each with its 95% half-interval.

    None when a stratum has fewer than two judged points.
    """
    if any((h == s).sum() < 2 for s in W):
        return None
    ones = np.ones(len(h))
    out = {"overall_accuracy": rounded(ratio((mapped == ref).astype(float), ones, h, W)),
           "users_accuracy": {}, "producers_accuracy": {}, "error_adjusted_area_km2": {}}
    for c in REF_CLASSES:
        agree = ((mapped == c) & (ref == c)).astype(float)
        out["users_accuracy"][c] = rounded(ratio(agree, (mapped == c).astype(float), h, W))
    for c in REF_CLASSES:
        agree = ((mapped == c) & (ref == c)).astype(float)
        out["producers_accuracy"][c] = rounded(ratio(agree, (ref == c).astype(float), h, W))
    for c in REF_CLASSES + ["paddy"]:
        y = paddy if c == "paddy" else (ref == c).astype(float)
        p, ci = ratio(y, ones, h, W)
        out["error_adjusted_area_km2"][c] = [round(A * p, 1), round(A * ci, 1)] if p is not None else None
    return out


def time_series_pass(pts):
    """Pass 2 per stratum: share flooded in the wet season and share with two or more crop cycles, with counts."""
    flood, cyc = defaultdict(list), defaultdict(list)
    for s, _, r in pts:
        if r.get("wet_season_flooding") in ("yes", "no"):
            flood[s].append(r["wet_season_flooding"] == "yes")
        if r.get("crop_cycles") in ("0", "1", "2", "3"):
            cyc[s].append(int(r["crop_cycles"]))
    flooded = {}
    for s, v in flood.items():
        flooded[s] = [round(100 * float(np.mean(v)), 1), len(v)]
    cropped = {}
    for s, v in cyc.items():
        cropped[s] = [round(100 * float(np.mean(np.array(v) >= 2)), 1), len(v)]
    return {"note": "Interpreted from the same Sentinel-1 and Sentinel-2 series the maps use; consistency, not accuracy.",
            "flooded_in_wet_season_pct": flooded,
            "two_or_more_crop_cycles_pct": cropped}


def evaluate(key, rows):
    """Accuracy, error-adjusted areas and the pass 2 consistency check of one site, from its label rows."""
    site = Site(key)
    design = site.load_json("validation_design")
    stratum = read_strata(site)
    area = {s: design["stratum_area_km2"].get(s, 0.0) for s in STRATA}
    A = sum(area.values())                         # study area, km2
    W = {s: a / A for s, a in area.items() if a > 0}  # area share of each stratum

    pts = [(stratum[r["id"]], int(r["cover_code"]), r) for r in rows if r["id"] in stratum]
    late = [r["id"] for _, _, r in pts if labelled_late(r)]
    judged = [(s, c, r) for s, c, r in pts if c in COVER and r["id"] not in late]
    h = np.array([s for s, _, _ in judged])                  # stratum of each judged point
    ref = np.array([COVER[c] for _, c, _ in judged])         # its class from the imagery
    paddy = np.array([c == PADDY for _, c, _ in judged], float)
    mapped = np.array([MAP_CLASS[s] for s, _, _ in judged])  # its class on our maps
    s2 = np.array([s2_viewed(r) for _, _, r in judged])
    res = {"site": site.meta["name"], "points_labelled": len(pts), "points_judged": len(judged),
           "points_rejected_pass1_after_pass2": late, "points_with_sentinel2_viewed": int(s2.sum()),
           "points_by_stratum": {s: int((h == s).sum()) for s in STRATA}, "study_area_km2": round(A, 1),
           "basis": "Pass 1, land cover from imagery only; independent of the maps."}
    short = [s for s in W if (h == s).sum() < 2]
    if short:
        res["note"] = f"No estimates yet: label at least two judged points in every stratum (missing: {', '.join(short)})."
        site.save_json("sample_accuracy", res)
        return res

    res.update(assess(h, ref, mapped, paddy, W, A))
    if s2.any():
        keep = ~s2
        res["without_sentinel2_viewed"] = assess(h[keep], ref[keep], mapped[keep], paddy[keep], W, A)
    else:
        res["without_sentinel2_viewed"] = "same as above"

    mapped_area = {}
    for c in REF_CLASSES:
        mapped_area[c] = round(sum([a for s, a in area.items() if MAP_CLASS[s] == c]), 1)
    paddy_share, cover_share = {}, {}
    for s in STRATA:
        in_s = h == s
        if not in_s.sum():
            continue
        paddy_share[s] = round(100 * float(paddy[in_s].mean()), 1)
        cover_share[s] = {c: round(100 * float((ref[in_s] == c).mean()), 1) for c in REF_CLASSES}
    res["mapped_area_km2"] = mapped_area
    res["paddy_share_by_stratum_pct"] = paddy_share
    res["cover_by_stratum_pct"] = cover_share
    res["time_series_pass"] = time_series_pass(pts)
    site.save_json("sample_accuracy", res)
    return res


def read_labels(path):
    """Label rows of the labels CSV, grouped by site."""
    by_site = defaultdict(list)
    with open(path) as fh:
        for r in csv.DictReader(fh):
            by_site[r["site"]].append(r)
    return by_site


if __name__ == "__main__":
    labs = read_labels(sys.argv[1])
    for key in sys.argv[2:] or [k for k in config.SITES if labs.get(k)]:
        r = evaluate(key, labs[key])
        print(key, r["points_judged"], "points judged, overall accuracy", r.get("overall_accuracy"))
