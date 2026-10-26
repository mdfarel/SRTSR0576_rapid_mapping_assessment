"""Stratified random validation sample per site (design only; labels to be added by an interpreter).

This is the documented next step of the rapid assessment, not part of its results: the sample is
drawn from the Waterfields maps alone and no point has been labelled. Strata are the mapped
classes, so that user's and producer's accuracy and error-adjusted areas with confidence intervals
can be computed once each point has been labelled (Olofsson et al. 2014). Labelling is planned in
two passes: land cover first, then the seasonal behaviour from the monthly series, so that the
independent label is never shaped by the data the maps were made from. Points are drawn with a
fixed seed, so the sample can be regenerated exactly, and written in random order, so any number
labelled from the start of the list is still a stratified random sample.

Inputs: water_dry_10m.tif, water_change_20m.tif and field_class_20m.tif.
Outputs: products/validation_sample.csv and .geojson, and products/validation_design.json.

    python -m pipeline.sample [site ...]
"""
import csv
import json
import sys

import numpy as np
from pyproj import Transformer

import config
from pipeline.layers import Site, block_mean

TOTAL = 300        # points per site
MIN_PER_STRATUM = 25  # every mapped class gets enough points for a user's accuracy; the rest follow area
SEED = 20261023
CSV_HEADER = ["id", "map_stratum", "lon", "lat", "reference_label", "interpreter", "imagery_used", "notes"]


def strata(site):
    """The strata of one site as 20 m masks, in order of priority."""
    inside = site.mask(20)
    # a 20 m cell is open water when at least half of its four 10 m pixels are
    water = (block_mean(site.load("water_dry_10m").astype(np.float32), 2) >= 0.5) & inside
    wc = site.load("water_change_20m")
    fields = site.load("field_class_20m")
    # strata must not overlap: a cell belongs to the first stratum that claims it, so seasonal
    # inundation here is wet-season water outside open water and outside the mapped fields
    order = [
        ("open water, dry season", water),
        ("flood-recession fields", fields == 1),
        ("irrigated or multi-cropped", fields == 2),
        ("rainfed rice landscape", fields == 3),
        ("open agricultural land", fields == 4),
        ("seasonal inundation", wc == 2),
        ("other land", inside),
    ]
    taken = ~inside
    s = {}
    for name, m in order:
        s[name] = m & ~taken
        taken |= s[name]
    return s


def site_seed(key):
    """The fixed seed, shifted by the site key so each site draws a different sample."""
    return SEED + sum(ord(ch) for ch in key)


def allocate(sizes):
    """Points per stratum: MIN_PER_STRATUM each (or every cell of a smaller stratum), the rest in proportion to area."""
    total_px = max(sum(sizes.values()), 1)
    alloc = {n: min(sz, MIN_PER_STRATUM) for n, sz in sizes.items()}
    spare = TOTAL - sum(alloc.values())
    for n, sz in sizes.items():
        alloc[n] = min(sz, alloc[n] + int(round(max(spare, 0) * sz / total_px)))
    return alloc


def draw_points(st, alloc, grid, rng):
    """Random 20 m cell centres from each stratum, as (stratum, lon, lat), stratum by stratum."""
    to_ll = Transformer.from_crs(grid.crs, "EPSG:4326", always_xy=True)
    points = []
    for name, m in st.items():
        idx = np.flatnonzero(m)
        if idx.size == 0:
            continue
        pick = rng.choice(idx, size=alloc[name], replace=False)
        r, c = np.unravel_index(pick, m.shape)
        x = grid.transform.c + (c + 0.5) * grid.res
        y = grid.transform.f - (r + 0.5) * grid.res
        lon, lat = to_ll.transform(x, y)
        points += [(name, round(a, 6), round(b, 6)) for a, b in zip(lon, lat)]
    return points


def write_sample(site, key, points):
    """The sample as a CSV for labelling and as GeoJSON for the map, numbered in the given order."""
    with open(site.products / "validation_sample.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(CSV_HEADER)
        for i, (name, lon, lat) in enumerate(points, start=1):
            w.writerow([f"{key}-{i:04d}", name, lon, lat, "", "", "", ""])
    features = []
    for i, (name, lon, lat) in enumerate(points, start=1):
        features.append({"type": "Feature", "properties": {"id": f"{key}-{i:04d}", "map_stratum": name},
                         "geometry": {"type": "Point", "coordinates": [lon, lat]}})
    gj = {"type": "FeatureCollection", "features": features}
    (site.products / "validation_sample.geojson").write_text(json.dumps(gj))


def run(key):
    """Draw the validation sample of one site and write it with its design."""
    site = Site(key)
    rng = np.random.default_rng(site_seed(key))
    st = strata(site)
    sizes = {n: int(m.sum()) for n, m in st.items()}
    alloc = allocate(sizes)
    area_km2 = {n: round(float(m.sum() * site.g20.pixel_area_km2), 2) for n, m in st.items()}
    points = draw_points(st, alloc, site.g20, rng)
    # number the points only after shuffling, so an id says nothing about the stratum it came from
    points = [points[i] for i in rng.permutation(len(points))]
    write_sample(site, key, points)
    site.save_json("validation_design", {"total": TOTAL, "min_per_stratum": MIN_PER_STRATUM, "allocation": alloc,
                                         "seed": SEED, "stratum_area_km2": area_km2,
                                         "points": len(points), "reference": "Olofsson et al. 2014, RSE 148:42-57"})
    print(key, len(points), "points in", len(area_km2), "strata", flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
