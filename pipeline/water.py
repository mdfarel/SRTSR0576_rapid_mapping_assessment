"""Water bodies and hydraulic features.

Dry-season open water is mapped at 10 m from the Sentinel-2 median composite
(NDWI above a scene-adaptive Otsu threshold, confirmed by MNDWI > 0 to reject
built-up surfaces, slopes above 15 degrees masked). Each connected water body is
then described by its shape: length along its principal axis, mean width,
elongation, how completely it fills its oriented bounding rectangle, whether it
encloses dry land, and how closely it aligns with the cardinal directions.
Rectilinearity is the four-fold coherence of edge directions along each body's
outline: 1 when every edge follows one pair of perpendicular directions, near 0
for a meandering or rounded outline. It holds for a reservoir that is only part
full, because its straight embankments still set the edge directions, and the
same pair of directions gives the body's deviation from the cardinal axes.
Those measures sort features into rivers, straight canals, enclosing moats,
rectilinear reservoirs and ponds, and irregular lakes and ponds. Shape alone
cannot date a feature, so the classes are morphological candidates, not
archaeological identifications.

Seasonal behaviour comes from monthly Sentinel-1 VV backscatter (open water is
radar-dark), which sees through wet-season cloud.

Inputs: the Sentinel-2 dry-season composites (10 m and 20 m), the DEM, the
Sentinel-1 monthly stack, JRC surface-water occurrence and the study-area polygon.
Outputs (data/<site>/products): water_dry_10m, water_class_10m,
water_cardinal_10m, s1_water_months_20m, s1_water_dry_20m, s1_water_wet_20m and
pond_density_20m (GeoTIFFs), water_stats.json, and water_features.csv with one
row per water body.

    python -m pipeline.water [site ...]
"""
import csv
import sys

import numpy as np
from pyproj import Transformer
from scipy import ndimage

import config
from pipeline.layers import Site, edge_otsu, grid_convergence_deg, ndi, slope_deg, upsample

CLASSES = {
    1: "River or stream",
    2: "Straight canal (candidate)",
    3: "Enclosing moat (candidate)",
    4: "Rectilinear reservoir or moat segment (10 ha or more)",
    5: "Rectilinear pond",
    6: "Lake or irregular reservoir",
    7: "Irregular pond",
    8: "Small water body (shape not resolved)",
}
RECTILINEAR = (2, 3, 4, 5)  # straight-sided classes, whose orientation is tested
POND_CLASSES = (5, 7, 8)    # counted in the pond density

MIN_PX = 10           # 0.1 ha minimum mapping unit
SHAPE_MIN_PX = 20     # below 0.2 ha the outline is too blocky to judge shape
CARDINAL_TOL = 7.5    # degrees either side of north-south and east-west
RECTILINEAR_R = 0.40  # edge-direction coherence; about the top 15% of bodies at Angkor
MAX_SLOPE_S2 = 15     # degrees; steeper ground is left out of the optical water map
MAX_SLOPE_S1 = 10     # degrees; steeper ground is left out of the radar water map

# form classes (see classify); areas in 10 m pixels, lengths and widths in metres
LARGE_PX = 1000                  # 10 ha: reservoirs and lakes rather than ponds
MOAT_MIN_ENCLOSED = 1.5          # the closed ring holds at least half its own area again as land
MOAT_MAX_WIDTH_M = 300
CHANNEL_MIN_ELONGATION = 5       # length at least five times the breadth
CHANNEL_MIN_LENGTH_M = 500
CANAL_MAX_WIDTH_M = 60
RIVER_MIN_LENGTH_M = 1000        # long, irregular bodies that fill little of their bounding rectangle
RIVER_MAX_RECTANGULARITY = 0.35

# ring test for enclosing moats (see enclosure_ratio)
RING_GAP_ALL = 3         # closing iterations on all water: bridges gaps of up to about 60 m
RING_GAP_SEGMENTS = 10   # closing iterations on straight-sided bodies: bridges gaps of up to about 200 m
SEGMENT_MIN_PX = 100     # 1 ha: smallest straight-sided body in the wide-gap ring test
RING_MIN_PX = 200        # smallest closed ring scored (2 ha) ...
RING_MIN_EXTENT_M = 300  # ... and it must span at least 300 m
BORDER_PX = 3            # a body must lie within 3 pixels of the land its ring encloses


def dry_indices(site):
    """Dry-season NDWI (10 m), MNDWI (20 m, repeated onto the 10 m grid) and slope in degrees (10 m)."""
    s10 = site.s2_composite("dry", 10)
    s20 = site.s2_composite("dry", 20)
    ndwi = ndi(s10["B03"], s10["B08"])
    mndwi = upsample(ndi(s20["B03"], s20["B11"]), 2)
    slope = slope_deg(site.static10("dem"), 10)
    return ndwi, mndwi, slope


def water_mask(ndwi, mndwi, slope, inside, t, mndwi_t=0.0):
    """Open water: NDWI above t and MNDWI above mndwi_t on gentle slopes, as bodies of MIN_PX pixels or more."""
    # no morphological opening: it erases rivers and canals one or two pixels wide;
    # bodies under the 0.1 ha minimum mapping unit are dropped instead
    water = (ndwi > t) & (mndwi > mndwi_t) & (slope < MAX_SLOPE_S2) & inside
    lab, _ = ndimage.label(water, np.ones((3, 3)))
    sizes = np.bincount(lab.ravel())
    keep = sizes >= MIN_PX
    keep[0] = False  # label 0 is the background
    return keep[lab]


def dry_water_sensitivity(ndwi, mndwi, slope, inside, t_ndwi, px_km2):
    """Dry-season water area (km2) with the NDWI threshold moved by 0.05 and the MNDWI threshold by 0.1."""
    areas = []
    for dt in (-0.05, 0.0, 0.05):
        for mndwi_t in (-0.1, 0.0, 0.1):
            mask = water_mask(ndwi, mndwi, slope, inside, t_ndwi + dt, mndwi_t)
            areas.append(float(mask.sum() * px_km2))
    return areas


def principal_axis(r, c, lid, n):
    """Area (pixels), centroid, major-axis angle and extent along and across that axis of each body.

    r and c are the row and column of every water pixel and lid its body label (1 to n; index 0 is
    unused). The axis comes from the second moments of the pixel positions, and the two extents (in
    pixels) are measured along and across it, so their product is the oriented bounding rectangle.
    """
    x, y = c.astype(np.float64), -r.astype(np.float64)  # y points north
    area = np.bincount(lid, minlength=n + 1).astype(np.float64)
    a = np.maximum(area, 1)
    mx, my = np.bincount(lid, x, n + 1) / a, np.bincount(lid, y, n + 1) / a
    dx, dy = x - mx[lid], y - my[lid]
    cxx = np.bincount(lid, dx * dx, n + 1) / a
    cyy = np.bincount(lid, dy * dy, n + 1) / a
    cxy = np.bincount(lid, dx * dy, n + 1) / a
    theta = 0.5 * np.arctan2(2 * cxy, cxx - cyy)  # major axis, radians from east
    # position of every pixel along (u) and across (v) the major axis of its body
    ux, uy = np.cos(theta), np.sin(theta)
    pu = dx * ux[lid] + dy * uy[lid]
    pv = -dx * uy[lid] + dy * ux[lid]
    umax, vmax = np.full(n + 1, -1e9), np.full(n + 1, -1e9)
    umin, vmin = np.full(n + 1, 1e9), np.full(n + 1, 1e9)
    np.maximum.at(umax, lid, pu)
    np.minimum.at(umin, lid, pu)
    np.maximum.at(vmax, lid, pv)
    np.minimum.at(vmin, lid, pv)
    length_px = umax - umin + 1
    breadth_px = vmax - vmin + 1
    return area, mx, my, theta, length_px, breadth_px


def outline_rectilinearity(water, lab, n):
    """Four-fold coherence of the outline edge directions of each body, and the angle of its main axis pair.

    Each outline pixel votes with its edge direction at four times the angle, so directions 90 degrees
    apart reinforce each other. The length of the summed vote over the total weight is 1 when every edge
    follows one pair of perpendicular directions and near 0 for a rounded or meandering outline. The
    angle of the axis pair is returned in degrees counter-clockwise from image east.
    """
    smooth = ndimage.gaussian_filter(water.astype(np.float32), 1.0)
    gx, gy = ndimage.sobel(smooth, axis=1), -ndimage.sobel(smooth, axis=0)  # rows run south, so flip y
    er, ec = np.nonzero(water & ~ndimage.binary_erosion(water))  # outline pixels
    eid = lab[er, ec]
    phi = np.arctan2(gy[er, ec], gx[er, ec])
    wgt = np.hypot(gx[er, ec], gy[er, ec])
    s4 = np.bincount(eid, wgt * np.cos(4 * phi), n + 1) + 1j * np.bincount(eid, wgt * np.sin(4 * phi), n + 1)
    rectilinear = np.abs(s4) / np.maximum(np.bincount(eid, wgt, n + 1), 1e-9)
    return rectilinear, np.degrees(np.angle(s4) / 4)


def enclosure_ratio(mask, lab, n, res, gap_iterations):
    """Ring test: for each body, the area of its closed ring plus the land inside, over the ring's own area.

    `mask` is closed over gaps of about 2 x gap_iterations pixels (causeways, bridges). Each closed
    component of at least RING_MIN_PX pixels that spans RING_MIN_EXTENT_M is scored by its filled area
    over its own area, so 1.5 means it encloses half its own area again as land. A water body (labelled
    in `lab`) takes the score of the component it falls in, but only if it lies within BORDER_PX of the
    enclosed land; otherwise it scores 0.
    """
    r, c = np.nonzero(lab)
    lid = lab[r, c]
    closed = ndimage.binary_closing(mask, ndimage.generate_binary_structure(2, 1), iterations=gap_iterations)
    land = ndimage.binary_fill_holes(closed) & ~closed
    borders = np.bincount(lid, ndimage.binary_dilation(land, iterations=BORDER_PX)[r, c], n + 1) > 0
    clab, cn = ndimage.label(closed, np.ones((3, 3)))
    ratio = np.zeros(cn + 1)
    carea = np.bincount(clab.ravel(), minlength=cn + 1)
    for j, sl in enumerate(ndimage.find_objects(clab), start=1):
        if sl is None or carea[j] < RING_MIN_PX:
            continue
        h, w = sl[0].stop - sl[0].start, sl[1].stop - sl[1].start
        if max(h, w) * res >= RING_MIN_EXTENT_M:
            ratio[j] = ndimage.binary_fill_holes(clab[sl] == j).sum() / carea[j]
    owner = np.zeros(n + 1, int)  # the closed component each body falls in
    owner[lid] = clab[r, c]
    return np.where(borders, ratio[owner], 0)


def describe(water, res, convergence):
    """Shape measures of every connected water body, as arrays indexed by body label (0 is the background)."""
    lab, n = ndimage.label(water, np.ones((3, 3)))  # 8-connected bodies
    r, c = np.nonzero(lab)
    area, mx, my, theta, length_px, breadth_px = principal_axis(r, c, lab[r, c], n)
    rect = area / np.maximum(length_px * breadth_px, 1)  # share of the oriented bounding rectangle filled
    azimuth = (90 - np.degrees(theta) + convergence) % 180  # major axis, degrees from true north
    elong = length_px / np.maximum(breadth_px, 1)

    # four-fold coherence of outline directions, and the deviation of that axis pair from cardinal;
    # the convergence is the true bearing of grid north, so it is added to bearings measured clockwise
    # from grid north (as for the azimuth above) and subtracted from angles counter-clockwise from east
    rectilinear, axis_from_east = outline_rectilinearity(water, lab, n)
    axis = axis_from_east - convergence
    dev = np.abs(((axis + 45) % 90) - 45)

    # enclosure: causeways and bridges cut most moats, so the ring test runs twice, on all water
    # closed over gaps of up to about 60 m, and on straight-sided bodies of 1 ha or more closed over
    # gaps of up to about 200 m (so that a moat cut in two by its causeways still closes into a ring);
    # a body inside a closed ring that holds at least half its own area again as land is enclosing,
    # provided it borders the enclosed land itself (a pond just outside the ring is swept into the
    # closed component but does not touch the land inside)
    segments = (rectilinear >= RECTILINEAR_R) & (area >= SEGMENT_MIN_PX)
    segments[0] = False
    ring_all = enclosure_ratio(water, lab, n, res, RING_GAP_ALL)
    ring_segments = enclosure_ratio(segments[lab], lab, n, res, RING_GAP_SEGMENTS)
    enclosed = np.maximum(ring_all, np.where(segments, ring_segments, 0))
    return dict(lab=lab, n=n, area=area, length=length_px * res, width=area * res / np.maximum(length_px, 1),
                elong=elong, rect=rect, rectilinear=rectilinear, enclosed=enclosed, azimuth=azimuth,
                deviation=dev, col=mx, row=-my)


def classify(f):
    """Form class (a key of CLASSES) and cardinal flag of each body described by `f`.

    The first rule that applies wins: too small to judge (8); a straight-sided ring around land (3);
    long and narrow, a canal if straight and narrow (2), otherwise a river (1); long, irregular and
    filling little of its bounding rectangle (1); other straight-sided bodies by size (4 or 5); the rest
    by size (6 or 7). A straight-sided body is cardinal when its axis pair lies within CARDINAL_TOL of
    north-south and east-west.
    """
    n = f["n"]
    cls = np.zeros(n + 1, np.uint8)
    cardinal = np.zeros(n + 1, bool)
    for i in range(1, n + 1):
        area_px, length, width = f["area"][i], f["length"][i], f["width"][i]
        straight = f["rectilinear"][i] >= RECTILINEAR_R
        if area_px < SHAPE_MIN_PX:
            cls[i] = 8
        elif f["enclosed"][i] >= MOAT_MIN_ENCLOSED and straight and width <= MOAT_MAX_WIDTH_M:
            cls[i] = 3  # straight-sided ring around land; rivers around islands fail the straightness test
        elif f["elong"][i] >= CHANNEL_MIN_ELONGATION and length >= CHANNEL_MIN_LENGTH_M:
            cls[i] = 2 if (straight and width <= CANAL_MAX_WIDTH_M) else 1
        elif length >= RIVER_MIN_LENGTH_M and not straight and f["rect"][i] < RIVER_MAX_RECTANGULARITY:
            cls[i] = 1
        elif straight:
            cls[i] = 4 if area_px >= LARGE_PX else 5
        elif area_px >= LARGE_PX:
            cls[i] = 6
        else:
            cls[i] = 7
        cardinal[i] = cls[i] in RECTILINEAR and f["deviation"][i] <= CARDINAL_TOL
    return cls, cardinal


def rectilinear_deviations(mask, res):
    """Deviation from the image axes (degrees, 0 to 45) of every straight-sided body in a water mask."""
    f = describe(mask, res, 0.0)
    cls, _ = classify(f)
    sel = np.isin(cls, RECTILINEAR)
    sel[0] = False
    return f["deviation"][sel]


def share_within_tol(dev, near):
    """Percentage of deviations within CARDINAL_TOL of `near` degrees, or None when there are none."""
    if dev.size == 0:
        return None
    return round(100 * float((np.abs(dev - near) <= CARDINAL_TOL).mean()), 1)


def rotation_control(water, res=10, angle=30.0):
    """Pixel-grid control for water-body orientation, as for field edges: rotate a central window
    by 30 degrees and measure the share of rectilinear bodies near the image axes (grid bias) and
    near the rotated ground axes (real pattern). Raster outlines staircase along the image axes,
    so a share that stays on the image axes after rotation would be an artefact."""
    h, w = water.shape
    half = min(h, w) // 3
    window = water[h // 2 - half: h // 2 + half, w // 2 - half: w // 2 + half]
    rotated = ndimage.rotate(window.astype(np.float32), angle, reshape=True, order=1) > 0.5
    dev_window = rectilinear_deviations(window, res)
    dev_rotated = rectilinear_deviations(rotated, res)
    n_window, n_rotated = int(dev_window.size), int(dev_rotated.size)
    image_axes = share_within_tol(dev_rotated, 0.0)
    ground_axes = share_within_tol(dev_rotated, angle)
    if n_window >= 30 and n_rotated >= 30:
        passed = bool(ground_axes > image_axes + 10)
    else:
        passed = None  # too few bodies to judge
    return {"window_rectilinear_bodies": n_window, "window_share_pct": share_within_tol(dev_window, 0.0),
            "after_30deg_rotation_image_axes_pct": image_axes, "after_30deg_rotation_original_axes_pct": ground_axes,
            "passed": passed,
            "rule": "passes when, after rotation, the share on the ground axes exceeds the share on the image axes by "
                    "more than 10 points; not judged with fewer than 30 rectilinear bodies"}


def s1_water(site):
    """Radar water from monthly Sentinel-1 VV (open water is radar-dark), on the 20 m grid.

    Returns the VV threshold and its bimodality, the number of months with water (255 where never
    observed or outside the study area), dry- and wet-season water, and the wet-season and
    inundation areas in pixels with the threshold moved by 1 dB either way.
    """
    vv, _ = site.s1_monthly()
    slope20 = slope_deg(site.static10("dem")[::2, ::2], 20)
    dry_idx = [m - 1 for m in site.season_months("dry")]  # months 1-12 to array index 0-11
    wet_idx = [m - 1 for m in site.season_months("wet")]
    with np.errstate(all="ignore"):
        vv_dry = np.nanmedian(vv[dry_idx], axis=0)
    inside = site.mask(20)
    gentle = slope20 < MAX_SLOPE_S1
    # Otsu on pixels near water edges, first split at -18 dB, held within [-24, -13] dB,
    # and the conventional -18 dB if the edge sample is not bimodal
    t, eta = edge_otsu(vv_dry, inside & gentle, -18.0, -24.0, -13.0, water_high=False, fallback=-18.0)
    ok = gentle & inside
    monthly = (vv < t) & ok
    valid = np.isfinite(vv)
    months = np.where(valid.any(0) & inside, (monthly & valid).sum(0), 255).astype(np.uint8)
    water_dry = (vv_dry < t) & ok
    water_wet = monthly[wet_idx].any(0)
    # sensitivity: the same maps with the threshold moved by 1 dB either way
    sens_wet, sens_inundation = [], []
    for dt in (-1.0, 1.0):
        wet_dt = ((vv < t + dt) & ok)[wet_idx].any(0)
        dry_dt = (vv_dry < t + dt) & ok
        sens_wet.append(wet_dt.sum())
        sens_inundation.append((wet_dt & ~dry_dt).sum())
    return dict(threshold=t, eta=eta, months=months, dry=water_dry, wet=water_wet,
                sens_wet=sens_wet, sens_inundation=sens_inundation)


def pond_density(f, ponds, grid):
    """Pond-class bodies per km2 (by centroid) in a 1 km moving window, sampled onto the 20 m grid."""
    pts = np.zeros(f["lab"].shape, np.float32)
    rr = np.clip(f["row"][1:].round().astype(int), 0, grid.height - 1)
    cc = np.clip(f["col"][1:].round().astype(int), 0, grid.width - 1)
    np.add.at(pts, (rr[ponds[1:]], cc[ponds[1:]]), 1)
    density = ndimage.uniform_filter(pts, size=100) * 100 * 100  # mean over 100 x 100 pixels (1 km2), as a count
    return density[::2, ::2].astype(np.float32)


def km2_range(pixel_counts, px_km2):
    """[smallest, largest] of a list of pixel counts, in km2 to two decimals."""
    return [round(float(min(pixel_counts) * px_km2), 2), round(float(max(pixel_counts) * px_km2), 2)]


def ndwi_method(t, eta):
    """How the NDWI threshold was set, in words."""
    if eta < 0.3:
        return "conventional 0.0 (edge sample not bimodal)"
    if t <= -0.0999:
        return "edge Otsu, held at the -0.10 lower bound"
    return "edge Otsu"


def s1_method(t, eta):
    """How the Sentinel-1 VV threshold was set, in words."""
    if eta < 0.3:
        return "conventional -18 dB (edge sample not bimodal)"
    if t in (-24.0, -13.0):
        return "edge Otsu, held at a bound"
    return "edge Otsu"


def write_features_csv(site, f, cls, cardinal):
    """water_features.csv: one row per water body with its class, shape measures and centroid in lon/lat."""
    g10 = site.g10
    to_ll = Transformer.from_crs(g10.crs, "EPSG:4326", always_xy=True)
    xs = g10.transform.c + (f["col"][1:] + 0.5) * 10
    ys = g10.transform.f - (f["row"][1:] + 0.5) * 10
    lon, lat = to_ll.transform(xs, ys)
    with open(site.products / "water_features.csv", "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["id", "class", "area_ha", "length_m", "mean_width_m", "elongation", "rectangularity",
                    "rectilinearity", "deviation_from_cardinal_deg", "encloses_land_ratio", "azimuth_deg",
                    "cardinal", "lon", "lat"])
        for i in range(1, f["n"] + 1):  # area in pixels of 0.01 ha
            w.writerow([i, CLASSES[cls[i]], round(f["area"][i] * 0.01, 2), round(f["length"][i]),
                        round(f["width"][i], 1), round(f["elong"][i], 2), round(f["rect"][i], 2),
                        round(f["rectilinear"][i], 2), round(f["deviation"][i], 1), round(f["enclosed"][i], 2),
                        round(f["azimuth"][i], 1), int(cardinal[i]),
                        round(lon[i - 1], 5), round(lat[i - 1], 5)])


def run(key):
    """Water maps, water_stats.json and water_features.csv for one site."""
    site = Site(key)
    g10, g20 = site.g10, site.g20
    inside10 = site.mask(10)

    # 1. dry-season open water at 10 m (Sentinel-2). The NDWI threshold is Otsu on pixels near water
    # edges, first split at 0.0, held within [-0.10, 0.25], and the conventional 0.0 if not bimodal
    ndwi, mndwi, slope = dry_indices(site)
    t_ndwi, eta_ndwi = edge_otsu(ndwi, inside10 & (slope < MAX_SLOPE_S2), 0.0, -0.10, 0.25, fallback=0.0)
    water = water_mask(ndwi, mndwi, slope, inside10, t_ndwi)
    valid10 = np.isfinite(ndwi) & inside10
    sens_dry = dry_water_sensitivity(ndwi, mndwi, slope, inside10, t_ndwi, g10.pixel_area_km2)
    del ndwi, mndwi, slope  # free memory before the shape analysis

    # 2. shape and form class of each water body
    f = describe(water, 10, grid_convergence_deg(key))
    cls, cardinal = classify(f)
    lab = f["lab"]
    site.save("water_dry_10m", water.astype(np.uint8), g10, nodata=None, dtype="uint8")
    site.save("water_class_10m", cls[lab], g10, nodata=0, dtype="uint8")
    site.save("water_cardinal_10m", (cardinal[lab] & (lab > 0)).astype(np.uint8), g10, nodata=0, dtype="uint8")

    # 3. seasonal water at 20 m (Sentinel-1)
    s1 = s1_water(site)
    site.save("s1_water_months_20m", s1["months"], g20, nodata=255, dtype="uint8")
    site.save("s1_water_dry_20m", s1["dry"].astype(np.uint8), g20, dtype="uint8")
    site.save("s1_water_wet_20m", s1["wet"].astype(np.uint8), g20, dtype="uint8")

    # 4. pond density
    ponds = np.isin(cls, POND_CLASSES)
    ponds[0] = False
    site.save("pond_density_20m", pond_density(f, ponds, g10), g20, dtype="float32")

    # 5. agreement with JRC surface water present in at least half of all observations (above 100 is no data)
    jrc = site.static10("jrc_occurrence")
    jrc_perm = (jrc >= 50) & (jrc <= 100) & inside10
    inter, union = (water & jrc_perm).sum(), (water | jrc_perm).sum()

    # 6. statistics
    px10, px20 = g10.pixel_area_km2, g20.pixel_area_km2
    area_km2 = float(inside10.sum() * px10)
    counts = {CLASSES[k]: int((cls[1:] == k).sum()) for k in CLASSES}
    areas = {CLASSES[k]: round(float(f["area"][1:][cls[1:] == k].sum() * px10), 3) for k in CLASSES}
    card = {CLASSES[k]: int((cardinal[1:] & (cls[1:] == k)).sum()) for k in RECTILINEAR}
    rect_total = int(np.isin(cls[1:], RECTILINEAR).sum())
    irregular = cls[1:] == 7
    if irregular.any():
        irregular_cardinal_pct = round(100 * float((f["deviation"][1:][irregular] <= CARDINAL_TOL).mean()), 1)
    else:
        irregular_cardinal_pct = None
    months = s1["months"]
    stats = {
        "site": site.meta["name"],
        "study_area_km2": round(area_km2, 1),
        "dry_season_water_km2": round(float(water.sum() * px10), 2),
        "dry_season_water_share_pct": round(float(water.sum() / inside10.sum() * 100), 2),
        "s1_dry_water_km2": round(float(s1["dry"].sum() * px20), 2),
        "s1_wet_max_water_km2": round(float(s1["wet"].sum() * px20), 2),
        "s1_permanent_water_km2": round(float(((months >= 10) & (months <= 12)).sum() * px20), 2),
        "s1_intermittent_water_km2": round(float(((months >= 1) & (months <= 9)).sum() * px20), 2),
        "seasonal_inundation_km2": round(float((s1["wet"] & ~s1["dry"]).sum() * px20), 2),
        "feature_count": int(f["n"]),
        "features_per_100km2": round(f["n"] / area_km2 * 100, 1),
        "ponds_per_100km2": round(int(ponds.sum()) / area_km2 * 100, 1),
        "count_by_class": counts,
        "area_km2_by_class": areas,
        "cardinal_rectilinear_by_class": card,
        "cardinal_share_of_rectilinear_pct": round(100 * sum(card.values()) / max(rect_total, 1), 1),
        "cardinal_share_of_irregular_ponds_pct": irregular_cardinal_pct,
        "orientation_control": rotation_control(water),
        "jrc_agreement": {
            "jrc_occurrence_50pct_km2": round(float(jrc_perm.sum() * px10), 2),
            "iou": round(float(inter / max(union, 1)), 3),
        },
        "threshold_method": {
            "ndwi": ndwi_method(t_ndwi, eta_ndwi),
            "s1_vv": s1_method(s1["threshold"], s1["eta"])},
        # s1_vv_db_otsu repeats s1_vv_db_edge_otsu under the older name that pipeline.fields reads
        "thresholds": {"ndwi_edge_otsu": round(t_ndwi, 3), "ndwi_bimodality": eta_ndwi,
                       "s1_vv_db_edge_otsu": round(s1["threshold"], 2), "s1_vv_bimodality": s1["eta"],
                       "s1_vv_db_otsu": round(s1["threshold"], 2)},
        "sensitivity": {
            "dry_season_water_km2_range": [round(min(sens_dry), 2), round(max(sens_dry), 2)],
            "dry_season_water_note": "NDWI threshold +/-0.05 and MNDWI threshold -0.1, 0, +0.1",
            "s1_wet_max_water_km2_range": km2_range(s1["sens_wet"], px20),
            "seasonal_inundation_km2_range": km2_range(s1["sens_inundation"], px20),
            "s1_note": "Sentinel-1 VV threshold +/-1 dB",
        },
        "s2_dry_clear_obs_median": int(np.median(site.s2_count("dry", 10)[inside10])),
        "s2_dry_valid_share_pct": round(float(valid10.sum() / inside10.sum() * 100), 1),
    }
    site.save_json("water_stats", stats)
    write_features_csv(site, f, cls, cardinal)

    print(key, {k: stats[k] for k in ("dry_season_water_km2", "s1_wet_max_water_km2",
                                       "feature_count", "cardinal_share_of_rectilinear_pct")})
    print("  ", counts)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
