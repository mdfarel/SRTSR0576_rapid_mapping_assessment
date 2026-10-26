"""Agricultural field systems, their relation to water, and field-boundary orientation.

Classes (20 m), assigned in this order of precedence inside the agricultural domain:

1. Flood-recession fields: inundated for two or more wet-season months (Sentinel-1)
   and green again in the dry season (Sentinel-2 NDVI 0.35 or more).
2. Irrigated or multi-cropped: green in the dry season (NDVI 0.45 or more, with a
   seasonal NDVI amplitude of 0.25 or more and at least one crop cycle, so that
   evergreen orchards and plantations do not count) or cropped two or more times in
   the year (pipeline.phenology), outside the flood zone. In a single-monsoon
   climate either needs water from storage, canals, groundwater or residual moisture.
3. Rice-field landscape (rainfed): the paddy radar signature, a flooded low in VH
   backscatter followed by a rise of several dB as the crop grows.
4. Open agricultural land: the rest of the agricultural domain (upland crops,
   fallow, fields without a detectable paddy cycle).

The agricultural domain is ESA WorldCover 2021 cropland, extended to grassland and
herbaceous wetland pixels that carry the paddy signature, since WorldCover often
files paddies under those classes, and to land outside built-up areas and
permanent water that shows a crop calendar (two or more cycles, or one cycle over
shrubland, with an NDVI amplitude of 0.30 or more), which brings in the small
sawah of Central Java.
Slopes of 20 degrees or more are excluded.

Field-boundary orientation uses the structure tensor of the 10 m dry-season NDVI:
at each pixel it gives the dominant direction of edges and how coherent that
direction is. Over cropland this measures the grain of the field pattern.

Inputs: ESA WorldCover, the DEM, the Sentinel-1 monthly stack, the Sentinel-2
dry-season composites (10 m and 20 m) and monthly NDVI, and the products of
pipeline.water (water_dry_10m, s1_water_*, water_stats.json) and
pipeline.phenology (crop_cycles_20m, ndvi_amplitude_20m).
Outputs (data/<site>/products): field_class_20m, paddy_signature_20m,
distance_to_perennial_water_20m, ag_intensity_20m, field_grain_deviation_100m and
field_grain_coherence_100m (GeoTIFFs) and field_stats.json.

    python -m pipeline.fields [site ...]
"""
import sys

import numpy as np
from rasterio.transform import from_origin
from scipy import ndimage

import config
from pipeline.common import Grid, read_tif
from pipeline.layers import Site, block_mean, grid_convergence_deg, ndi, slope_deg, upsample
from pipeline.phenology import MIN_CLEAR_MONTHS, NOT_JUDGED, PEAK_MIN, PROMINENCE, fill_gaps

CLASSES = {
    1: "Flood-recession fields",
    2: "Irrigated or multi-cropped",
    3: "Rice-field landscape (rainfed)",
    4: "Open agricultural land",
}
RICE_VH_MIN = -19.0      # dB, flooded or bare paddy before the crop closes
RICE_VH_RANGE = 5.0      # dB, rise from flooded low to crop peak
NDVI_IRRIGATED = 0.45    # dry-season NDVI of a growing crop
AMP_CROP = 0.30          # seasonal NDVI amplitude of a crop cycle
AMP_IRRIGATED = 0.25     # minimum amplitude for dry-season greenness to count as cropping
NDVI_RECESSION = 0.35    # dry-season NDVI of a crop sown as the flood recedes
MAX_SLOPE = 20           # degrees; steeper ground is left out of the agricultural domain
CARDINAL_TOL = 7.5       # degrees either side of north-south and east-west
MIN_COHERENCE = 0.5      # structure-tensor coherence of an edge used in the orientation statistics

# ESA WorldCover 2021 class codes
WC_SHRUB = 20
WC_GRASS = 30
WC_CROP = 40
WC_BUILT = 50
WC_WATER = 80     # permanent water bodies
WC_WETLAND = 90   # herbaceous wetland


def classify(site):
    """Field class (a key of CLASSES, 0 outside agriculture) at 20 m, and the layers and tests reported with it."""
    # land cover: share of each 20 m cell (four 10 m pixels) in the WorldCover classes
    wc = site.static10("worldcover")
    crop = block_mean((wc == WC_CROP).astype(np.float32), 2)
    grassy = block_mean(np.isin(wc, (WC_GRASS, WC_WETLAND)).astype(np.float32), 2)
    shrubby = block_mean(np.isin(wc, (WC_SHRUB, WC_GRASS, WC_WETLAND)).astype(np.float32), 2)
    built = block_mean(np.isin(wc, (WC_BUILT, WC_WATER)).astype(np.float32), 2)  # built-up or permanent water

    # radar: the paddy signature (a low in VH, then a rise as the crop grows) and wet-season flooding
    vv, vh = site.s1_monthly()
    with np.errstate(all="ignore"):
        vh_min, vh_max = np.nanmin(vh, 0), np.nanmax(vh, 0)
    months = site.load("s1_water_months_20m")  # months with radar water in the year (255 = no data)
    wet_idx = [m - 1 for m in site.season_months("wet")]  # months 1-12 to array index 0-11
    t_vv = site.load_json("water_stats")["thresholds"]["s1_vv_db_otsu"]
    wet_flood_months = (vv[wet_idx] < t_vv).sum(0)
    # water in six months or fewer, so that permanent water is not taken for paddy
    rice = (vh_min <= RICE_VH_MIN) & (vh_max - vh_min >= RICE_VH_RANGE) & (months <= 6)

    # optical: dry-season greenness, and the crop calendar from pipeline.phenology
    dry = site.s2_composite("dry", 20)
    ndvi_dry = ndi(dry["B08"], dry["B04"])
    dry_ok = site.s2_count("dry", 20) >= 1  # at least one clear dry-season view
    slope20 = slope_deg(site.static10("dem")[::2, ::2], 20)
    cyc = site.load("crop_cycles_20m")
    amp = np.nan_to_num(site.load("ndvi_amplitude_20m"))
    judged = cyc != NOT_JUDGED
    multi = (cyc >= 2) & judged

    # phenology extension of the domain: double-cropped land with a strong NDVI cycle is farmland whatever
    # WorldCover calls it (Java's sawah is often labelled tree cover); single-cycle land is added only from
    # grass, shrub or herbaceous wetland, so that dry-season leaf fall in deciduous forest is not counted
    cropping = (multi & (amp >= AMP_CROP)) | ((cyc == 1) & (amp >= AMP_CROP) & (shrubby >= 0.5))
    domain = (crop >= 0.5) | ((grassy >= 0.5) & rice) | (cropping & (built < 0.5))
    ag = domain & (slope20 < MAX_SLOPE) & site.mask(20)

    # classes, in order of precedence. Green in the dry season counts as irrigated only with a crop cycle
    # behind it: perennial tree crops and home gardens are green all year with a flat NDVI curve
    green = (dry_ok & (ndvi_dry >= NDVI_IRRIGATED) & (amp >= AMP_IRRIGATED) & (cyc >= 1) & judged) | multi
    recession = ag & (wet_flood_months >= 2) & dry_ok & (ndvi_dry >= NDVI_RECESSION) & (months <= 9)
    irrigated = ag & green & ~recession
    rainfed = ag & rice & ~recession & ~irrigated
    open_ag = ag & ~recession & ~irrigated & ~rainfed
    cls = np.zeros(ag.shape, np.uint8)
    cls[recession] = 1
    cls[irrigated] = 2
    cls[rainfed] = 3
    cls[open_ag] = 4

    # sensitivity of the paddy signature and of dry-season greenness to their thresholds
    base_ag = ((crop >= 0.5) | (grassy >= 0.5)) & (slope20 < MAX_SLOPE) & site.mask(20)
    px = site.g20.pixel_area_km2
    paddy_km2 = []
    for vmin in (-18.0, -19.0, -20.0, -21.0):
        for rise in (5.0, 6.5, 8.0):
            paddy = base_ag & (vh_min <= vmin) & (vh_max - vh_min >= rise) & (months <= 6)
            paddy_km2.append(float(paddy.sum() * px))
    green_km2 = {}
    for t in (0.35, 0.45, 0.55):
        green_km2[str(t)] = round(float((ag & dry_ok & (ndvi_dry >= t)).sum() * px), 1)
    sensitivity = {"paddy_signature_km2_range": [round(min(paddy_km2), 1), round(max(paddy_km2), 1)],
                   "paddy_note": "VH low -18 to -21 dB, VH rise 5, 6.5 or 8 dB",
                   "dry_season_green_km2_by_ndvi": green_km2}
    return cls, dict(paddy=ag & rice, cycles=cyc, dry_green=ag & dry_ok & (ndvi_dry >= NDVI_IRRIGATED),
                     multi=ag & multi, sensitivity=sensitivity)


def water_relation(site, cls):
    """Distance from each field pixel to perennial (dry-season) water."""
    w10 = site.load("water_dry_10m").astype(np.float32)
    perennial = (block_mean(w10, 2) > 0) | (site.load("s1_water_dry_20m") == 1)
    dist = ndimage.distance_transform_edt(~perennial) * 20.0  # pixels to metres
    land = ~perennial & site.mask(20)
    base_1km = float((dist[land] <= 1000).mean())
    out = {"all_land": {"median_distance_m": float(np.median(dist[land])),
                        "within_500m_pct": round(100 * float((dist[land] <= 500).mean()), 1),
                        "within_1km_pct": round(100 * base_1km, 1)}}
    for k, name in CLASSES.items():
        d = dist[cls == k]
        if d.size == 0:
            continue
        share = float((d <= 1000).mean())
        out[name] = {"median_distance_m": float(np.median(d)),
                     "within_500m_pct": round(100 * float((d <= 500).mean()), 1),
                     "within_1km_pct": round(100 * share, 1),
                     "enrichment_vs_all_land_1km": round(share / max(base_1km, 1e-6), 2)}
    return dist.astype(np.float32), out


def tensor_orientation(jxx, jyy, jxy):
    """Coherence (0 to 1) and boundary azimuth of a structure tensor.

    The azimuth is in degrees clockwise from grid north and is not yet wrapped into 0 to 180.
    """
    with np.errstate(all="ignore"):
        coherence = np.sqrt((jxx - jyy) ** 2 + 4 * jxy ** 2) / (jxx + jyy)
    grad_angle = 0.5 * np.degrees(np.arctan2(2 * jxy, jxx - jyy))  # dominant gradient, from east, counter-clockwise
    # a boundary runs at right angles to the gradient (+90), and azimuths run clockwise from north (90 - angle)
    return coherence, 90 - (grad_angle + 90)


def tensor(ndvi):
    """Structure tensor of an image: its three components, edge energy, coherence and boundary azimuth (0 to 180)."""
    gx = ndimage.sobel(ndvi, axis=1)
    gy = -ndimage.sobel(ndvi, axis=0)  # rows run south, so flip to make y point north
    jxx = ndimage.gaussian_filter(gx * gx, 2)
    jyy = ndimage.gaussian_filter(gy * gy, 2)
    jxy = ndimage.gaussian_filter(gx * gy, 2)
    energy = jxx + jyy
    coherence, azimuth = tensor_orientation(jxx, jyy, jxy)
    return jxx, jyy, jxy, energy, coherence, azimuth % 180  # boundary line, from grid north


def deviation_from_axes(az):
    """Angle in degrees (0 to 45) between an azimuth and the nearest of north-south and east-west."""
    return np.minimum(az % 90, 90 - az % 90)


def strong_edges(energy, coherence, region):
    """Coherent edges (MIN_COHERENCE or more) in the stronger half of edge energy within a region."""
    return region & (coherence >= MIN_COHERENCE) & (energy > np.percentile(energy[region], 50))


def cardinal_share(edge_azimuth, energy, coherence, region, offset=0.0):
    """Energy-weighted share of strong, coherent edges within CARDINAL_TOL of the axes rotated by `offset`.

    Returns the share (0 to 1), the 5-degree histogram of edge azimuths and the number of pixels used.
    """
    if not region.any():
        return float("nan"), None, 0
    sel = strong_edges(energy, coherence, region)
    w = energy[sel]
    az = (edge_azimuth[sel] - offset) % 180
    dev = deviation_from_axes(az)
    hist, _ = np.histogram(az, bins=36, range=(0, 180), weights=w)
    share = float((w * (dev <= CARDINAL_TOL)).sum() / max(w.sum(), 1e-9))
    return share, hist / max(hist.sum(), 1e-9), int(sel.sum())


def block_bootstrap(edge_azimuth, energy, coherence, region, block=200, n=1000, seed=20261023):
    """95% interval of the cardinal share from resampling 2 km blocks (spatially correlated edges)."""
    if not region.any():
        return None
    sel = strong_edges(energy, coherence, region)
    r, c = np.nonzero(sel)
    az = edge_azimuth[r, c] % 180
    w = energy[r, c]
    hit = w * (deviation_from_axes(az) <= CARDINAL_TOL)
    # number the blocks, and total the edge weight and the near-cardinal weight in each
    bid = (r // block) * (edge_azimuth.shape[1] // block + 1) + c // block
    _, inv = np.unique(bid, return_inverse=True)
    w_block, hit_block = np.bincount(inv, w), np.bincount(inv, hit)
    if w_block.size < 10:
        return None
    # n draws, each of as many blocks as there are, picked with replacement
    pick = np.random.default_rng(seed).integers(0, w_block.size, (n, w_block.size))
    shares = hit_block[pick].sum(1) / w_block[pick].sum(1)
    return [round(100 * float(np.percentile(shares, 2.5)), 1), round(100 * float(np.percentile(shares, 97.5)), 1)]


def rotation_test(ndvi, ag10, edge_azimuth, energy, coherence):
    """Control 2: cardinal share in a central window, and again after rotating the image by 30 degrees."""
    h, w = ndvi.shape
    c0, r0, half = w // 2, h // 2, min(h, w) // 3
    win = (slice(r0 - half, r0 + half), slice(c0 - half, c0 + half))
    rot = ndimage.rotate(ndvi[win], 30, reshape=False, order=1, mode="nearest")
    rot_ag = ndimage.rotate(ag10[win].astype(np.float32), 30, reshape=False, order=0) > 0.5
    # keep the centre only: rotation pads the corners of the window
    core = np.zeros_like(rot_ag)
    margin = int(half * 0.3)
    core[margin:-margin, margin:-margin] = True
    region = rot_ag & core
    _, _, _, e2, c2, a2 = tensor(rot)
    rot_image_axes, _, _ = cardinal_share(a2, e2, c2, region)
    # the sign of the rotation depends on the array orientation, so take the better of the two
    rot_minus = cardinal_share(a2, e2, c2, region, offset=-30.0)[0]
    rot_plus = cardinal_share(a2, e2, c2, region, offset=30.0)[0]
    unrot_window, _, _ = cardinal_share(edge_azimuth[win], energy[win], coherence[win], ag10[win])
    return {"window_share_pct": round(100 * unrot_window, 1),
            "after_30deg_rotation_image_axes_pct": round(100 * rot_image_axes, 1),
            "after_30deg_rotation_original_axes_pct": round(100 * max(rot_minus, rot_plus), 1)}


def local_grain(jxx, jyy, jxy, ag10, convergence):
    """Field grain on a 100 m grid: the tensor of agricultural pixels summed over a 1 km window.

    Returns the deviation of the dominant boundary direction from cardinal (degrees) and its coherence;
    the deviation is blank (NaN) where under 30% of the window is agricultural.
    """
    k = 100  # 1 km at 10 m
    sxx = ndimage.uniform_filter(np.where(ag10, jxx, 0), k)[::10, ::10]
    syy = ndimage.uniform_filter(np.where(ag10, jyy, 0), k)[::10, ::10]
    sxy = ndimage.uniform_filter(np.where(ag10, jxy, 0), k)[::10, ::10]
    agf = ndimage.uniform_filter(ag10.astype(np.float32), k)[::10, ::10]
    coh, azimuth = tensor_orientation(sxx, syy, sxy)
    az = (azimuth + convergence) % 180
    local_dev = deviation_from_axes(az)
    local_dev[(agf < 0.3) | ~np.isfinite(coh)] = np.nan
    return local_dev.astype(np.float32), coh.astype(np.float32)


def orientation(site, ag20, convergence):
    """Structure-tensor orientation of field boundaries on 10 m dry-season NDVI, with two controls.

    Control 1 repeats the statistic on coherent edges outside the agricultural
    domain. Control 2 rotates a central window of the image by 30 degrees and
    measures the share around the image axes: a method free of pixel-grid bias
    finds the same low share as for random directions, while the real pattern
    moves to the rotated axes.
    """
    s10 = site.s2_composite("dry", 10)
    ndvi = ndi(s10["B08"], s10["B04"])
    ndvi = np.where(np.isfinite(ndvi), ndvi, np.nanmedian(ndvi)).astype(np.float32)  # gaps get the median
    jxx, jyy, jxy, energy, coherence, edge_azimuth = tensor(ndvi)
    edge_azimuth = (edge_azimuth + convergence) % 180  # from true north

    ag10 = upsample(ag20, 2)[: ndvi.shape[0], : ndvi.shape[1]]
    share, hist, used = cardinal_share(edge_azimuth, energy, coherence, ag10)
    expected = 4 * CARDINAL_TOL / 180  # the share if edge directions were random
    other = site.mask(10) & ~ag10
    control_share, _, _ = cardinal_share(edge_azimuth, energy, coherence, other)
    controls = {"non_agricultural_edges_pct": round(100 * control_share, 1),
                "non_agricultural_ci95_pct": block_bootstrap(edge_azimuth, energy, coherence, other),
                "rotation_test": rotation_test(ndvi, ag10, edge_azimuth, energy, coherence)}

    local_dev, coh = local_grain(jxx, jyy, jxy, ag10, convergence)
    stats = {"histogram_5deg": [round(float(x), 4) for x in hist],
             "cardinal_share_pct": round(100 * share, 1),
             "cardinal_share_ci95_pct": block_bootstrap(edge_azimuth, energy, coherence, ag10),
             "uniform_expectation_pct": round(100 * expected, 1),
             "cardinal_ratio": round(share / expected, 2),
             "pixels_used": used,
             "controls": controls}
    return local_dev, coh, stats


def class_ndvi_curves(site, cls):
    """Mean gap-filled NDVI of each field class for each month (twelve None for a class that is absent)."""
    filled = fill_gaps(read_tif(site.s2 / f"ndvi_monthly_{config.YEAR}.tif").astype(np.float32))
    curves = {}
    for k, name in CLASSES.items():
        in_class = cls == k
        if in_class.any():
            curves[name] = [round(float(np.nanmean(filled[m][in_class])), 3) for m in range(12)]
        else:
            curves[name] = [None] * 12
    return curves


def run(key):
    """Field classes, distance to water, field grain and field_stats.json for one site."""
    site = Site(key)
    cls, aux = classify(site)
    site.save("field_class_20m", cls, site.g20, nodata=0, dtype="uint8")
    site.save("paddy_signature_20m", aux["paddy"].astype("uint8"), site.g20, dtype="uint8")
    dist, relation = water_relation(site, cls)
    site.save("distance_to_perennial_water_20m", dist, site.g20, dtype="float32")

    ag = cls > 0
    win = 50  # 1 km at 20 m
    ag_density = ndimage.uniform_filter(ag.astype(np.float32), win)
    irr_density = ndimage.uniform_filter((cls == 2).astype(np.float32), win)
    site.save("ag_intensity_20m", ag_density, site.g20, dtype="float32")

    local_dev, coh, orient = orientation(site, ag, grid_convergence_deg(key))
    g = site.g10
    g100 = Grid(g.epsg, from_origin(g.transform.c, g.transform.f, 100, 100),
                local_dev.shape[1], local_dev.shape[0], 100)
    site.save("field_grain_deviation_100m", local_dev, g100, nodata=np.nan, dtype="float32")
    site.save("field_grain_coherence_100m", coh, g100, nodata=np.nan, dtype="float32")

    px = site.g20.pixel_area_km2
    inside = site.mask(20)
    total = float(inside.sum() * px)
    n_ag = max(int(ag.sum()), 1)

    # cropping intensity: share of judged agricultural land with 0, 1, 2 or 3 crop cycles
    cyc = aux["cycles"]
    judged = ag & (cyc != NOT_JUDGED)
    n_judged = max(int(judged.sum()), 1)
    intensity = {str(c): round(100 * float(((cyc == c) & judged).sum()) / n_judged, 1) for c in range(4)}
    stats = {
        "site": site.meta["name"],
        "agricultural_km2": round(float(ag.sum() * px), 1),
        "agricultural_share_pct": round(float(ag.sum() / inside.sum() * 100), 1),
        "cropping_cycles_share_pct": intensity,
        "multi_cropped_km2": round(float(aux["multi"].sum() * px), 1),
        "dry_season_green_km2": round(float(aux["dry_green"].sum() * px), 1),
        "ndvi_monthly_mean_by_class": class_ndvi_curves(site, cls),
        "area_km2_by_class": {n: round(float((cls == k).sum() * px), 1) for k, n in CLASSES.items()},
        "share_of_agriculture_pct": {n: round(100 * float((cls == k).sum()) / n_ag, 1) for k, n in CLASSES.items()},
        "concentrated_ag_km2": round(float((ag_density >= 0.8).sum() * px), 1),
        "concentrated_irrigated_km2": round(float((irr_density >= 0.5).sum() * px), 1),
        "paddy_signature_km2": round(float(aux["paddy"].sum() * px), 1),
        "worldcover_cropland_km2": round(float(((site.static10("worldcover") == WC_CROP) & site.mask(10)).sum()
                                               * site.g10.pixel_area_km2), 1),
        "relation_to_perennial_water": relation,
        "field_orientation": orient,
        "sensitivity": aux["sensitivity"],
        "thresholds": {"rice_vh_min_db": RICE_VH_MIN, "rice_vh_range_db": RICE_VH_RANGE,
                       "ndvi_irrigated": NDVI_IRRIGATED, "ndvi_recession": NDVI_RECESSION,
                       "amplitude_crop": AMP_CROP, "amplitude_irrigated": AMP_IRRIGATED},
        "study_area_km2": round(total, 1),
        "crop_cycle_rule": {"peak_ndvi_min": PEAK_MIN, "prominence_min": PROMINENCE,
                            "min_clear_months": MIN_CLEAR_MONTHS},
    }
    site.save_json("field_stats", stats)
    print(key, stats["agricultural_km2"], stats["area_km2_by_class"], "cardinal ratio", orient["cardinal_ratio"])


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        run(k)
