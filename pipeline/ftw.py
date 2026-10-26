"""Fields of The World 2024 field boundaries inside each study area.

The global FTW map (Robinson et al. 2026, CC-BY-4.0) is a polygonised 10 m
prediction of the PRUE model. Polygons are read straight from the per-country
GeoParquet on Source Cooperative, filtered by bounding box so that only the
row groups covering the study area are fetched.

Polygon outlines follow the 10 m pixel staircase, which is aligned with the
UTM grid and therefore close to north-south and east-west. Orientation is
therefore measured from each field's minimum rotated rectangle, and only for
fields of 2 ha or more that are at least 1.5 times longer than wide, where
the rectangle is set by the field's shape rather than by its pixels.

Inputs: the FTW GeoParquet (downloaded once to products/ftw_fields_2024.gpkg)
and our field classes (field_class_20m.tif).
Outputs: products/ftw_fields_20m.tif (FTW fields inside the polygon) and
products/ftw_stats.json.

    python -m pipeline.ftw [site ...]
"""
import sys

import duckdb
import numpy as np
import shapely
from pyproj import Transformer
from rasterio.features import rasterize

import config
from pipeline.common import aoi_utm, bbox_wgs84, site_grid
from pipeline.layers import Site, grid_convergence_deg

BASE = "https://data.source.coop/ftw/global-data/predictions/vectors/alpha/results-by-admin-conf"
COUNTRY = {"angkor": "KH/Cambodia", "srei_santhor": "KH/Cambodia", "sukhothai": "TH/Thailand",
           "prambanan": "ID/Indonesia"}
YEAR = 2024
CARDINAL_TOL = 7.5    # degrees either side of north-south or east-west that count as cardinal
MIN_FIELD_HA = 2      # smaller fields are too few pixels across for their rectangle to show their shape
MIN_ELONGATION = 1.5  # long side over short side; a near-square rectangle has no clear long axis
HIST_BINS = 36        # 5-degree bins over 0 to 180 degrees


def fetch(key):
    """Download the FTW 2024 polygons that overlap the site grid to a GeoPackage, once."""
    site = Site(key)
    out = site.products / "ftw_fields_2024.gpkg"
    if out.exists():
        return out
    cc, name = COUNTRY[key].split("/")
    url = f"{BASE}/admin:country_code={cc}/{name}.parquet"
    x0, y0, x1, y1 = bbox_wgs84(site_grid(key, 10))
    con = duckdb.connect()
    con.sql("INSTALL httpfs; LOAD httpfs; INSTALL spatial; LOAD spatial; SET geometry_always_xy = true;")
    tmp = out.with_suffix(".partial.gpkg")
    tmp.unlink(missing_ok=True)
    con.sql(f"""
        COPY (
            SELECT id, "metrics:area" AS area_m2, confidence, "determination:datetime" AS observed, geometry
            FROM read_parquet('{url}')
            WHERE bbox.xmax >= {x0} AND bbox.xmin <= {x1} AND bbox.ymax >= {y0} AND bbox.ymin <= {y1}
              AND year("determination:datetime") = {YEAR}
        ) TO '{tmp}' WITH (FORMAT GDAL, DRIVER 'GPKG', LAYER_NAME 'ftw_fields_2024')""")
    tmp.replace(out)  # only a complete download gets the final name
    return out


def lonlat_to_utm(geoms, crs):
    """Reproject shapely geometries from lon/lat to the site's UTM coordinates."""
    to_utm = Transformer.from_crs("EPSG:4326", crs, always_xy=True)

    def project(xy):
        return np.column_stack(to_utm.transform(xy[:, 0], xy[:, 1]))

    return shapely.transform(geoms, project)


def read_fields(key, crs, aoi):
    """FTW polygons (UTM) whose interior point lies in the study-area polygon, and their confidence."""
    path = fetch(key)
    con = duckdb.connect()
    con.sql("INSTALL spatial; LOAD spatial; SET geometry_always_xy = true;")
    rows = con.sql(f"SELECT area_m2, confidence, ST_AsWKB(geom) FROM ST_Read('{path}')").fetchall()
    geoms = shapely.from_wkb([r[2] for r in rows])
    conf = np.array([np.nan if r[1] is None else r[1] for r in rows], float)
    utm = lonlat_to_utm(geoms, crs)
    inside = shapely.intersects(shapely.point_on_surface(utm), aoi)
    return utm[inside], conf[inside]


def long_axis(utm):
    """Long-side vector and elongation (long over short side) of each field's minimum rotated rectangle."""
    long_side, e1, e2 = [], [], []
    for rect in shapely.minimum_rotated_rectangle(utm):
        c = np.asarray(rect.exterior.coords)
        side1, side2 = c[1] - c[0], c[2] - c[1]  # two adjacent sides, as (dx, dy)
        len1, len2 = np.hypot(*side1), np.hypot(*side2)
        long_side.append(side1 if len1 >= len2 else side2)
        e1.append(len1)
        e2.append(len2)
    e1, e2 = np.array(e1), np.array(e2)
    elong = np.maximum(e1, e2) / np.maximum(np.minimum(e1, e2), 1e-6)
    return np.array(long_side), elong


def analyse(key):
    """Field counts, sizes, orientation and overlap with our agricultural land for one site."""
    site = Site(key)
    aoi = shapely.Polygon(aoi_utm(key))
    utm, conf = read_fields(key, site.g10.crs, aoi)
    area_ha = shapely.area(utm) / 1e4

    # orientation of the long axis, in degrees clockwise from north (0 to 180), corrected for grid convergence
    long_vec, elong = long_axis(utm)
    azimuth = (np.degrees(np.arctan2(long_vec[:, 0], long_vec[:, 1])) + grid_convergence_deg(key)) % 180
    use = (area_ha >= MIN_FIELD_HA) & (elong >= MIN_ELONGATION)
    dev = np.minimum(azimuth % 90, 90 - azimuth % 90)  # degrees away from the nearest cardinal direction
    hist, _ = np.histogram(azimuth[use], bins=HIST_BINS, range=(0, 180))
    hist = hist / max(hist.sum(), 1)

    ag = site.load("field_class_20m") > 0
    ftw_mask = (rasterize(((g, 1) for g in utm), out_shape=ag.shape, transform=site.g20.transform, fill=0,
                          dtype="uint8") == 1) & site.mask(20)
    site.save("ftw_fields_20m", ftw_mask.astype(np.uint8), site.g20, nodata=0, dtype="uint8")
    px = site.g20.pixel_area_km2
    stats = {
        "source": "Fields of The World 2024 (Robinson et al. 2026), PRUE predictions, CC-BY-4.0",
        "fields": int(len(utm)),
        "fields_per_km2": round(len(utm) / (aoi.area / 1e6), 2),
        "field_area_km2": round(float(ftw_mask.sum() * px), 1),
        "median_field_ha": round(float(np.median(area_ha)), 2) if len(area_ha) else None,
        "field_ha_quartiles": [round(float(q), 2) for q in np.percentile(area_ha, [25, 50, 75])] if len(area_ha) else None,
        "confidence_median": None if np.isnan(conf).all() else round(float(np.nanmedian(conf)), 1),
        "confidence_null_share_pct": round(100 * float(np.isnan(conf).mean()), 1) if len(conf) else None,
        "share_of_our_agriculture_covered_pct": round(100 * float((ftw_mask & ag).sum()) / max(int(ag.sum()), 1), 1),
        "share_of_ftw_inside_our_agriculture_pct": round(100 * float((ftw_mask & ag).sum()) / max(int(ftw_mask.sum()), 1), 1),
        "orientation": {"fields_used": int(use.sum()),
                        "cardinal_share_pct": round(100 * float((dev[use] <= CARDINAL_TOL).mean()), 1) if use.any() else None,
                        "uniform_expectation_pct": round(100 * 2 * CARDINAL_TOL / 90, 1),
                        "histogram_5deg": [round(float(h), 4) for h in hist]},
    }
    site.save_json("ftw_stats", stats)
    print(key, {k: stats[k] for k in ("fields", "median_field_ha", "field_area_km2", "share_of_our_agriculture_covered_pct")},
          "cardinal", stats["orientation"]["cardinal_share_pct"], flush=True)


if __name__ == "__main__":
    for k in sys.argv[1:] or list(config.SITES):
        analyse(k)
