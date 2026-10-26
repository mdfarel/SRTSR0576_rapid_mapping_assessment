"""Study areas, seasons and grids for the rapid landscape-signature scoping task.

Srei Santhor, Sukhothai and Prambanan use the study-area polygons supplied by
A/Prof Scott Hawken with the task ("SRT AOI.kml", 22 October 2026). The file
has no Angkor polygon, so Angkor uses a 44.72 km square (about 2,000 km²)
centred on 13.400 N, 103.880 E. Each site is analysed on the UTM grid
that encloses its polygon, and every statistic is computed inside the polygon.

Seasons follow the regional monsoon: the mainland sites are dry in February to
April and wettest in late August to October; Central Java is the reverse.

Every pipeline step reads its settings from here; nothing is run by this file.
"""
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
OUT = ROOT / "out"

SIDE_M = 44_720  # Angkor square: 44.72 km x 44.72 km = 1,999.9 km²
YEAR = 2024      # most recent complete year with strong wet/dry contrast at all four sites

# Season windows (start, end), inclusive. The dates must fall within YEAR.
MAINLAND_SEASONS = {
    "dry": ("2024-02-01", "2024-04-30"),
    "wet": ("2024-08-15", "2024-10-31"),
}
JAVA_SEASONS = {
    "wet": ("2024-01-01", "2024-03-31"),
    "dry": ("2024-07-01", "2024-09-30"),
}

# One entry per site, keyed by the name used on the command line and in data/<site>/.
#   lat, lon    site centre: the Angkor square is centred here, and grid north is compared with true north here
#   epsg        UTM zone of the analysis grid
#   seasons     MAINLAND_SEASONS or JAVA_SEASONS
#   aoi         study-area polygon as (lon, lat) corners, or None for the SIDE_M square; aoi_source says where from
SITES = {
    "angkor": {
        "name": "Angkor",
        "country": "Cambodia",
        # centred so that the large rectilinear reservoirs and the lake margin fall inside the square
        "lat": 13.400, "lon": 103.880,
        "epsg": 32648,
        "seasons": MAINLAND_SEASONS,
        "aoi": None,  # square of SIDE_M around lat/lon
        "aoi_source": "44.72 km square centred on 13.400 N, 103.880 E (no Angkor polygon in SRT AOI.kml)",
    },
    "srei_santhor": {
        "name": "Srei Santhor and region",
        "country": "Cambodia",
        "lat": 11.778, "lon": 105.026,
        "epsg": 32648,
        "seasons": MAINLAND_SEASONS,
        "aoi": [(104.6830894435153, 11.9353694052719), (104.6813657268019, 11.62069779894017),
                (105.370056541653, 11.61621412145414), (105.3665092403897, 11.93993162683726)],
        "aoi_source": "SRT AOI.kml: Srei Santhor and Region",
    },
    "sukhothai": {
        "name": "Sukhothai region",
        "country": "Thailand",
        "lat": 17.215, "lon": 99.745,
        "epsg": 32647,
        "seasons": MAINLAND_SEASONS,
        "aoi": [(99.58013037293878, 17.49853846049007), (99.58038049316266, 16.92147366375391),
                (99.9157425296379, 16.93014175318412), (99.89931091849226, 17.5104670511704)],
        "aoi_source": "SRT AOI.kml: Sukhothai Region",
    },
    "prambanan": {
        "name": "Prambanan region",
        "country": "Indonesia",
        "lat": -7.668, "lon": 110.356,
        "epsg": 32749,
        "seasons": JAVA_SEASONS,
        "aoi": [(110.14053611424, -7.53753460798359), (110.1362311631844, -7.788271573056234),
                (110.5756584443073, -7.802011974810339), (110.5735036097574, -7.534079681937389)],
        "aoi_source": "SRT AOI.kml: supplied polygon",
    },
}

# Sentinel-2 bands in the seasonal composites: blue, green, red, near infrared and the two shortwave
# infrared bands at 20 m; green, red and near infrared at 10 m.
S2_BANDS_20 = ["B02", "B03", "B04", "B08", "B11", "B12"]
S2_BANDS_10 = ["B03", "B04", "B08"]
# Scene classification classes kept as clear: dark area (2), vegetation (4), bare (5), water (6),
# unclassified (7) and snow (11). Dark, still water such as a narrow moat is often labelled 2 or 7,
# so dropping those classes cuts moats and canals into pieces. Saturated or defective (1), cloud
# shadow (3), cloud (8, 9) and cirrus (10) pixels are dropped; the per-pixel median removes the rest.
S2_VALID_SCL = (2, 4, 5, 6, 7, 11)
