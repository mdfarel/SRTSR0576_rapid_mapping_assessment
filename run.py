"""SRTSR0576 rapid mapping assessment: run the whole process, from the downloads to the web GIS, for one or more sites.

    python3 run.py                                  all four sites, every step, then the comparative dataset
    python3 run.py angkor                           one site, then the comparative dataset
    python3 run.py --from water                     rerun the analysis only, downloads kept
    python3 run.py angkor --from "crop phenology"   one site, from that step on
    python3 run.py --from "comparative dataset"     only the cross-site CSV and charts
    python3 run.py --from maps --gis                maps onwards, then the QGIS projects

The steps (listed at the end of this help) run in order for each site, and the comparative
dataset follows once every site is done. The first four steps download the satellite data and
skip whatever is already in data/, so an interrupted run resumes; the analysis steps always
rerun. The Fields of The World step also downloads its polygons once and keeps them in data/.
The web GIS in webgis/ is a static page with no build of its own: the "web tiles" step writes
each site's layers and numbers into webgis/tiles/ and webgis/data/, which the page loads, and
the comparative dataset step copies its CSV there too. Step names with a space need quotes
after --from.

--gis builds the QGIS project of each site in out/gis/<site>/ with the Python that has PyQGIS,
after the comparative dataset.
"""
import argparse
import subprocess

import config
from pipeline import (accuracy, change, dataset, fetch_s1, fetch_s2, fetch_s2_monthly, fetch_static, fields, ftw,
                      maps, phenology, sample, sources, tiles, water)

# (name, function called with the site key), in the order they run for each site.
STEPS = [
    ("static layers", fetch_static.fetch),
    ("Sentinel-2 composites", fetch_s2.run),
    ("Sentinel-1 monthly", fetch_s1.run),
    ("Sentinel-2 monthly NDVI", fetch_s2_monthly.run),
    ("water", water.run),
    ("crop phenology", phenology.run),
    ("fields", fields.run),
    ("change", change.run),
    ("water sources", sources.run),
    ("Fields of The World", ftw.analyse),
    ("agreement with reference maps", accuracy.run),
    ("validation sample design", sample.run),
    ("maps", maps.run),
    ("web tiles", tiles.run),
]
DATASET_STEP = "comparative dataset"  # runs once, after every site: the cross-site CSV and charts
STEP_NAMES = [name for name, _ in STEPS] + [DATASET_STEP]
QGIS_PYTHON = "python3"  # the system Python, which has PyQGIS (pipeline/gis_project.py)


def parse_args():
    """Site keys to process, the index of the first step to run, and the --gis flag."""
    parser = argparse.ArgumentParser(
        prog="python3 run.py", description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="sites: " + ", ".join(config.SITES) + "\nsteps: " + ", ".join(STEP_NAMES))
    parser.add_argument("sites", nargs="*", metavar="site", help="site keys to process (default: all)")
    parser.add_argument("--from", dest="first_step", metavar="STEP", default=STEP_NAMES[0],
                        help="start at this step instead of the first")
    parser.add_argument("--gis", action="store_true", help="also build the QGIS project of each site")
    args = parser.parse_intermixed_args()  # sites may come before or after --from
    for key in args.sites:
        if key not in config.SITES:
            parser.error(f"unknown site {key!r}; choose from {', '.join(config.SITES)}")
    if args.first_step not in STEP_NAMES:
        parser.error(f"unknown step {args.first_step!r}; choose from {', '.join(STEP_NAMES)}")
    return args.sites or list(config.SITES), STEP_NAMES.index(args.first_step), args.gis


def main():
    """Every site step from the chosen one on, then the comparative dataset and the optional extras."""
    sites, start, gis = parse_args()
    for key in sites:
        for name, step in STEPS[start:]:
            print(f"[{key}] {name}", flush=True)
            step(key)
    print(f"[all sites] {DATASET_STEP}", flush=True)
    dataset.run()
    tiles.copy_comparative()  # the web GIS offers the CSV for download
    if gis:
        print("[QGIS projects]", ", ".join(sites), flush=True)
        subprocess.run([QGIS_PYTHON, "-m", "pipeline.gis_project", *sites], cwd=config.ROOT, check=True)


if __name__ == "__main__":
    main()
