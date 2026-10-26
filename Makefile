# Waterfields: one target per job. Each target runs the command shown under it, so `make run`
# and `python3 run.py` do the same thing. Use another Python with, for example,
# `make run PYTHON=.venv/bin/python`.

PYTHON ?= python3
SITE ?= angkor

.PHONY: help install run site all gis dataset serve deploy

# List the targets (what `make` alone does)
help:
	@echo "make install | run | site SITE=angkor | all | gis | dataset | serve | deploy"

# Python packages from requirements.txt (the GDAL tools and QGIS are installed separately)
install:
	$(PYTHON) -m pip install -r requirements.txt

# All four sites, then the cross-site dataset; downloads are cached in data/ and resume after a break
run:
	$(PYTHON) run.py

# One site, for example `make site SITE=sukhothai`
site:
	$(PYTHON) run.py $(SITE)

# Everything: all four sites, the cross-site dataset and the QGIS projects
all:
	$(PYTHON) run.py --gis

# QGIS project of each site in out/gis/ (needs the system Python that has PyQGIS)
gis:
	python3 -m pipeline.gis_project

# Cross-site CSV and comparison charts, from the sites already processed
dataset:
	$(PYTHON) run.py --from "comparative dataset"

# Web GIS (webgis/) at http://localhost:8765
serve:
	$(PYTHON) serve.py 8765

# Publish webgis/ to waterfields.geoai.id as a Cloudflare Worker (NODE_OPTIONS makes Node try IPv4 first,
# which this network needs)
deploy:
	NODE_OPTIONS=--dns-result-order=ipv4first npx wrangler deploy
