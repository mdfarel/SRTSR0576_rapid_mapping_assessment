# SRTSR0576 rapid mapping assessment

This repository holds the code of the SRTSR0576 task assignment, a rapid mapping of fields and water bodies from open satellite data. It maps water bodies, agricultural field systems and wet/dry-season change around four archaeological landscapes in Southeast Asia. The mapping, which is rapid and exploratory, uses open satellite data only. The web GIS is available at <https://waterfields.geoai.id>.

| Site | Country | Study area | Seasons used (2024) |
|---|---|---|---|
| Angkor | Cambodia | 44.72 km square around 13.400 N, 103.880 E (1,999.9 km²), used because the supplied KML had no Angkor polygon | dry Feb to Apr, wet mid Aug to Oct |
| Srei Santhor and region | Cambodia | supplied polygon (2,638 km²) | dry Feb to Apr, wet mid Aug to Oct |
| Sukhothai region | Thailand | supplied polygon (2,229 km²) | dry Feb to Apr, wet mid Aug to Oct |
| Prambanan region | Indonesia | supplied polygon (1,380 km²) | wet Jan to Mar, dry Jul to Sep |

Each site is analysed on the UTM grid that encloses its polygon, and every statistic is computed inside that polygon. The question is whether distinctive landscape signatures can be detected from publicly available data at all. The project is a scoping exercise, not a survey.

## Data

| Dataset | Use | Licence |
|---|---|---|
| Sentinel-2 L2A (Microsoft Planetary Computer) | seasonal and monthly median composites, water indices, NDVI, true-colour background | Copernicus open licence |
| Sentinel-1 RTC (Planetary Computer) | monthly VV and VH, seasonal water, paddy cycle | CC-BY 4.0 (Planetary Computer RTC; contains modified Copernicus Sentinel data) |
| JRC Global Surface Water v1.3 (1984 to 2020) | long-term water history, agreement check (water) | Copernicus / EC, open |
| ESA WorldCover 2021 v200 | starting point of the agricultural domain, agreement check (water, cropland) | CC-BY 4.0 |
| Copernicus DEM GLO-30 | slope masks, hillshade, height of farmland below reservoirs | Copernicus DEM licence |
| Fields of The World 2024 (Robinson et al. 2026) | field size and orientation, agreement check (fields) | CC-BY 4.0 |

These are the datasets named in the task brief, together with Fields of The World, which the brief itself suggests. No other dataset enters the analysis, the products, the QGIS projects or the web data. The web GIS opens on the dry-season Sentinel-2 true colour of the analysis, which carries no labels. It also offers Light, Dark and Streets backgrounds (OpenFreeMap vector tiles, © OpenStreetMap contributors) with their own place labels, so that viewers can find their way, but these backgrounds are not used as data.

## Method in brief

- **Water.** Dry-season open water is mapped at 10 m where Sentinel-2 NDWI exceeds an Otsu threshold, which is set on the pixels around water edges. Each water pixel must also be confirmed by MNDWI > 0, and steep slopes are masked. Where the edge sample is not bimodal, a conventional threshold is used instead, and this choice is flagged. Each water body is then described by its length, width, elongation, enclosure and rectilinearity, and these measures sort the bodies into rivers, straight canals, enclosing moats, rectilinear reservoirs and ponds, and irregular lakes and ponds. Enclosure is measured with a ring test. Before the test, gaps of about 60 m are bridged, as are gaps of about 200 m between straight segments, because causeways cut most moats, and a body must border the land that it encloses. Rectilinearity, the four-fold coherence of the outline directions, also gives the deviation of each body from the cardinal axes. Water mapped monthly from Sentinel-1 VV gives the wet-season extent and the number of months with water. Every area is reported with the range that results from varying its thresholds.
- **Fields.** The agricultural domain starts from WorldCover cropland and is extended to grass and wetland pixels that show the paddy radar cycle and to land with two crop cycles a year. Within it, four classes are assigned in order: flood-recession fields, irrigated or multi-cropped land, the rainfed rice-field landscape, and other open agricultural land. Crop cycles are identified as peaks in gap-filled monthly NDVI.
- **Orientation.** The direction of field boundaries is derived from the structure tensor of 10 m NDVI and compared with the orientation of Fields of The World polygons. The share of boundary energy within 7.5° of north-south or east-west is checked in two ways: against edges outside farmland, and by rotating the image by 30°, since a share that stays on the image axes after rotation reflects the pixel grid rather than the ground. Rectilinear water bodies are tested in the same way. A share that fails the rotation test is reported but not interpreted.
- **Change.** Water is mapped from Sentinel-1 VV with one threshold for both seasons, comparing the dry-season median with every wet-season month. Water present in the wet season only is counted as expansion (seasonal inundation), and water present in the dry season only as contraction. Both are measured on this radar basis and are never set against the optical dry-season water. Vegetation change is the Sentinel-2 NDVI difference, wet minus dry, where both seasons have a clear view. Greening and browning, defined as a change of more than 0.15, are reported for the whole study area and for agricultural land alone, in km² and as a share of farmland, because the brief asks for the change associated with cultivation. VH backscatter change is given as a cloud-free complement, and the long-term water history from 1984 to 2020 comes from the JRC transitions.
- **Water sources.** To see which kind of water lies next to today's irrigated land, every agricultural pixel is placed in one of four groups. Farmland in reservoir reach lies within 5 km of a rectilinear reservoir or moat segment of 10 ha or more and at least 1 m below its surface on the Copernicus DEM, so that gravity could carry water to it. Farmland near natural water lies within 1 km of a river, stream, lake or irregular reservoir, or was flooded for at least one month in 2024 (Sentinel-1). The two remaining groups hold farmland with both kinds of source and farmland with neither. For each group, the share that is irrigated or multi-cropped is compared with the share over all farmland. The three distances are then varied together over 27 settings (reach of 3, 5 or 8 km, natural water within 0.5, 1 or 2 km, and a drop of 0, 1 or 3 m), and the number of settings in which reservoir reach ranks above natural water is reported, so that the dependence on these choices can be seen. The result describes where irrigated land sits in 2024. It cannot show what past reservoirs did, because modern canals, pumps and wells are not visible at this resolution and one year may not be typical, so it can frame, but not settle, the long debate over whether large reservoirs supplied rice fields (van Liere 1980). The groups are written to `data/<site>/products/water_sources.json` and the cross-site CSV, and compared across the sites in `out/figures/compare_water_sources.png`.
- **Agreement.** Each result is compared with the datasets of the brief, and the scores are recorded in `data/<site>/products/accuracy_stats.json`. Dry-season water is compared with JRC occurrence and WorldCover water, seasonal inundation with JRC seasonal water, the agricultural domain with WorldCover cropland and Fields of The World, and radar water with Sentinel-2. The comparison with WorldCover cropland is not independent, because the domain starts from it. In any case, agreement between maps is not accuracy. As the next step, a stratified sample of about 300 points per site (Olofsson et al. 2014) has been drawn from the maps in `pipeline/sample.py`. No point has been labelled yet, and the estimators that would score the labels (Stehman 2014) are in `pipeline/sample_accuracy.py`.
- **Interpretation.** Because the assessment is rapid, it interprets the mapped classes only. The maps and statistics name no place, label no individual feature and delineate no special area. A class name such as "Enclosing moat (candidate)" describes a shape, not an identification.

All thresholds are set in the code, next to the step that uses them. Every number on a slide is written by the pipeline to `data/<site>/products/*.json`. The cross-site table is saved as `out/dataset/waterfields_comparative_2024.csv`, and its note column flags the values that are not read: radar water at a site where Sentinel-1 cannot see the water, and orientation shares that fail the rotation test.

## Limits

The imagery shows only today's surface. Buried and forested features, and sub-pixel features such as bunds, minor canals and low embankments, therefore need airborne LiDAR or very high resolution imagery. The shape and orientation classes are morphological candidates and do not date a feature. At Prambanan the radar cannot see the open water in 2024: Sentinel-1 VV over the water mapped from Sentinel-2 stays at about -11 dB, well above the usual open-water range, so radar water, expansion and contraction there are reported as not measurable. The Prambanan farmland total (250.7 km², 18.2% of the area) is a lower bound, because 10 m land cover labels many small sawah and home-garden plots as tree cover or built-up land, and the paddy radar signature fails on such small plots. In the same area, 41 small dark surfaces above 500 m on the upper volcano slope (0.21 km² of the 1.1 km² of dry-season water) may be shadow or volcanic deposit rather than water. Wet-season cloud leaves gaps in the NDVI change, so the share of each area with a clear view is reported with it. On the static maps, water bodies under 1 ha are drawn enlarged so that they can be seen at slide size, whereas the products and statistics keep the mapped outlines. A single year cannot show whether 2024 was typical. Finally, agreement with other maps is not the same as accuracy, and until the validation sample is labelled, the totals are mapped areas rather than error-adjusted estimates.

## Structure

```
run.py        the single entry point: every step for each site, then the comparative dataset
pipeline/     one module per step, called by run.py in order (each can also be run on its own)
webgis/       the web GIS: a static MapLibre page that loads the results that run.py writes into it
config.py     study areas, seasons and grids
```

## Run it

The repository holds the code only. The folders `data/`, `out/`, `webgis/data/` and `webgis/tiles/` are generated by the pipeline and are not committed, so in a fresh clone the web GIS served by `make serve` stays empty until the pipeline has run. The published web GIS, built from the results described here, is at <https://waterfields.geoai.id>.

One file runs the whole process. `python3 run.py` downloads and caches the data, runs the analysis, draws the maps and charts, writes the layers and numbers of the web GIS and builds the cross-site dataset, site by site. The flag `--gis` adds the optional QGIS project of each site. `python3 run.py --help` lists the steps, and `--from` restarts at any of them.

```bash
python3 run.py                                 # all four sites, every step, then the comparative dataset
python3 run.py angkor                          # one site
python3 run.py --from water                    # rerun the analysis only, downloads kept
python3 run.py --gis                          # everything, including the QGIS projects
```

The code needs Python 3.11 or later and the GDAL command-line tools, which are used for the web tiles. QGIS is needed only for the optional QGIS projects. Each job also has a `make` target. The table shows the command that each target runs, and either the target or the command can be typed.

**Tested with.** The published results were produced with Python 3.14.4 on Ubuntu 26.04 under WSL 2, with numpy 2.3.5, scipy 1.18.1, rasterio 1.5.1 (GDAL 3.12), pyproj 3.7.2 (PROJ 9.7), shapely 2.1.2, pystac-client 0.9.0, planetary-computer 1.0.0, duckdb 2.0.0, matplotlib 3.10.7, pmtiles 3.8.1, Pillow 12.2.0 and QGIS 3.40. `requirements-lock.txt` pins these packages, and every package that they depend on, at the versions used, whereas `requirements.txt` gives lower bounds only. The pinned versions need Python 3.12 or later, so on Python 3.11 the packages should be installed from `requirements.txt` instead.

| Target | Runs | Does |
|---|---|---|
| `make install` | `python3 -m pip install -r requirements.txt` | installs the Python packages |
| `make run` | `python3 run.py` | all four sites, then the cross-site dataset |
| `make site SITE=angkor` | `python3 run.py angkor` | one site |
| `make all` | `python3 run.py --gis` | everything, including the QGIS projects |
| `make gis` | `python3 -m pipeline.gis_project` | QGIS project per site (needs a Python with PyQGIS) |
| `make dataset` | `python3 run.py --from "comparative dataset"` | cross-site CSV and comparison charts only |
| `make serve` | `python3 serve.py 8765` | the web GIS on this computer |
| `make deploy` | `NODE_OPTIONS=--dns-result-order=ipv4first npx wrangler deploy` | owner's step: publishes `webgis/` to waterfields.geoai.id |

The packages can also be kept in a virtual environment, which is needed where the system Python refuses `pip install`. To use one, create it once with `python3 -m venv .venv` and add its Python to each target, for example `make install PYTHON=.venv/bin/python`.

**Time and disk.** Measured in decimal units (1 GB = 10⁹ bytes) after the published run, `data/` holds about 23.5 GB for the four sites (Angkor 5.7 GB, Srei Santhor 8.2 GB, Sukhothai 6.1 GB and Prambanan 3.4 GB), most of it cached Sentinel-1 and Sentinel-2 reads. The outputs add about 1.0 GB in `out/`, mainly the maps and charts (0.73 GB), the QGIS projects (0.10 GB) and their zips (0.07 GB), and the web GIS adds 53 MB of tiles in `webgis/tiles/`. A full run therefore needs about 25 GB of free disk. The run time is dominated by the downloads, so it depends mainly on the speed of the connection to Microsoft Planetary Computer, which serves the Sentinel, JRC, WorldCover and DEM layers, and to Source Cooperative, which serves Fields of The World. Once the downloads were cached, rerunning the analysis for all four sites with `python3 run.py --from water` took about 13 minutes on the author's computer.

**Deploy.** `make deploy` is the owner's step. It publishes `webgis/` to waterfields.geoai.id and needs the owner's Cloudflare login, so reviewers do not need to run it.

Downloads are cached under `data/` and resume after an interruption. The outputs are organised as follows:

```
data/<site>/products/          GeoTIFF products and *.json statistics
out/figures/<site>/            maps and charts
out/dataset/                   cross-site comparative CSV
out/gis/<site>/                QGIS project (<site>.qgz) with all of its layers
webgis/data/, webgis/tiles/    the results as the web GIS loads them (JSON, CSV and PMTiles)
```

Each folder in `out/gis/` is a self-contained QGIS project. The project file refers to its layers by relative paths, so the folder can be zipped, sent and opened in QGIS on another computer. It holds the open data of the brief, clipped to the study area (Sentinel-2 true colour and Sentinel-1 VV for both seasons, JRC occurrence, WorldCover and a DEM hillshade), and the water, field and change layers in the colours of the maps. It also holds the water-body polygons and table, the Fields of The World fields (as polygons at Prambanan and as a 10 m outline raster at the other sites, which keeps each zip small), the study area and the comparative CSV.

`make serve` supports the byte-range requests that PMTiles needs. Once the server is running, the web GIS is at <http://localhost:8765>.

## Layout

```
run.py                     runs every step in order (python3 run.py --help lists them)
config.py                  study areas, seasons, grids
pipeline/common.py         grid, polygons, STAC access, cached reads
pipeline/layers.py         per-site access to the cached layers, indices and thresholds
pipeline/fetch_*.py        static layers, Sentinel-2 composites and monthly NDVI, Sentinel-1 monthly
pipeline/water.py          water bodies, form classes, seasonal water, threshold sensitivity
pipeline/phenology.py      crop cycles from monthly NDVI
pipeline/fields.py         field systems, water relation, orientation and its controls
pipeline/change.py         wet against dry season (water and vegetation), JRC long-term transitions
pipeline/sources.py        water sources next to irrigated farmland (reservoir reach, natural water) and its 27-setting sweep
pipeline/ftw.py            Fields of The World polygons inside each study area
pipeline/accuracy.py       agreement with the datasets of the brief
pipeline/sample.py         stratified validation sample (design only)
pipeline/sample_accuracy.py  accuracy and error-adjusted areas, once labels exist
pipeline/maps.py, tiles.py, dataset.py, gis_project.py   outputs
webgis/                    MapLibre GL JS page that loads webgis/data/ and webgis/tiles/
serve.py, worker.js        local server and Cloudflare Worker, both answering byte-range requests
Makefile                   the commands above, one target each
```

## How it was built

I designed the questions, the data choices, the thresholds and the checks, and judged the results. The code was written with AI coding assistance under my direction. It was then run, tested and checked against the data and the literature.

## Citation

Citation details are in `CITATION.cff`. The code is licensed under the MIT licence, and derived data inherit the licences of their sources.

Muhammad Dita Farel
