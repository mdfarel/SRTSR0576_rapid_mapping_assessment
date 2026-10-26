# Classification reference

This page lists, in one place, how every water, field and seasonal class is assigned: the bands read, the indices computed from them, the thresholds and the order of the rules. It also maps each item of the task brief to the class or statistic that answers it. The rules themselves live in the code. Every threshold below names the constant (or, where the value is written inline, the function) that holds it, so that the page can be checked against the code and does not drift from it. Thresholds that are set from the data, such as the Otsu values, differ by site and are read from `data/<site>/products/*.json` rather than repeated here.

No class is learned from training data. Every class is assigned by explicit thresholds and shape rules. Two inputs are themselves products of trained models: ESA WorldCover, where the agricultural domain starts, and Fields of The World, which is used for checks and for field size and orientation only.

## 1. Inputs

| Input | Bands or values used | Processing | Code |
|---|---|---|---|
| Sentinel-2 L2A, seasonal composites | B03 green, B04 red, B08 near infrared at 10 m; B02 blue, B03, B04, B08, B11 and B12 shortwave infrared at 20 m | Scenes up to 60% (dry) or 90% (wet) cloud, the clearest per MGRS tile; pixels kept where the scene classification (SCL) is 2, 4, 5, 6, 7 or 11; the processing-baseline offset of -1000 removed; per-pixel median per season | `config.S2_BANDS_10`, `config.S2_BANDS_20`, `config.S2_VALID_SCL`, `fetch_s2.MAX_CLOUD`, `fetch_s2.BOA_OFFSET` |
| Sentinel-2 L2A, monthly NDVI | B04, B08 at 20 m | Up to 4 clearest scenes per tile and month (80% cloud or less), same cloud mask, per-pixel median NDVI per month | `fetch_s2_monthly.PER_TILE`, `fetch_s2_monthly.MAX_CLOUD` |
| Sentinel-1 RTC | VV and VH gamma0 at 20 m | Every 2024 acquisition, monthly median in dB | `fetch_s1` |
| Copernicus DEM GLO-30 | Elevation | Slope in degrees | `layers.slope_deg` |
| ESA WorldCover 2021 | Class codes 20 shrub, 30 grass, 40 cropland, 50 built-up, 80 permanent water, 90 herbaceous wetland | Share of each 20 m cell in each class | `fields.WC_*` |
| JRC Global Surface Water v1.3 | Occurrence, seasonality, transitions | Checks of the water maps; long-term transitions | `accuracy`, `change.GROUPS` |
| Fields of The World 2024 | Field polygons | Fields of 2 ha or more and 1.5 times longer than wide, for orientation | `ftw` |

Seasons (`config.MAINLAND_SEASONS`, `config.JAVA_SEASONS`), all in 2024 (`config.YEAR`):

| Sites | Dry season | Wet season |
|---|---|---|
| Angkor, Srei Santhor, Sukhothai | 1 February to 30 April | 15 August to 31 October |
| Prambanan | 1 July to 30 September | 1 January to 31 March |

## 2. Indices

| Index | Formula (Sentinel-2 bands) | Why | Used for |
|---|---|---|---|
| NDWI (McFeeters 1996) | (B03 - B08) / (B03 + B08) | Water reflects green and absorbs near infrared | Dry-season open water |
| MNDWI (Xu 2006) | (B03 - B11) / (B03 + B11) | Shortwave infrared is absorbed even more strongly by water; rejects built-up surfaces that pass NDWI | Confirmation of open water |
| NDVI | (B08 - B04) / (B08 + B04) | Green vegetation reflects near infrared and absorbs red | Crop cycles, dry-season greenness, seasonal change, field edges |

A normalised difference is left undefined where the two reflectances sum to 0.02 or less (`layers.MIN_REFLECTANCE_SUM`).

Radar values are read directly. Calm open water reflects the radar pulse away from the sensor, so VV is low over water. VH responds to vegetation volume: a flooded paddy is dark, and VH rises as the crop grows.

## 3. Thresholds set from the data

Both water thresholds are Otsu thresholds computed on the pixels within 2 pixels of a first water-land edge (`layers.edge_otsu`, after Donchyts et al. 2016 and Markert et al. 2020), so that water and land are sampled in similar numbers even where water covers a few percent of the area.

| Threshold | First split | Held within | Fallback when the edge sample is not bimodal (eta under 0.3) | Set in | Value per site |
|---|---|---|---|---|---|
| NDWI, dry season | 0.0 | -0.10 to 0.25 | 0.0 | `water.run` | `water_stats.json`: `thresholds.ndwi_edge_otsu` |
| Sentinel-1 VV | -18 dB | -24 to -13 dB | -18 dB | `water.s1_water` | `water_stats.json`: `thresholds.s1_vv_db_edge_otsu` |

How each value was reached (Otsu, Otsu held at a bound, or the fallback) is written to `water_stats.json` under `threshold_method`.

## 4. Water bodies

### 4.1 Water pixels

| Layer | Rule | Constants |
|---|---|---|
| Dry-season open water, 10 m | NDWI above its Otsu threshold, MNDWI above 0, slope under 15°, connected bodies of 10 pixels (0.1 ha) or more | `water.water_mask`, `MAX_SLOPE_S2`, `MIN_PX` |
| Radar water, 20 m, per month | VV below its Otsu threshold, slope under 10° | `water.s1_water`, `MAX_SLOPE_S1` |
| Months with water | Count of months with radar water (0 to 12) | `s1_water_months_20m` |
| Dry-season radar water | Median VV over the dry-season months below the threshold | `s1_water_dry_20m` |
| Wet-season radar water | VV below the threshold in any wet-season month | `s1_water_wet_20m` |
| Radar not measurable | Median monthly VV over the optically mapped water above -15 dB | `maps.RADAR_BLIND_DB`, `maps.radar_cannot_see_water` |

### 4.2 Shape measures

Each 8-connected body of dry-season open water is described by (`water.describe`):

- area, and length and breadth along and across its principal axis (second moments of the pixel positions);
- mean width (area over length) and elongation (length over breadth);
- rectangularity, the share of its oriented bounding rectangle that it fills;
- rectilinearity, the four-fold coherence of its outline edge directions: 1 when every edge follows one pair of perpendicular directions, near 0 for a rounded or meandering outline (`water.outline_rectilinearity`);
- enclosure, from a ring test: water closed over gaps of about 60 m (and straight segments of 1 ha or more over gaps of about 200 m, for causeways), scored by the area of the filled ring over the ring's own area (`water.enclosure_ratio`);
- deviation of its axis pair from north-south and east-west, corrected for grid convergence.

### 4.3 Form classes

The first rule that applies wins (`water.classify`). "Straight" means rectilinearity of 0.40 or more (`RECTILINEAR_R`).

| Order | Class | Rule | Constants |
|---|---|---|---|
| 1 | Small water body (shape not resolved) | Under 20 pixels (0.2 ha) | `SHAPE_MIN_PX` |
| 2 | Enclosing moat (candidate) | Enclosure 1.5 or more, straight, mean width 300 m or less | `MOAT_MIN_ENCLOSED`, `MOAT_MAX_WIDTH_M` |
| 3 | Straight canal (candidate) | Elongation 5 or more, length 500 m or more, straight, width 60 m or less | `CHANNEL_MIN_ELONGATION`, `CHANNEL_MIN_LENGTH_M`, `CANAL_MAX_WIDTH_M` |
| 3 | River or stream | Elongation 5 or more and length 500 m or more, but not straight or wider than 60 m | as above |
| 4 | River or stream | Length 1 km or more, not straight, rectangularity under 0.35 | `RIVER_MIN_LENGTH_M`, `RIVER_MAX_RECTANGULARITY` |
| 5 | Rectilinear reservoir or moat segment (10 ha or more) | Straight, 1,000 pixels (10 ha) or more | `LARGE_PX` |
| 5 | Rectilinear pond | Straight, under 10 ha | `LARGE_PX` |
| 6 | Lake or irregular reservoir | Not straight, 10 ha or more | `LARGE_PX` |
| 6 | Irregular pond | Not straight, under 10 ha | `LARGE_PX` |

A straight-sided body (canal, moat, reservoir or moat segment, rectilinear pond) is cardinal when its axis pair lies within 7.5° of north-south and east-west (`CARDINAL_TOL`). If directions were random, 16.7% would be.

Ponds for the density map and counts are the rectilinear pond, irregular pond and small water body classes (`POND_CLASSES`). Large ponds are rectilinear and irregular ponds of 1 to 10 ha, counted from `water_features.csv` (`dataset.count_large_ponds`).

### 4.4 Orientation control for water

A central window of the dry-season water is rotated by 30° and classified again (`water.rotation_control`). The cardinal share passes when, after rotation, the share on the rotated ground axes exceeds the share on the image axes by more than 10 points. It is not judged with fewer than 30 rectilinear bodies in the window. A share that fails or is not judged is reported but not interpreted.

## 5. Agricultural field systems

### 5.1 Crop cycles

From the twelve monthly NDVI values (`phenology`):

1. months without a clear view are filled by linear interpolation between the nearest clear months (`fill_gaps`);
2. the series is smoothed 1-2-1 over neighbouring months, with December next to January;
3. a peak counts as a crop cycle when it reaches NDVI 0.45 and rises at least 0.15 above the lowest value in the three months before and in the three months after it (`PEAK_MIN`, `PROMINENCE`, `WINDOW`);
4. pixels with fewer than 6 clear months are not judged (`MIN_CLEAR_MONTHS`).

The amplitude is the range of the smoothed series over the year.

### 5.2 Agricultural domain (20 m)

A 20 m cell is agricultural when its slope is under 20° (`fields.MAX_SLOPE`), it lies inside the study area, and any of the following holds (`fields.classify`):

- half or more of it is WorldCover cropland;
- half or more of it is WorldCover grass or herbaceous wetland and it shows the paddy radar signature (5.3);
- it has two or more crop cycles with an amplitude of 0.30 or more (`AMP_CROP`), or one cycle with that amplitude over half or more shrub, grass or wetland, and under half of it is built-up or permanent water.

### 5.3 Field classes

Assigned in this order of precedence, so each pixel takes the first class it qualifies for.

| Order | Class | Rule | Constants |
|---|---|---|---|
| 1 | Flood-recession fields | Radar water in 2 or more wet-season months, water in 9 months or fewer of the year, a clear dry-season view, dry-season NDVI 0.35 or more | `NDVI_RECESSION` |
| 2 | Irrigated or multi-cropped | Dry-season NDVI 0.45 or more with amplitude 0.25 or more and at least one judged crop cycle; or two or more crop cycles | `NDVI_IRRIGATED`, `AMP_IRRIGATED` |
| 3 | Rice-field landscape (rainfed) | Paddy signature: lowest monthly VH -19 dB or less, VH rising 5 dB or more over the year, radar water in 6 months or fewer | `RICE_VH_MIN`, `RICE_VH_RANGE` |
| 4 | Open agricultural land | The rest of the domain | |

Because irrigated rice qualifies for class 2 first, class 3 holds rainfed rice only, and all rice is larger than class 3. The paddy signature is flagged as failing at a site where it covers under 1% of the mapped farmland (`dataset.PADDY_BLIND_SHARE`).

Sensitivity: the paddy signature area is recomputed for a VH low of -18 to -21 dB and a rise of 5, 6.5 or 8 dB, and dry-season greenness for NDVI 0.35, 0.45 and 0.55 (`field_stats.json`: `sensitivity`).

### 5.4 Field statistics

| Statistic | Rule | Code |
|---|---|---|
| Distance to perennial water | Distance from each pixel to dry-season water, either optical (any of the four 10 m pixels of a cell) or radar | `fields.water_relation` |
| Relation of fields to water | Share of each class within 500 m and 1 km, against the share of all land | `field_stats.json`: `relation_to_perennial_water` |
| Concentrated agricultural activity | 80% or more of a 1 km moving window is agricultural | `fields.run`, `concentrated_ag_km2` |
| Cropping intensity | Share of judged farmland with 0, 1, 2 or 3 crop cycles | `cropping_cycles_share_pct` |

### 5.5 Field-edge orientation

The structure tensor of the 10 m dry-season NDVI (Sobel gradients, Gaussian smoothing of 2 pixels) gives each pixel a dominant edge direction and a coherence (`fields.tensor`). Edges with coherence 0.5 or more (`MIN_COHERENCE`) in the stronger half of edge energy are kept, and the energy-weighted share within 7.5° of the cardinal axes is reported with a 95% interval from resampling 2 km blocks (`block_bootstrap`). Two controls follow:

- the same share for edges outside farmland;
- a central window rotated by 30°: the share passes when, after rotation, the share on the image axes falls below the share before rotation (`dataset.passes_rotation_test`).

## 6. Seasonal change

| Measure | Data | Rule | Code |
|---|---|---|---|
| Water in both seasons | Sentinel-1 VV | Radar water in the dry-season median and in at least one wet-season month | `change.run`, class 1 |
| Expansion (wet season only) | Sentinel-1 VV | Radar water in a wet-season month but not in the dry-season median | class 2 |
| Contraction (dry season only) | Sentinel-1 VV | Radar water in the dry-season median but in no wet-season month | class 3 |
| Vegetation change | Sentinel-2 NDVI | Wet-season composite minus dry-season composite, where both have a clear view; greener above +0.15, browner below -0.15, strongly so beyond ±0.3 | `change.NDVI_CHANGE`, `maps.NDVI_STRONG` |
| Vegetation change on farmland | as above | The same, inside the agricultural domain, as km² and as a share of farmland | `change.farmland_ndvi_change` |
| Radar vegetation change | Sentinel-1 VH | Wet-season median minus dry-season median, blank over radar water; rise above 3 dB reported on farmland | `change.VH_RISE_DB` |
| Long-term water | JRC transitions | Water since 1984, new, lost and ephemeral | `change.trajectories` |

The same VV threshold is used in both seasons, and expansion and contraction are never set against the optical dry-season water.

## 7. Sensitivity ranges

| Area | Variation | Field in the products |
|---|---|---|
| Dry-season open water | NDWI threshold ±0.05, MNDWI threshold -0.1, 0 and +0.1 | `water_stats.json`: `sensitivity.dry_season_water_km2_range` |
| Wet-season radar water, expansion | VV threshold ±1 dB | `sensitivity.s1_wet_max_water_km2_range`, `seasonal_inundation_km2_range` |
| Paddy signature, dry-season greenness | see 5.3 | `field_stats.json`: `sensitivity` |

## 8. Items of the brief

### 8.1 Water bodies and hydraulic features

| Brief | Class or measure | Gaps |
|---|---|---|
| Rivers | River or stream | Rivers narrower than about 20 m break into pieces or fall into the pond classes |
| Lakes | Lake or irregular reservoir | Lakes and reservoirs are separated by shape, not by origin |
| Reservoirs | Rectilinear reservoir or moat segment, and Lake or irregular reservoir | A reservoir that is dry in the dry season is not mapped |
| Large ponds | Rectilinear and irregular ponds of 1 to 10 ha (count, no map class) | Ponds under 0.1 ha are not mapped and many under 0.5 ha are missed |
| Moats | Enclosing moat (candidate) | A dry, vegetated or broken moat fails the ring test and falls into moat segments or ponds |
| Major canals | Straight canal (candidate) | Canals narrower than about 20 m or dry in the dry season are missed |
| Total water area | Dry-season open water; wet-season radar water | |
| Density and distribution | Bodies and ponds per 100 km²; pond density map | Lower bounds |
| Strong seasonal fluctuation | Expansion and contraction; months with water | Not measurable where radar does not see the water |

All form classes are read from dry-season optical water only.

### 8.2 Agricultural and field systems

| Brief | Class or measure | Gaps |
|---|---|---|
| Open agricultural land | Open agricultural land | |
| Rice-field landscapes | Rice-field landscape (rainfed) | Rainfed rice only; irrigated rice is in the next row |
| Irrigated cultivation zones | Irrigated or multi-cropped | Inferred from dry-season greenness or multiple cycles, not from irrigation works |
| Other field systems | Flood-recession fields; field-edge orientation; Fields of The World field size | No class for field layouts; bunds and plot patterns need very high resolution imagery or LiDAR |
| Total agricultural area | Agricultural domain | A lower bound where small plots are labelled tree cover or built-up |
| Relation between fields and water | Distance to perennial water by class | Numbers only on the deck |
| Concentrated agricultural activity | 80% or more of a 1 km window | Numbers only on the deck |

### 8.3 Annual and seasonal change

| Brief | Class or measure | Gaps |
|---|---|---|
| Expansion and contraction of water bodies | Radar water change classes 2 and 3 | Monthly medians damp short floods |
| Changes in vegetation cover associated with cultivation | NDVI change on farmland; VH change | Wet-season cloud leaves gaps; the clear share is reported |
| Areas affected by seasonal inundation | Expansion (class 2); flood-recession fields | Reported as the same radar area as expansion |
