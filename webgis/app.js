/* Web map of the SRTSR0576 task assignment: MapLibre GL JS + PMTiles, static, no build step.

   Reads, from the same site: data/sites.json (the published sites), data/<site>.json (statistics and
   study-area outline), data/<site>_water_features.csv (one row per water body) and tiles/<site>/<layer>.pmtiles
   (one raster archive per layer, written by pipeline/tiles.py). The statistics panel also links to
   data/waterfields_comparative_<year>.csv, the cross-site table.
   Preview locally with `python serve.py` (or `make serve`). */
"use strict";

/* Labels are the study-area names of config.py; the page names no other place. */
const SITES = [
  { key: "angkor", label: "Angkor, Cambodia", epsg: 32648 },
  { key: "srei_santhor", label: "Srei Santhor and region, Cambodia", epsg: 32648 },
  { key: "sukhothai", label: "Sukhothai region, Thailand", epsg: 32647 },
  { key: "prambanan", label: "Prambanan region, Indonesia", epsg: 32749 },
];

const ACCENT = "#c2410c";
const SEASONAL = "#86b6e6";  // seasonal inundation (wet season only) in the water layer
const FTW_URL = "https://data.source.coop/ftw/global-data/predictions/vectors/alpha/2024_with_confidence.pmtiles";

/* The layer list, in panel order. Every raster is drawn in the exact colours of its legend (pipeline/tiles.py),
   so the identify tool names a pixel by its legend colour. `radar` marks a layer drawn from Sentinel-1 water:
   "all" for the whole layer and "seasonal" for its seasonal-inundation class alone, which the page flags at a
   site where the radar cannot see the water. */
const LAYERS = [
  { id: "water", group: "Analysis", radar: "seasonal", title: "Water bodies by form", sub: "Sentinel-2 10 m dry season; Sentinel-1 wet season", on: true, legend: [
    ["#3b7fc4", "River or stream"], ["#1f4e89", "Lake or irregular reservoir"], ["#0b7a75", "Straight canal (candidate)"],
    [ACCENT, "Enclosing moat (candidate)"], ["#14a39c", "Rectilinear reservoir or moat, 10 ha+"], ["#5cc8c2", "Rectilinear pond"],
    ["#6a3d9a", "Irregular pond"], ["#a75cc0", "Small water body, shape not resolved"], [SEASONAL, "Seasonal inundation, wet only"]] },
  { id: "cardinal", group: "Analysis", title: "Orientation of rectilinear water", sub: "Within 7.5 degrees of north-south or east-west", on: false, legend: [
    [ACCENT, "Rectilinear, cardinally aligned"], ["#3a5f7d", "Rectilinear, other orientation"], ["#9fb8cf", "Other water bodies"]] },
  { id: "density", group: "Analysis", title: "Density of ponds and small water bodies", sub: "Pond-class bodies (all under 10 ha) per km², 1 km window", on: false, legend: [
    ["#c7e9b4", "1 per km²"], ["#7fcdbb", "2 to 3 per km²"], ["#1d91c0", "4 to 7 per km²"], ["#081d58", "8 or more per km²"]] },
  { id: "fields", group: "Analysis", title: "Agricultural field systems", sub: "Sentinel-1 paddy cycle, Sentinel-2 dry-season NDVI, WorldCover", on: false, legend: [
    ["#178c96", "Flood-recession fields"], ["#2e8540", "Irrigated or multi-cropped"],
    ["#b5c64b", "Rice-field landscape (rainfed)"], ["#e8dccb", "Open agricultural land"]] },
  { id: "cycles", group: "Analysis", title: "Cropping intensity", sub: "Crop cycles in 2024 from monthly Sentinel-2 NDVI", on: false, legend: [
    ["#b5c64b", "One crop cycle"], ["#2e8540", "Two crop cycles"], ["#0f4d2a", "Three crop cycles"], ["#a3876a", "No clear crop cycle"]] },
  { id: "grain", group: "Analysis", title: "Field-boundary grain", sub: "Dominant edge direction over cropland, 1 km window", on: false, legend: [
    [ACCENT, "Within 7.5 degrees of cardinal"], ["#e8a27c", "7.5 to 15 degrees"], ["#c9c3bc", "15 to 30 degrees"], ["#8e9aaf", "30 to 45 degrees"]] },
  { id: "months", group: "Seasonal change", radar: "all", title: "Months with open water", sub: "Sentinel-1 monthly medians, 2024", on: false, legend: [
    ["#c6dbef", "1 to 2 months"], ["#6baed6", "3 to 5 months"], ["#2171b5", "6 to 9 months"], ["#08306b", "10 to 12 months"]] },
  { id: "water_change", group: "Seasonal change", radar: "all", title: "Water, dry versus wet season", sub: "Sentinel-1 VV monthly medians", on: false, legend: [
    ["#1f4e89", "Water in both seasons"], [SEASONAL, "Seasonal inundation, wet only"], ["#e07b39", "Water in dry season only"]] },
  { id: "ndvi_change", group: "Seasonal change", title: "Vegetation, wet minus dry season", sub: "Sentinel-2 NDVI difference; gaps where the wet season stayed cloudy", on: false, legend: [
    ["#01665e", "Strongly greener, +0.3 or more"], ["#80cdc1", "Greener, +0.15 to +0.3"], ["#f5f5f5", "Little change, −0.15 to +0.15"],
    ["#dfc27d", "Browner, −0.15 to −0.3"], ["#8c510a", "Strongly browner, −0.3 or less"]] },
  { id: "ftw", group: "Reference", vector: true, ink: "ftw", title: "Fields of The World 2024", sub: "Predicted modern field boundaries, CC-BY-4.0; zoom 11 and closer", on: false, legend: [
    ["#2f6f3a", "Predicted field boundary"]] },
  { id: "jrc", group: "Reference", title: "JRC water occurrence", sub: "Share of 1984 to 2020 observations with water", on: false, legend: [
    ["#c6dbef", "Under 25% of observations"], ["#6baed6", "25 to 49%"], ["#2171b5", "50 to 74%"], ["#08306b", "75% or more"]] },
];

/* Order in which the layers are drawn, bottom first: the layers that cover whole areas lie under the water
   bodies, and the Fields of The World lines lie on top, so that switching on an area layer hides no outline. */
const DRAW_ORDER = ["ndvi_change", "jrc", "months", "water_change", "density", "cycles", "fields", "grain", "water", "cardinal", "ftw"];

/* Note shown with a radar layer at a site where the radar cannot see the water (see radarBlind). */
const RADAR_NOTE = { all: "Not measurable at this site: radar does not see the water here",
  seasonal: "Seasonal inundation not shown at this site: radar does not see the water here" };

const BASEMAPS = [
  { id: "light", label: "Light", style: "https://tiles.openfreemap.org/styles/positron" },
  { id: "dark", label: "Dark", style: "https://tiles.openfreemap.org/styles/dark" },
  { id: "streets", label: "Streets", style: "https://tiles.openfreemap.org/styles/liberty" },
  { id: "imagery", label: "Sentinel-2", style: null },
];
/* The page opens on the Sentinel-2 dry-season true colour that the analysis itself uses. */
const DEFAULT_BASEMAP = "imagery";

/* Width at which the panels lie over the bottom of the map, one at a time (the max-width of the narrow layout in
   style.css), and the height under which a narrow screen, such as a phone held sideways, opens on the map alone. */
const NARROW = matchMedia("(max-width: 1023px)");
const SHORT = matchMedia("(max-height: 450px)");

const $ = (s) => document.querySelector(s);

/* A number with d decimals and thousands separators, or "n/a" when missing. */
function fmt(v, d = 1) {
  if (v === null || v === undefined) return "n/a";
  return Number(v).toLocaleString("en-AU", { maximumFractionDigits: d, minimumFractionDigits: d });
}

/* Remembered choices (basemap, theme); the page works without them when storage is blocked. */
const store = {
  get(k, d) { try { return localStorage.getItem(k) ?? d; } catch (e) { return d; } },
  set(k, v) { try { localStorage.setItem(k, v); } catch (e) { /* storage unavailable */ } },
};

const siteInfo = (key) => SITES.find((s) => s.key === key);

const state = {
  site: null, data: {}, features: {}, styleReady: false, pick: null, popup: null, identifyRun: 0,
  on: Object.fromEntries(LAYERS.map((l) => [l.id, l.on])),
  opacity: Object.fromEntries(LAYERS.map((l) => [l.id, 0.9])),
  basemap: store.get("waterfields-basemap", DEFAULT_BASEMAP),
};

const savedTheme = store.get("waterfields-theme", "auto");
if (savedTheme !== "auto") document.documentElement.dataset.theme = savedTheme;

/* The chosen basemap, or the default when the remembered choice is not offered (any more). */
function basemapId() {
  return BASEMAPS.some((b) => b.id === state.basemap) ? state.basemap : DEFAULT_BASEMAP;
}
function styleUrl() {
  const b = BASEMAPS.find((x) => x.id === basemapId());
  // The Sentinel-2 basemap is a raster layer drawn over the light vector style under either theme, with the
  // style's place labels hidden (applyLabels), so that the default view carries no place name.
  return b.style || BASEMAPS[0].style;
}

/* Visibility of each label (symbol) layer of the loaded vector style, as the style itself sets it. */
let labelDefaults = {};
function rememberLabels() {
  labelDefaults = {};
  for (const l of map.getStyle().layers || []) {
    if (l.type === "symbol") labelDefaults[l.id] = (l.layout && l.layout.visibility) || "visible";
  }
}
/* Hide the vector style's labels under the Sentinel-2 basemap, and show them as the style sets them under the
   Light, Dark and Streets backgrounds, which keep their own labels to help viewers find their way. */
function applyLabels() {
  const imagery = basemapId() === "imagery";
  for (const [id, vis] of Object.entries(labelDefaults)) {
    if (map.getLayer(id)) map.setLayoutProperty(id, "visibility", imagery ? "none" : vis);
  }
}
/* Ink for the outline and field boundaries, from the basemap rather than the page theme:
   light on the Dark and Sentinel-2 basemaps, dark on Light and Streets. */
function mapInk() {
  return ["dark", "imagery"].includes(basemapId())
    ? { ink: "#ececee", halo: "#111", ftw: "#8fd19e" }
    : { ink: "#1c1c1e", halo: "#fff", ftw: "#2f6f3a" };
}

/* PMTiles archives are shared between the map and the identify tool */
const protocol = new pmtiles.Protocol();
maplibregl.addProtocol("pmtiles", protocol.tile);
const archives = {};
function archive(url) {
  if (!archives[url]) {
    archives[url] = new pmtiles.PMTiles(url);
    protocol.add(archives[url]);
  }
  return archives[url];
}
const tileUrl = (key, id) => `${location.origin}/tiles/${key}/${id}.pmtiles`;

const map = new maplibregl.Map({
  container: "map",
  style: styleUrl(),
  center: [103.88, 13.4],
  zoom: 9,
  attributionControl: false,  // credits live in the status bar, so they never cover the map
});
map.addControl(new maplibregl.NavigationControl({ showCompass: false }), "top-right");
map.addControl(new maplibregl.ScaleControl({ unit: "metric" }), "bottom-right");

/* The first layer of the basemap's closing run of symbol layers, so that the analysis layers lie over every road,
   building and boundary and under the place labels. The dark and streets styles draw roads and bridges after
   their first label layer, so the first symbol layer of any kind is too low. */
function labelsStartId() {
  const layers = map.getStyle().layers || [];
  let i = layers.length;
  while (i > 0 && layers[i - 1].type === "symbol") i -= 1;
  return i < layers.length ? layers[i].id : undefined;
}

/* Add the site's layers to the map, replacing those of the previous site. Raster layer r-<id> reads source s-<id>. */
function mountSite() {
  const key = state.site;
  if (!key || !state.styleReady) return;
  for (const id of ["truecolor", ...LAYERS.map((l) => l.id)]) {
    if (map.getLayer(`r-${id}`)) map.removeLayer(`r-${id}`);
    if (map.getSource(`s-${id}`)) map.removeSource(`s-${id}`);
  }
  for (const id of ["aoi-casing", "aoi-line", "pick"]) {
    if (map.getLayer(id)) map.removeLayer(id);
  }
  for (const id of ["aoi", "pick"]) {
    if (map.getSource(id)) map.removeSource(id);
  }

  const before = labelsStartId();  // analysis layers go over the basemap's roads and under its place labels
  const tc = tileUrl(key, "truecolor");
  archive(tc);
  map.addSource("s-truecolor", { type: "raster", url: `pmtiles://${tc}`, tileSize: 256 });
  map.addLayer({ id: "r-truecolor", type: "raster", source: "s-truecolor",
    layout: { visibility: basemapId() === "imagery" ? "visible" : "none" }, paint: { "raster-resampling": "linear" } }, before);

  for (const l of DRAW_ORDER.map((id) => LAYERS.find((x) => x.id === id))) {
    if (l.vector) {
      map.addSource(`s-${l.id}`, { type: "vector", url: `pmtiles://${FTW_URL}`, attribution: "Fields of The World (Robinson et al. 2026), CC-BY-4.0" });
      map.addLayer({ id: `r-${l.id}`, type: "line", source: `s-${l.id}`, "source-layer": "2024", minzoom: 11,
        layout: { visibility: state.on[l.id] ? "visible" : "none" },
        paint: { "line-color": mapInk().ftw, "line-width": 0.7, "line-opacity": state.opacity[l.id] } }, before);
      continue;
    }
    const url = tileUrl(key, l.id);
    archive(url);
    map.addSource(`s-${l.id}`, { type: "raster", url: `pmtiles://${url}`, tileSize: 256 });
    map.addLayer({ id: `r-${l.id}`, type: "raster", source: `s-${l.id}`,
      layout: { visibility: state.on[l.id] ? "visible" : "none" },
      paint: { "raster-opacity": state.opacity[l.id], "raster-resampling": "nearest" } }, before);
  }
  if (state.data[key]) mountOutline(state.data[key]);
  // ring around the point last clicked with the identify tool, kept across basemap changes
  map.addSource("pick", { type: "geojson", data: pickData() });
  map.addLayer({ id: "pick", type: "circle", source: "pick", paint: { "circle-radius": 6, "circle-color": "rgba(0,0,0,0)", "circle-stroke-color": ACCENT, "circle-stroke-width": 2 } });
}
const pickData = () => (state.pick
  ? { type: "Feature", geometry: { type: "Point", coordinates: state.pick } }
  : { type: "FeatureCollection", features: [] });

/* Study-area boundary: dashed, in the basemap's ink, on a casing in the halo colour so that it shows against
   imagery and basemap alike. The page draws no other outline, label or point of its own. */
function mountOutline(d) {
  const { ink, halo } = mapInk();
  map.addSource("aoi", { type: "geojson", data: { type: "Feature", geometry: { type: "Polygon", coordinates: [d.boundary] } } });
  map.addLayer({ id: "aoi-casing", type: "line", source: "aoi", paint: { "line-color": halo, "line-width": 3, "line-opacity": 0.7 } });
  map.addLayer({ id: "aoi-line", type: "line", source: "aoi", paint: { "line-color": ink, "line-width": 1.2, "line-dasharray": [3, 2] } });
}

/* Recolour the outline and field boundaries after a basemap change that keeps the same style. */
function applyInk() {
  const { ink, halo, ftw } = mapInk();
  const paint = (id, prop, v) => { if (map.getLayer(id)) map.setPaintProperty(id, prop, v); };
  paint("aoi-casing", "line-color", halo);
  paint("aoi-line", "line-color", ink);
  paint("r-ftw", "line-color", ftw);
}

/* Swatches drawn in the basemap's ink (the Fields of The World line) take the colour of the current basemap. */
function inkSwatches() {
  const ink = mapInk();
  for (const i of document.querySelectorAll("i[data-ink]")) i.style.background = ink[i.dataset.ink];
}

/* layer list and legend */
const chevron = '<svg viewBox="0 0 16 16" aria-hidden="true"><path d="M4 6l4 4 4-4" fill="none" stroke="currentColor" stroke-width="1.5"/></svg>';
/* Legend rows. A layer drawn in the basemap's ink names its key in mapInk() as `ink`; its swatch shows that colour
   and carries data-ink, so that inkSwatches() can update it after a basemap change. */
const swatches = (legend, ink) => `<div class="swatches">${legend.map(([c, t]) => (ink
  ? `<div><i data-ink="${ink}" style="background:${mapInk()[ink]}"></i>${t}</div>`
  : `<div><i style="background:${c}"></i>${t}</div>`)).join("")}</div>`;

/* True when the radar cannot see the water at the site shown (see radarBlind). */
const siteBlind = () => Boolean(state.site && state.data[state.site] && radarBlind(state.data[state.site]));
/* A layer's legend rows at the site shown: where the radar cannot see the water, the water layer does not draw
   seasonal inundation (pipeline/tiles.py), so its legend leaves that row out. */
const legendOf = (l) => (l.radar === "seasonal" && siteBlind() ? l.legend.filter(([c]) => c !== SEASONAL) : l.legend);
/* The radar note of a layer at the site shown, or "". */
const radarNote = (l) => (l.radar && siteBlind() ? RADAR_NOTE[l.radar] : "");

function renderLayers() {
  const box = $("#layers");
  box.innerHTML = "";
  let group = null;
  for (const l of LAYERS) {
    if (l.group !== group) {
      group = l.group;
      box.insertAdjacentHTML("beforeend", `<div class="group-title">${group}</div>`);
    }
    const row = document.createElement("div");
    row.className = "layer" + (state.on[l.id] ? " open" : "");
    const note = radarNote(l);
    row.innerHTML = `
      <div class="layer-row">
        <input type="checkbox" id="cb-${l.id}" ${state.on[l.id] ? "checked" : ""}>
        <label for="cb-${l.id}"><span class="title">${l.title}</span><span class="sub">${l.sub}</span>${note ? `<span class="note">${note}</span>` : ""}</label>
        <button type="button" class="expand" aria-label="Show details">${chevron}</button>
      </div>
      <div class="layer-detail">
        <div class="opacity"><span>Opacity</span><input type="range" min="0.1" max="1" step="0.05" value="${state.opacity[l.id]}"></div>
        ${swatches(legendOf(l), l.ink)}
      </div>`;
    row.querySelector("input[type=checkbox]").addEventListener("change", (e) => {
      state.on[l.id] = e.target.checked;
      if (e.target.checked) row.classList.add("open");
      if (map.getLayer(`r-${l.id}`)) map.setLayoutProperty(`r-${l.id}`, "visibility", e.target.checked ? "visible" : "none");
    });
    row.querySelector(".expand").addEventListener("click", () => row.classList.toggle("open"));
    row.querySelector("input[type=range]").addEventListener("input", (e) => {
      state.opacity[l.id] = Number(e.target.value);
      if (map.getLayer(`r-${l.id}`)) map.setPaintProperty(`r-${l.id}`, l.vector ? "line-opacity" : "raster-opacity", state.opacity[l.id]);
    });
    box.appendChild(row);
  }
}

function renderBasemaps() {
  const box = $("#basemaps");
  box.innerHTML = BASEMAPS.map((b) => `<button type="button" data-id="${b.id}" aria-pressed="${b.id === basemapId()}">${b.label}</button>`).join("");
  inkSwatches();
  for (const btn of box.querySelectorAll("button")) btn.addEventListener("click", () => {
    const prev = styleUrl();
    state.basemap = btn.dataset.id;
    store.set("waterfields-basemap", state.basemap);
    renderBasemaps();
    if (styleUrl() !== prev) { setStyle(styleUrl()); return; }
    if (map.getLayer("r-truecolor")) map.setLayoutProperty("r-truecolor", "visibility", basemapId() === "imagery" ? "visible" : "none");
    applyLabels();
    applyInk();
  });
}

/* statistics panel: one section per topic, built as HTML text */
/* A label and value row; cls "stack" puts the value under the label, for long values and plain phrases.
   `note` is a short line in small type under the label, such as a threshold range. */
const kv = (label, value, cls = "", note = "") => `<div class="kv${cls ? ` ${cls}` : ""}">`
  + `<span>${label}${note ? `<small>${note}</small>` : ""}</span><span>${value}</span></div>`;

/* A paragraph of fine print inside a section. */
const fine = (text) => `<div class="fine tight">${text}</div>`;

/* Horizontal bars, each row [label, value, colour], scaled so that `total` fills the track. Values under 1 get
   two decimals, so that a small class does not read as 0.0. */
function bars(rows, total, unit) {
  return `<div class="bars">${rows.map(([label, v, c]) => `
    <div class="bar-row"><span>${label}</span><span class="val">${fmt(v, v > 0 && v < 1 ? 2 : 1)}${unit}</span>
    <div class="track"><div class="fill" style="width:${Math.min(100, (100 * v) / Math.max(total, 1e-9))}%;background:${c}"></div></div></div>`).join("")}</div>`;
}
const FIELD_COLOURS = { "Flood-recession fields": "#178c96", "Irrigated or multi-cropped": "#2e8540",
  "Rice-field landscape (rainfed)": "#b5c64b", "Open agricultural land": "#e8dccb" };
const WATER_COLOURS = { "River or stream": "#3b7fc4", "Straight canal (candidate)": "#0b7a75", "Enclosing moat (candidate)": ACCENT,
  "Rectilinear reservoir or moat segment (10 ha or more)": "#14a39c", "Rectilinear pond": "#5cc8c2", "Lake or irregular reservoir": "#1f4e89",
  "Irregular pond": "#6a3d9a", "Small water body (shape not resolved)": "#a75cc0" };

/* Radar cannot see a site's water when the median monthly VV over the optically mapped water stays above this
   (dB), the test of pipeline/maps.py (radar_cannot_see_water); radar water figures are then not shown. */
const RADAR_BLIND_DB = -15;

/* Least number of straight-sided bodies in the central window for the water rotation test to be judged
   (pipeline/water.py, rotation_control). */
const ROTATION_MIN_BODIES = 30;

/* Median of the monthly VV over optically mapped water (dB), or null when it was not measured. */
function vvOverWater(d) {
  const monthly = ((d.accuracy || {}).s1_vv_over_optical_water_db || {}).monthly_median || [];
  const v = monthly.filter((x) => x !== null).sort((a, b) => a - b);
  if (!v.length) return null;
  const mid = Math.floor(v.length / 2);
  return v.length % 2 ? v[mid] : (v[mid - 1] + v[mid]) / 2;
}
const radarBlind = (d) => { const vv = vvOverWater(d); return vv !== null && vv > RADAR_BLIND_DB; };

/* A threshold range as "range 131.3 to 140.8 km²", or nothing. */
const range = (r, unit = " km²") => (r ? `range ${fmt(r[0])} to ${fmt(r[1])}${unit}` : "");

/* A 95% interval as "95% interval 29.9 to 33.8%", or nothing. */
const ci = (r) => (r ? `95% interval ${fmt(r[0])} to ${fmt(r[1])}%` : "");

/* A share that may be missing (null), as "35.7%" or "Not measured". */
const pct = (v) => (v === null || v === undefined ? "Not measured" : `${fmt(v)}%`);

/* A signed area, "+213.1 km²". */
const signedKm2 = (v) => `${v > 0 ? "+" : ""}${fmt(v)} km²`;

function renderStats() {
  const d = state.data[state.site];
  if (!d) return;
  $("#statsSite").textContent = d.name;
  const w = d.water, f = d.fields, s = d.seasons;
  const rel = f.relation_to_perennial_water || {};
  const irr = rel["Irrigated or multi-cropped"];
  // classes whose area rounds to 0.0 km² are left out of both sets of bars
  const wa = Object.entries(w.area_km2_by_class).filter(([, v]) => v > 0);
  const fa = Object.entries(f.area_km2_by_class).filter(([, v]) => v > 0);
  const cyc = f.cropping_cycles_share_pct;
  $("#stats").innerHTML = `
    <div class="section"><h3>Study area</h3>
      ${kv("Area", `${fmt(w.study_area_km2, 0)} km²`)}
      ${kv("Dry season", `${s.dry[0]} to ${s.dry[1]}`)}
      ${kv("Wet season", `${s.wet[0]} to ${s.wet[1]}`)}
      ${kv("Analysis grid", `EPSG:${siteInfo(state.site).epsg}`)}
      ${fine(d.aoi_source)}
    </div>
    <div class="section"><h3>Water, dry season (Sentinel-2)</h3>
      ${kv("Open water", `${fmt(w.dry_season_water_km2)} km²`, "", range(w.sensitivity.dry_season_water_km2_range))}
      ${kv("Share of the study area", `${fmt(w.dry_season_water_share_pct)}%`)}
      ${kv("Water bodies", fmt(w.feature_count, 0), "", "0.1 ha or larger")}
      ${kv("Water bodies per 100 km²", fmt(w.features_per_100km2))}
      ${kv("Ponds per 100 km²", fmt(w.ponds_per_100km2), "", "ponds and small water bodies, all under 10 ha")}
      ${bars(wa.map(([k, v]) => [k, v, WATER_COLOURS[k] || "#999"]), Math.max(...wa.map(([, v]) => v), 1e-9), " km²")}
    </div>
    <div class="section"><h3>Fields</h3>
      ${kv("Agricultural land", `${fmt(f.agricultural_km2)} km²`, "", `${fmt(f.agricultural_share_pct)}% of the study area`)}
      ${bars(fa.map(([k, v]) => [k, v, FIELD_COLOURS[k]]), f.agricultural_km2, " km²")}
      ${kv("Cropped twice or more", `${fmt(cyc["2"] + cyc["3"])}%`, "", "of farmland with six or more clear months")}
      ${kv("Concentrated farmland", `${fmt(f.concentrated_ag_km2)} km²`, "", "80% or more farmland within 1 km")}
      ${irr ? kv("Irrigated within 1 km of water", `${fmt(irr.within_1km_pct)}%`, "", "dry-season water") : ""}
      ${kv("All land within 1 km of water", pct(rel.all_land ? rel.all_land.within_1km_pct : null), "", "for comparison")}
    </div>
    ${radarWaterSection(d)}
    ${vegetationSection(d.change)}
    ${orientationSection(d)}
    ${d.ftw ? ftwSection(d.ftw) : ""}
    ${d.accuracy ? agreementSection(d.accuracy, radarBlind(d)) : ""}
    <div class="fine">All figures come from maps made with open data. A validation sample has been designed but not yet labelled, so the areas are mapped areas rather than error-adjusted estimates. Ranges come from moving each threshold, and class names describe shapes rather than identifying or dating features.<br>
      <a href="data/${d.key}_water_features.csv" download>Water bodies table (CSV)</a> · <a href="data/${d.key}.json" download>Site statistics (JSON)</a> · <a href="data/waterfields_comparative_${d.year}.csv" download>Comparative dataset, four sites (CSV)</a></div>`;
}

/* Expansion and contraction of water on the one radar basis, or why they cannot be measured here. */
function radarWaterSection(d) {
  const r = d.change.radar_water_change_km2, ws = d.water.sensitivity;
  if (radarBlind(d)) {
    return `<div class="section"><h3>Seasonal water (Sentinel-1)</h3>
      ${kv("Radar water", "Not measurable at this site", "stack")}
      ${fine(`Radar backscatter (VV) over the water mapped from Sentinel-2 stays at about ${fmt(vvOverWater(d))} dB in every month, well above the open-water range, so the radar does not see this water. Expansion, contraction and seasonal inundation are therefore not reported.`)}
    </div>`;
  }
  return `<div class="section"><h3>Seasonal water (Sentinel-1)</h3>
    ${kv("Dry-season water", `${fmt(r.dry_season_water)} km²`)}
    ${kv("Wet-season maximum", `${fmt(r.wet_season_maximum_water)} km²`, "", range(ws.s1_wet_max_water_km2_range))}
    ${kv("Water in both seasons", `${fmt(r.water_in_both_seasons)} km²`)}
    ${kv("Expansion (seasonal inundation)", `${fmt(r.expansion_wet_season_only)} km²`, "", range(ws.seasonal_inundation_km2_range))}
    ${kv("Contraction (dry season only)", `${fmt(r.contraction_dry_season_only)} km²`)}
    ${kv("Net change, wet minus dry", signedKm2(r.net_change_wet_minus_dry))}
    ${fine("One VV threshold for both seasons: the dry-season median against any wet-season month.")}
  </div>`;
}

/* Vegetation change between the seasons, on farmland and over the whole study area. */
function vegetationSection(c) {
  return `<div class="section"><h3>Seasonal vegetation (Sentinel-2)</h3>
    ${kv("Farmland with a clear view", `${fmt(c.ag_ndvi_change_coverage_pct)}%`, "", "in both seasons")}
    ${kv("Farmland greener, +0.15 NDVI or more", `${fmt(c.ag_ndvi_greening_km2)} km²`, "", `${fmt(c.ag_ndvi_greening_share_pct)}% of all farmland`)}
    ${kv("Farmland browner, −0.15 NDVI or less", `${fmt(c.ag_ndvi_browning_km2)} km²`, "", `${fmt(c.ag_ndvi_browning_share_pct)}% of all farmland`)}
    ${kv("Whole area greener", `${fmt(c.ndvi_greening_km2)} km²`)}
    ${kv("Whole area browner", `${fmt(c.ndvi_browning_km2)} km²`)}
    ${kv("Whole area with a clear view", `${fmt(c.ndvi_change_coverage_pct)}%`)}
    ${kv("Farmland with radar VH rise over 3 dB", `${fmt(c.ag_vh_rise_over_3db_km2)} km²`, "", "cloud-free companion")}
  </div>`;
}

/* Share of rectilinear water and of field edges within 7.5 degrees of north-south or east-west, each with its
   30-degree rotation control, and the Fields of The World share for comparison. */
function orientationSection(d) {
  const w = d.water, fo = d.fields.field_orientation;
  const ftwShare = d.ftw && d.ftw.orientation ? d.ftw.orientation.cardinal_share_pct : null;
  return `<div class="section"><h3>Orientation, within 7.5° of north-south or east-west</h3>
    ${kv("Rectilinear water bodies", `${fmt(w.cardinal_share_of_rectilinear_pct)}%`)}
    ${waterRotationNote(w.orientation_control)}
    ${kv("Field-edge energy", `${fmt(fo.cardinal_share_pct)}%`, "", ci(fo.cardinal_share_ci95_pct))}
    ${kv("Edges outside farmland (control)", `${fmt(fo.controls.non_agricultural_edges_pct)}%`, "", ci(fo.controls.non_agricultural_ci95_pct))}
    ${fieldRotationNote(fo.controls.rotation_test)}
    ${ftwShare !== null ? kv("Fields of The World fields", pct(ftwShare), "", "fields of 2 ha or more") : ""}
    ${kv("Expected if random", `${fmt(fo.uniform_expectation_pct)}%`)}
  </div>`;
}

/* The 30° rotation control for water bodies, with the verdict recorded by pipeline/water.py. */
function waterRotationNote(o) {
  if (!o) return "";
  if (o.passed === null) {
    const n = o.window_rectilinear_bodies;
    return fine(n < ROTATION_MIN_BODIES
      ? `Rotation test for water: only ${fmt(n, 0)} straight-sided ${n === 1 ? "body lies" : "bodies lie"} in the central test window (${ROTATION_MIN_BODIES} are needed), so this share is not interpreted.`
      : `Rotation test for water: fewer than ${ROTATION_MIN_BODIES} straight-sided bodies remain in the rotated test window, so this share is not interpreted.`);
  }
  return fine(o.passed
    ? `Rotation test for water passed: after a 30° rotation, ${fmt(o.after_30deg_rotation_original_axes_pct)}% of rectilinear bodies follow the ground axes and ${fmt(o.after_30deg_rotation_image_axes_pct)}% the image axes.`
    : `Rotation test for water failed: after a 30° rotation, ${fmt(o.after_30deg_rotation_image_axes_pct)}% of rectilinear bodies follow the image axes and ${fmt(o.after_30deg_rotation_original_axes_pct)}% the ground axes, so this share is not interpreted.`);
}

/* The 30° rotation control for field edges: a share that stays on the image axes after rotation is pixel-grid
   bias (the rule of pipeline/dataset.py, passes_rotation_test). */
function fieldRotationNote(r) {
  if (!r) return "";
  const pass = r.after_30deg_rotation_image_axes_pct < r.window_share_pct;
  return fine(pass
    ? `Rotation test for field edges passed: after a 30° rotation, the share on the image axes falls from ${fmt(r.window_share_pct)}% to ${fmt(r.after_30deg_rotation_image_axes_pct)}%.`
    : `Rotation test for field edges failed: after a 30° rotation, the share on the image axes is ${fmt(r.after_30deg_rotation_image_axes_pct)}%, against ${fmt(r.window_share_pct)}% before, so this share is not interpreted.`);
}

function ftwSection(ftw) {
  return `<div class="section"><h3>Fields of The World 2024</h3>
    ${kv("Fields", fmt(ftw.fields, 0))}
    ${kv("Median field size", `${fmt(ftw.median_field_ha, 2)} ha`)}
    ${kv("Fields per km²", fmt(ftw.fields_per_km2))}
    ${kv("Mapped farmland covered by FTW", `${fmt(ftw.share_of_our_agriculture_covered_pct)}%`)}
  </div>`;
}

/* Agreement of each result with the reference datasets of the brief, from accuracy_stats.json. Radar rows are
   not shown as scores where the radar cannot see the water (blind). */
function agreementSection(a, blind) {
  function f1(x, radar) {
    if (radar && blind) return ["Not measurable", "radar does not see the water here"];
    if (x.f1 === null) return ["Too little to compare", ""];
    return [fmt(x.f1, 2), `P ${fmt(x.precision, 2)}, R ${fmt(x.recall, 2)}`];
  }
  // a reference missing from the file (not computed for this site) gets no row
  const row = (label, x, radar = false, extra = "") => {
    if (!x) return "";
    const [value, note] = f1(x, radar);
    return kv(label, value, "", radar && blind ? note : [note, extra].filter(Boolean).join("; "));
  };
  const dw = a.dry_season_water || {}, ww = a.wet_season_water || {}, ag = a.agricultural_domain || {};
  const si = a.seasonal_inundation || {};
  return `<div class="section"><h3>Agreement with reference maps (F1)</h3>
    ${fine("Agreement between maps, not accuracy against ground truth. P is the share of the mapped class that the reference also shows, and R is the share of the reference that is mapped.")}
    ${row("Dry-season water vs JRC", dw["JRC occurrence 50% or more (1984-2020, 30 m)"], false, "water in at least half of the 1984 to 2020 observations")}
    ${row("Dry-season water vs WorldCover", dw["ESA WorldCover 2021 water"], false, "2021 water class")}
    ${row("Wet-season water vs Sentinel-2", ww["Sentinel-2 wet-season MNDWI above 0 (clear pixels only)"], true, "clear pixels only")}
    ${row("Seasonal inundation vs JRC", si["JRC seasonal water (1 to 11 months in 2020)"], true, "seasonal water in 2020")}
    ${row("Permanent water vs JRC", a.permanent_water_vs_jrc_permanent, true, "10 or more months against 12 in 2020")}
    ${row("Farmland vs WorldCover cropland", ag["ESA WorldCover 2021 cropland"], false, "not independent")}
    ${row("Farmland vs Fields of The World", ag["Fields of The World 2024 fields"])}
  </div>`;
}

/* identify: read the clicked pixel from each visible layer's own tile, then name it by its legend colour */

/* Web Mercator tile x, y at zoom z holding a point, and the point's position inside that tile (0 to 1). */
function lonLatToTile(lon, lat, z) {
  const n = 2 ** z;
  const x = ((lon + 180) / 360) * n;
  const r = (lat * Math.PI) / 180;
  const y = ((1 - Math.log(Math.tan(r) + 1 / Math.cos(r)) / Math.PI) / 2) * n;
  return { x: Math.floor(x), y: Math.floor(y), fx: x % 1, fy: y % 1 };
}
/* RGB of one layer's tile pixel under a point, or null where the layer is transparent (not mapped). */
async function pixelAt(url, lon, lat) {
  const p = archive(url);
  const h = await p.getHeader();
  const z = Math.max(h.minZoom, Math.min(h.maxZoom, Math.round(map.getZoom())));
  const t = lonLatToTile(lon, lat, z);
  const res = await p.getZxy(z, t.x, t.y);
  if (!res || !res.data) return null;
  const bmp = await createImageBitmap(new Blob([res.data]));
  const canvas = new OffscreenCanvas(bmp.width, bmp.height);
  const ctx = canvas.getContext("2d", { willReadFrequently: true });
  ctx.drawImage(bmp, 0, 0);
  const d = ctx.getImageData(Math.floor(t.fx * bmp.width), Math.floor(t.fy * bmp.height), 1, 1).data;
  return d[3] === 0 ? null : [d[0], d[1], d[2]];
}
/* "#3b7fc4" -> [59, 127, 196] */
function hexToRgb(hex) {
  return [parseInt(hex.slice(1, 3), 16), parseInt(hex.slice(3, 5), 16), parseInt(hex.slice(5, 7), 16)];
}

/* The legend entry whose colour is closest to a pixel's colour. */
function nearestLegend(layer, rgb) {
  let best = null, bd = Infinity;
  for (const [c, t] of layer.legend) {
    const q = hexToRgb(c);
    const dd = (q[0] - rgb[0]) ** 2 + (q[1] - rgb[1]) ** 2 + (q[2] - rgb[2]) ** 2;
    if (dd < bd) { bd = dd; best = [c, t]; }
  }
  return best;
}
/* The site's water bodies table (one object per row, values as text), read once. */
async function waterFeatures(key) {
  if (!state.features[key]) {
    const text = await fetch(`data/${key}_water_features.csv`).then((r) => r.text());
    const [head, ...lines] = text.trim().split(/\r?\n/);  // the tables are written with CRLF line endings
    const cols = head.split(",").map((c) => c.trim());
    const rows = [];
    for (const line of lines) {
      const values = line.split(",");
      const row = {};
      cols.forEach((col, i) => { row[col] = (values[i] ?? "").trim(); });
      rows.push(row);
    }
    state.features[key] = rows;
  }
  return state.features[key];
}

/* Distance in metres between two [lon, lat] points (equirectangular, fine within a few km). */
function metres(a, b) {
  const R = 6371000, r = Math.PI / 180;
  const x = (b[0] - a[0]) * r * Math.cos(((a[1] + b[1]) / 2) * r), y = (b[1] - a[1]) * r;
  return Math.sqrt(x * x + y * y) * R;
}
/* Search radius for the nearest water body of any class, and the least reach allowed for a body of the clicked
   class (a body counts when the click lies within its length of its centroid). */
const NEAREST_ANY_M = 3000;
const SAME_CLASS_MIN_M = 300;

/* The water body to describe at a click: the nearest one of the clicked class whose length reaches the click,
   or else the nearest one of any class within NEAREST_ANY_M. Returns [row, distance in m, same class] or null. */
function nearestBody(feats, lnglat, legendColour) {
  const cls = Object.keys(WATER_COLOURS).find((k) => WATER_COLOURS[k] === legendColour);
  let same = null, sd = Infinity, any = null, ad = Infinity;
  for (const f of feats) {
    const dd = metres(lnglat, [Number(f.lon), Number(f.lat)]);
    if (dd < ad) { ad = dd; any = f; }
    if (f.class === cls && dd < sd && dd <= Math.max(Number(f.length_m) || 0, SAME_CLASS_MIN_M)) { sd = dd; same = f; }
  }
  if (same) return [same, sd, true];
  return any && ad < NEAREST_ANY_M ? [any, ad, false] : null;
}

/* True when a [lon, lat] point lies inside the study-area polygon (ray casting); the layers are clipped to it. */
function insideStudyArea(pt) {
  const d = state.data[state.site];
  if (!d) return false;
  const ring = d.boundary;
  let inside = false;
  for (let i = 0, j = ring.length - 1; i < ring.length; j = i++) {
    const [xi, yi] = ring[i], [xj, yj] = ring[j];
    if ((yi > pt[1]) !== (yj > pt[1]) && pt[0] < ((xj - xi) * (pt[1] - yi)) / (yj - yi) + xi) inside = !inside;
  }
  return inside;
}

/* "13.41234° N, 103.87654° E" */
const latLonText = (lat, lng) => `${Math.abs(lat).toFixed(5)}° ${lat < 0 ? "S" : "N"}, ${Math.abs(lng).toFixed(5)}° ${lng < 0 ? "W" : "E"}`;

/* Side and size of the popup at a point (MapLibre pixels from the map's top left). The popup opens on the side of
   the point with more room, above or below it, between the map buttons at the top right (basemap and zoom) and the
   part of the map that an open panel leaves visible (narrow screens), and `room` is the height it has there. Near
   the left or right edge it opens towards the middle. Below a point beside those buttons it opens to their left,
   or, where the map leaves no room there, a little lower than the point, clear of them. Its width is at most two
   thirds of the map. */
const POPUP_TIP = 10;
function popupPlace(pt) {
  const el = map.getContainer();
  const width = Math.min(300, Math.floor((el.clientWidth - 8) / 1.5));
  const visible = el.clientHeight - (NARROW.matches ? dockHeight() : 0);
  // left, right and bottom edge of each group of buttons (the map and the basemap buttons share a top left corner)
  const buttons = [$("#basemaps"), el.querySelector(".maplibregl-ctrl-top-right")]
    .map((b) => [b.offsetLeft, b.offsetLeft + b.offsetWidth, b.offsetTop + b.offsetHeight]);
  const span = { "-left": [pt.x, pt.x + width], "": [pt.x - width / 2, pt.x + width / 2], "-right": [pt.x - width, pt.x] };
  // lowest bottom edge of the buttons that a popup on this side would reach across, or 0
  const reach = (s) => Math.max(0, ...buttons.filter(([l, r]) => span[s][1] > l && span[s][0] < r).map(([, , b]) => b));
  let side = pt.x < width / 2 ? "-left" : pt.x > el.clientWidth - width / 2 ? "-right" : "";
  let dy = 0;
  if (pt.y + POPUP_TIP < reach(side)) {
    if (pt.x >= width && pt.y + POPUP_TIP >= reach("-right")) side = "-right";
    else dy = reach(side) + 2 - pt.y - POPUP_TIP;
  }
  const above = pt.y - reach(side) - 8, below = visible - pt.y - dy - 8;
  return {
    anchor: (below > above ? "top" : "bottom") + side,  // MapLibre names the part of the popup at the point
    offset: [0, dy],
    width,
    room: Math.max(above, below),
  };
}
/* Hold the popup's body to the room left under its tip, border and head as drawn (the head wraps on a narrow map):
   at least one line, which scrolls, so that the room holds on maps down to about 235 px tall (a phone held
   sideways), and at most 360 px. */
function fitPopupBody(popup, room) {
  const el = popup.getElement();
  const head = el.querySelector(".pop-head").offsetHeight;
  el.querySelector(".pop-body").style.maxHeight = `${Math.max(32, Math.min(360, room - POPUP_TIP - head - 6))}px`;
}
/* Least strip of map, in px, that the popup needs beside an open panel; on a shorter strip the map is shown alone. */
const POPUP_MIN_STRIP = 220;

/* Close the identify popup, if one is open; its close handler clears the ring. */
function closePopup() {
  if (state.popup) state.popup.remove();
}

/* Click on the map: a popup at that point with the class of every visible layer there and the nearest mapped
   water body. A click elsewhere replaces it; its close button or Escape closes it. */
async function identify(e) {
  const { lng, lat } = e.lngLat;
  const run = ++state.identifyRun;  // a later click wins over a slower earlier one
  closePopup();
  state.pick = [lng, lat];
  map.getSource("pick")?.setData(pickData());
  // a panel over a short map (a phone held sideways) leaves too thin a strip for the popup
  if (NARROW.matches && map.getContainer().clientHeight - dockHeight() < POPUP_MIN_STRIP) showView("map");
  const place = popupPlace(e.point);
  const head = `<div class="pop-head">${latLonText(lat, lng)}</div>`;
  const body = (html) => `${head}<div class="pop-body">${html}</div>`;
  const popup = new maplibregl.Popup({ className: "identify", anchor: place.anchor, offset: place.offset,
    maxWidth: `${place.width}px`, closeOnClick: false, focusAfterOpen: false })
    .setLngLat(e.lngLat)
    .setHTML(body('<div class="fine">Reading the layers at this point</div>'))
    .addTo(map);
  fitPopupBody(popup, place.room);
  popup.on("close", () => {
    if (state.popup !== popup) return;
    state.popup = null;
    state.pick = null;
    map.getSource("pick")?.setData(pickData());
  });
  state.popup = popup;

  const blocks = [];
  const shown = LAYERS.filter((l) => state.on[l.id] && !l.vector);
  if (!insideStudyArea([lng, lat])) {
    blocks.push(fine("Outside the study area. The layers are mapped inside it only."));
  } else if (!shown.length) {
    blocks.push(fine("Switch on an analysis layer to read its class here."));
  } else {
    let waterClass = null, waterColour = null;
    for (const l of shown) {
      let html;
      try {
        const rgb = await pixelAt(tileUrl(state.site, l.id), lng, lat);
        if (!rgb) html = '<div class="id-value"><span>Not mapped here</span></div>';
        else {
          const [c, t] = nearestLegend(l, rgb);
          if (l.id === "water") { waterClass = t; waterColour = c; }
          html = `<div class="id-value"><i style="background:${c}"></i><span>${t}</span></div>`;
          if (l.radar === "all" && siteBlind()) html += fine(RADAR_NOTE.all);
        }
      } catch (err) { html = '<div class="id-value"><span>Could not read this layer</span></div>'; }
      blocks.push(`<div class="id-block"><h4>${l.title}</h4>${html}</div>`);
    }
    if (waterClass && !waterClass.startsWith("Seasonal")) {
      const found = nearestBody(await waterFeatures(state.site), [lng, lat], waterColour);
      if (found) {
        const [best, bd, sameClass] = found;
        blocks.push(`<div class="id-block"><h4>Nearest mapped water body${sameClass ? " of this class" : ""}, centroid ${fmt(bd, 0)} m away</h4>
          ${kv("Class", best.class, "stack")}${kv("Area", `${fmt(best.area_ha)} ha`)}${kv("Length", `${fmt(best.length_m, 0)} m`)}
          ${kv("Mean width", `${fmt(best.mean_width_m, 0)} m`)}${kv("Rectangularity", fmt(best.rectangularity, 2))}
          ${kv("Azimuth of long axis", `${fmt(best.azimuth_deg)}°`)}${kv("Cardinally aligned", best.cardinal === "1" ? "Yes" : "No")}</div>`);
      }
    }
  }
  if (run !== state.identifyRun || state.popup !== popup) return;  // closed or replaced while the layers were read
  popup.setHTML(body(blocks.join("")));
  fitPopupBody(popup, place.room);
}

/* site, view, theme */
/* Switch to a site: read its statistics, draw its layers and zoom to its study area. */
async function loadSite(key) {
  closePopup();
  state.site = key;
  state.pick = null;
  $("#site").value = key;
  if (!state.data[key]) state.data[key] = await fetch(`data/${key}.json`).then((r) => r.json());
  const d = state.data[key];
  renderStats();
  renderLayers();  // the radar notes and the water legend depend on the site
  mountSite();
  $("#crs").textContent = `Display EPSG:3857 · analysis grid EPSG:${siteInfo(key).epsg}`;
  history.replaceState(null, "", `#${key}`);
  const lons = d.boundary.map((p) => p[0]), lats = d.boundary.map((p) => p[1]);
  const [w, s, e, n] = [Math.min(...lons), Math.min(...lats), Math.max(...lons), Math.max(...lats)];
  map.fitBounds([[w, s], [e, n]], { padding: fitPadding(), duration: 0 });
  // On a short map the fitted zoom can fall below the least zoom at which the tile archives are drawn, leaving the
  // study area empty, so the map then zooms in to it about the area's centre. Map zooms count 512 px tiles and
  // raster zooms are rounded, so the 256 px tiles of the archives' least zoom z are drawn from map zoom z - 1.5.
  const minZoom = await archive(tileUrl(key, "truecolor")).getHeader().then((h) => h.minZoom - 1.5, () => 0);
  if (state.site === key && map.getZoom() < minZoom) map.easeTo({ zoom: minZoom, around: [(w + e) / 2, (s + n) / 2], duration: 0 });
}

/* Padding for fitting a study area: on narrow screens the open panel lies over the bottom of the map, with the
   scale bar just above it, so the area is fitted into the part above them, below the basemap buttons. */
const FIT_MIN_PX = 96;  // least map height left for the study area on a short screen, such as a phone held sideways
function fitPadding() {
  if (!NARROW.matches) return 32;
  const above = 8 + 31 + 12;  // the scale bar sits 8 px above the panel (style.css) and is about 31 px tall
  const room = map.getContainer().clientHeight - 56 - FIT_MIN_PX;
  return { top: 56, bottom: Math.max(24, Math.min(dockHeight() + above, room)), left: 24, right: 24 };
}

function renderSiteOptions() {
  $("#site").innerHTML = SITES.map((s) => `<option value="${s.key}">${s.label}</option>`).join("");
}

/* Height of the panel that lies over the map on narrow screens (0 in the Map view); a hidden panel measures 0. */
const DOCKS = ["#layersDock", "#statsDock"].map((s) => $(s));
const dockHeight = () => Math.max(...DOCKS.map((d) => d.offsetHeight));

/* Narrow screens show one view at a time: the layer list, the map alone, or the statistics. Wide screens show
   both panels beside the map whatever the view. */
function showView(id) {
  document.body.dataset.view = id;
  for (const b of document.querySelectorAll(".views button")) b.setAttribute("aria-pressed", String(b.dataset.view === id));
}

map.on("style.load", () => {
  state.styleReady = true;
  rememberLabels();
  applyLabels();
  mountSite();
});
/* A full style reload: MapLibre 5 would otherwise diff the two styles, drop the site's layers and never
   fire 'style.load', which is where they are added back. */
function setStyle(url) {
  state.styleReady = false;
  map.setStyle(url, { diff: false });
}
map.on("click", identify);
map.on("mousemove", (e) => { $("#coords").textContent = `${e.lngLat.lat.toFixed(5)}, ${e.lngLat.lng.toFixed(5)}`; });
map.on("zoomend", () => { $("#zoom").textContent = `Zoom ${map.getZoom().toFixed(1)}`; });

renderSiteOptions();
$("#site").addEventListener("change", (e) => loadSite(e.target.value));
for (const b of document.querySelectorAll(".views button")) b.addEventListener("click", () => showView(b.dataset.view));
document.addEventListener("keydown", (e) => { if (e.key === "Escape") closePopup(); });
// the open panel's height, so that on narrow screens the scale bar can sit just above it
const dockObserver = new ResizeObserver(() => {
  document.body.style.setProperty("--dock-h", `${Math.round(dockHeight())}px`);
});
for (const d of DOCKS) dockObserver.observe(d);
// the credits end in an ellipsis on one line; the button or a click on them shows them all
function toggleCredits(open) {
  $("#status").classList.toggle("open", open);
  $("#creditsMore").setAttribute("aria-expanded", String(open));
  $("#creditsMore").setAttribute("aria-label", open ? "Show fewer credits" : "Show all credits");
}
$("#creditsMore").addEventListener("click", () => toggleCredits(!$("#status").classList.contains("open")));
$("#credits").addEventListener("click", (e) => { if (!e.target.closest("a")) toggleCredits(!$("#status").classList.contains("open")); });
$("#theme").value = savedTheme;
// the theme styles the page only; the basemap is chosen on its own
$("#theme").addEventListener("change", (e) => {
  const v = e.target.value;
  if (v === "auto") delete document.documentElement.dataset.theme;
  else document.documentElement.dataset.theme = v;
  store.set("waterfields-theme", v);
});

(async function init() {
  try {
    const ready = await fetch("data/sites.json").then((r) => (r.ok ? r.json() : null));
    if (Array.isArray(ready) && ready.length) {  // offer only the sites whose layers are published
      const keep = SITES.filter((x) => ready.includes(x.key));
      SITES.splice(0, SITES.length, ...keep);
      renderSiteOptions();
    }
  } catch (e) { /* no manifest: offer every site */ }
  renderLayers();
  renderBasemaps();
  if (NARROW.matches && SHORT.matches) showView("map");  // a phone held sideways opens on the map alone
  const start = SITES.some((s) => `#${s.key}` === location.hash) ? location.hash.slice(1) : SITES[0].key;
  await loadSite(start);
  $("#zoom").textContent = `Zoom ${map.getZoom().toFixed(1)}`;
})();
