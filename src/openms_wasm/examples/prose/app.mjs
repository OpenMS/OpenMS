// ProSE in the browser: a web application around the OpenMS search engine ProSE, compiled to
// WebAssembly. OpenMS runs in a Web Worker (worker.mjs); this page holds the settings form,
// generated from ProSE's parameters, and the views of the results.
//
// Served from the package directory: the module is ../../openms.mjs (or ?module=<url>, e.g.
// ../../index.mjs, the loader of the assembled package). More than one thread needs a
// cross-origin isolated page (COOP/COEP headers).

import { bars, histogram, ionType, mirrorSpectrum, scatter } from "./plots.mjs";

const $ = (id) => document.getElementById(id);
const status = (text) => { $("status").textContent = text; };
const params = new URLSearchParams(location.search);
const moduleUrl = new URL(params.get("module") ?? "../../openms.mjs", location.href).href;

/// The parameters shown first; all others are under "All parameters".
const COMMON = [
  "decoys", "enzyme", "peptide:enzyme_specificity", "peptide:missed_cleavages", "peptide:min_size", "peptide:max_size",
  "precursor:mass_tolerance_lower", "precursor:mass_tolerance_upper", "precursor:mass_tolerance_unit",
  "precursor:min_charge", "precursor:max_charge", "precursor:isotopes",
  "fragment:mass_tolerance", "fragment:mass_tolerance_unit",
  "modifications:fixed", "modifications:variable", "modifications:variable_max_per_peptide",
  "FDR:PSM", "FDR:protein", "report:top_hits",
];

//------------------------------------------------------------------------------
// The worker

let worker = null;
let nextRequest = 0;
const pending = new Map();
let onProgress = () => {};

function startWorker() {
  worker?.terminate();
  worker = new Worker(new URL("worker.mjs", import.meta.url), { type: "module" });
  worker.onmessage = ({ data }) => {
    if (data.type === "progress") {
      onProgress(data);
      return;
    }
    const request = pending.get(data.id);
    pending.delete(data.id);
    if (data.error !== undefined) request?.reject(new Error(data.error));
    else request?.resolve(data.result);
  };
  // An error the worker does not catch: OpenMS stopped, e.g. a WebAssembly trap on one of its threads.
  worker.onerror = (e) => {
    e.preventDefault();
    const error = Object.assign(new Error(e.message || "the worker stopped"), { trapped: true });
    for (const request of pending.values()) request.reject(error);
    pending.clear();
  };
}

function call(type, args = {}) {
  const id = nextRequest++;
  return new Promise((resolve, reject) => {
    pending.set(id, { resolve, reject });
    worker.postMessage({ id, type, ...args });
  });
}

//------------------------------------------------------------------------------
// Threads

const cores = navigator.hardwareConcurrency || 1;
const isolated = globalThis.crossOriginIsolated === true;
/// Default: all cores, at most 16. Measured with ProSE (from 4 to 32 threads), more threads hardly
/// shortened a search: decoy generation, parts of the index build and reading the spectra are serial.
const defaultThreads = Math.min(cores, 16);
for (let n = 1; n <= (isolated ? cores : 1); ++n) $("threads").append(new Option(String(n), String(n), false, n === defaultThreads));
$("threads-hint").textContent = isolated
  ? `${cores} logical cores reported by the browser.`
  : "One thread: the page is not cross-origin isolated (COOP/COEP headers), so threads cannot share memory.";

async function loadOpenMS() {
  status("Loading OpenMS…");
  $("run").disabled = true;
  startWorker();
  const numThreads = Number($("threads").value);
  const info = await call("init", { moduleUrl, numThreads });
  $("version").textContent = info.version;
  status(`OpenMS ${info.version} is ready: ${info.threading.threaded ? "threaded build" : "single-threaded build"}, ` +
         `${info.threading.numThreads} thread${info.threading.numThreads > 1 ? "s" : ""}.`);
  return info;
}

//------------------------------------------------------------------------------
// The settings form, generated from the parameter entries

let entries = [];

/// A form field for a parameter; @p fullName: captioned with its section (e.g. "fragment mass tolerance").
function inputFor(entry, fullName = false) {
  const label = document.createElement("label");
  label.title = entry.description;
  label.dataset.name = entry.name;
  if (entry.tags.includes("advanced")) label.classList.add("advanced");
  const caption = document.createTextNode((fullName ? entry.name : entry.name.split(":").pop()).replaceAll(/[:_]/g, " "));
  let input;
  if (entry.type === "bool") {
    label.classList.add("checkbox");
    input = Object.assign(document.createElement("input"), { type: "checkbox", checked: entry.value === true });
    label.append(input, caption);
  } else {
    label.append(caption);
    if (entry.type === "string" && entry.validStrings) {
      input = document.createElement("select");
      for (const s of entry.validStrings) input.append(new Option(s, s, false, s === entry.value));
    } else if (entry.type.endsWith("[]")) {
      input = document.createElement("textarea");
      input.rows = Math.min(4, Math.max(1, entry.value.length));
      input.value = entry.value.join("\n");
      input.placeholder = entry.type === "string[]" ? "one per line" : "values, one per line";
    } else if (entry.type === "int" || entry.type === "double") {
      input = Object.assign(document.createElement("input"), { type: "number", value: entry.value, step: entry.type === "int" ? 1 : "any" });
      if (entry.min !== undefined) input.min = entry.min;
      if (entry.max !== undefined) input.max = entry.max;
    } else {
      input = Object.assign(document.createElement("input"), { type: "text", value: entry.value ?? "" });
    }
    label.append(input);
  }
  input.dataset.name = entry.name;
  input.dataset.type = entry.type;
  return label;
}

function buildForm(parameterEntries) {
  entries = parameterEntries;
  const byName = new Map(entries.map((e) => [e.name, e]));
  $("common-parameters").replaceChildren(...COMMON.filter((n) => byName.has(n)).map((n) => inputFor(byName.get(n), true)));
  const sections = new Map();
  for (const e of entries) {
    if (COMMON.includes(e.name)) continue;
    const section = e.name.includes(":") ? e.name.slice(0, e.name.lastIndexOf(":")) : "general";
    if (!sections.has(section)) sections.set(section, []);
    sections.get(section).push(e);
  }
  const blocks = [];
  for (const [section, list] of sections) {
    const h = document.createElement("h3");
    h.textContent = section;
    const grid = document.createElement("div");
    grid.className = "parameters";
    grid.append(...list.map((e) => inputFor(e)));
    if (list.every((e) => e.tags.includes("advanced"))) h.classList.add("advanced-section");
    blocks.push(h, grid);
  }
  $("all-parameters").replaceChildren(...blocks);
}

$("show-advanced").addEventListener("change", (e) => document.body.classList.toggle("show-advanced", e.target.checked));

/// The values of the form as {name: value}.
function formValues() {
  const values = {};
  for (const input of document.querySelectorAll("#settings [data-type]")) {
    const { name, type } = input.dataset;
    let v;
    if (type === "bool") v = input.checked;
    else if (type === "int" || type === "double") v = Number(input.value);
    else if (type.endsWith("[]")) {
      const items = input.value.split("\n").map((s) => s.trim()).filter(Boolean);
      v = type === "string[]" ? items : items.map(Number);
    } else v = input.value;
    values[name] = v;
  }
  return values;
}

//------------------------------------------------------------------------------
// Files

const files = { spectra: null, database: null };

function setFile(kind, file) {
  files[kind] = file;
  $(`${kind}-name`).textContent = `${file.name} (${(file.size / 1e6).toFixed(1)} MB)`;
  $(`${kind}-drop`).classList.add("loaded");
  $("run").disabled = !(files.spectra && files.database);
}

for (const kind of ["spectra", "database"]) {
  const drop = $(`${kind}-drop`);
  $(`${kind}-file`).addEventListener("change", (e) => e.target.files[0] && setFile(kind, e.target.files[0]));
  drop.addEventListener("dragover", (e) => { e.preventDefault(); drop.classList.add("over"); });
  drop.addEventListener("dragleave", () => drop.classList.remove("over"));
  drop.addEventListener("drop", (e) => {
    e.preventDefault();
    drop.classList.remove("over");
    if (e.dataTransfer.files[0]) setFile(kind, e.dataTransfer.files[0]);
  });
}

//------------------------------------------------------------------------------
// Search

let results = null;

function showProgress(p) {
  if (p.event === "start" && p.label) $("progress-label").textContent = p.label.replace(/\.\.\.$/, "…");
  const fraction = p.end > p.begin ? (p.value - p.begin) / (p.end - p.begin) : 0;
  if (p.depth <= 1) $("progress-bar").style.width = `${Math.round(Math.min(1, Math.max(0, fraction)) * 100)}%`;
}

$("threads").addEventListener("change", async () => {
  results = null;
  currentPSM = null;
  $("results").hidden = true;
  $("download-idxml").hidden = $("download-tsv").hidden = true;
  try {
    await loadOpenMS();
    $("run").disabled = !(files.spectra && files.database);
  } catch (e) {
    status("OpenMS could not be loaded.");
    $("error").textContent = e.message;
    $("error").hidden = false;
  }
});

$("settings").addEventListener("submit", async (event) => {
  event.preventDefault();
  $("error").hidden = true;
  $("run").disabled = true;
  $("progress").hidden = false;
  $("progress-bar").style.width = "0%";
  const started = performance.now();
  const timer = setInterval(() => status(`Searching… ${((performance.now() - started) / 1000).toFixed(0)} s`), 500);
  onProgress = showProgress;
  try {
    results = await call("search", { spectraFile: files.spectra, databaseFile: files.database, parameters: formValues() });
    results.hasQValues = results.psms.some((p) => p.qValue !== null);
    clearInterval(timer);
    const s = results.summary;
    status(`${s.status === "EXECUTION_OK" ? "Done" : s.status}: ${results.psms.length} PSMs in ` +
           `${((performance.now() - started) / 1000).toFixed(1)} s.`);
    showResults();
  } catch (e) {
    clearInterval(timer);
    status("The search failed.");
    // A WebAssembly trap (e.g. more than the 4 GB the module can address) ends the module: start a new one.
    // A trap ends the module; after an OutOfMemory error it lives on, but with its memory grown to the limit
    // (WebAssembly memory does not shrink): start a new one in both cases.
    const trapped = e.trapped || /unreachable|out of memory|OutOfMemory|Aborted|RuntimeError/i.test(e.message);
    $("error").textContent = trapped
      ? `OpenMS failed (${e.message}), most likely because the search needs more than the 4 GB of memory a ` +
        "WebAssembly module can use. Set database:chunk_size (e.g. 2000) under All parameters to build and search the " +
        "index in parts, or use fewer variable modifications or a narrower enzyme specificity."
      : e.message;
    $("error").hidden = false;
    if (trapped) {
      results = null;
      $("results").hidden = true;
      $("download-idxml").hidden = $("download-tsv").hidden = true;
      // A new module in this page; the browser may not have released the memory of the old one yet:
      // then reloading the page is the way out.
      const restarted = await Promise.race([loadOpenMS().then(() => true, () => false),
                                            new Promise((resolve) => setTimeout(() => resolve(false), 30000))]);
      if (restarted) {
        status("The search failed: OpenMS ran out of memory and was restarted.");
      } else {
        worker?.terminate();
        status("The search failed: OpenMS ran out of memory and could not be restarted. Reload the page.");
        $("reload").hidden = false;
      }
    }
  } finally {
    $("progress").hidden = true;
    $("run").disabled = false;
  }
});

//------------------------------------------------------------------------------
// Results

const qThreshold = () => Number(params.get("q") ?? 0.01);
/// Whether the PSMs have q-values (the search had decoys); computed once per search.
const hasQValues = () => results.hasQValues;
/// Target PSMs that pass the q-value threshold (all targets without q-values).
const acceptedPSMs = () => results.psms.filter((p) => !p.decoy && (!hasQValues() || p.qValue <= qThreshold()));

function card(name, value) {
  const div = document.createElement("div");
  div.className = "card";
  div.innerHTML = '<div class="value"></div><div class="name"></div>';
  div.firstChild.textContent = value;
  div.lastChild.textContent = name;
  return div;
}

function showResults() {
  const { summary, psms, threading, memory } = results;
  const accepted = acceptedPSMs();
  const fdr = hasQValues() ? ` at ${qThreshold() * 100}% FDR` : " (no FDR: no decoys)";
  const peptides = new Set(accepted.map((p) => p.sequence));
  const proteins = new Set(accepted.flatMap((p) => p.proteins.split(";").filter(Boolean)));
  const seconds = (s) => `${s.toFixed(1)} s`;
  $("cards").replaceChildren(
    card("MS2 spectra", summary.run.ms2Spectra.toLocaleString()),
    card(`PSMs${fdr}`, accepted.length.toLocaleString()),
    card("peptides", peptides.size.toLocaleString()),
    card("proteins", proteins.size.toLocaleString()),
    card("indexed peptides", summary.search.indexedPeptides.toLocaleString()),
    card("index build", seconds(summary.search.secondsIndexBuild)),
    card("search total", seconds(summary.search.secondsTotal)),
    card("memory (peak)", `${(memory.heapSize / 2 ** 20).toFixed(0)} MB`),
    card("threads", String(threading.numThreads)),
  );

  const score = (p) => p.hyperscore;
  histogram($("score-plot"), [
    { name: "targets", cls: "bar-target", values: psms.filter((p) => !p.decoy).map(score) },
    { name: "decoys", cls: "bar-decoy", values: psms.filter((p) => p.decoy).map(score) },
  ].filter((s) => s.values.length > 0), { xLabel: "ln(hyperscore)" });
  histogram($("precursor-error-plot"), [{ name: "PSMs", cls: "bar-target", values: accepted.map((p) => p.precursorErrorPPM) }],
    { xLabel: "precursor error (ppm)" });

  const analysis = summary.modificationAnalysis;
  $("open-search").hidden = !analysis;
  if (analysis) {
    bars($("delta-mass-plot"), analysis.deltaMasses.entries.map((e) => ({
      x: e.deltaMass, height: e.count, cls: e.known ? "known" : "", label: e.modification || e.deltaMass.toFixed(2),
    })), { xLabel: "Δ mass (Da)", yLabel: "PSMs" });
    fillRows($("ptm-table"), analysis.ptms.entries.map((e) => [
      e.name, e.theoreticalMass.toFixed(4), e.observedMass.toFixed(4), e.count, e.uniquePeptides, e.percentage.toFixed(1), e.targetResidues,
    ]));
  }

  showPSMTable();
  showProteinTable();
  $("results").hidden = false;
  $("psm-detail").hidden = true;
  $("download-idxml").hidden = $("download-tsv").hidden = false;
}

function showProteinTable() {
  fillRows($("protein-table"), sortRows($("protein-table"), results.proteins.filter((p) => !p.decoy)).map((p) => [
    p.accession, fmt(p.score, 4), fmt(p.coverage, 1), p.description,
  ]));
}

const fmt = (v, digits) => (v === null || v === undefined || !Number.isFinite(v) ? "" : v.toFixed(digits));

function fillRows(table, rows, { onClick, classes } = {}) {
  table.tBodies[0].replaceChildren(...rows.map((cells, i) => {
    const tr = document.createElement("tr");
    cells.forEach((c, j) => {
      const td = document.createElement("td");
      td.textContent = c ?? "";
      if (classes?.[j]) td.className = classes[j];
      tr.append(td);
    });
    if (onClick) tr.addEventListener("click", () => onClick(i, tr));
    return tr;
  }));
}

/// Rows sorted by the column the user chose (a click on a header).
function sortRows(table, rows) {
  const th = table.querySelector("th[aria-sort]");
  if (!th) return rows;
  const { key, type } = th.dataset;
  const sign = th.getAttribute("aria-sort") === "ascending" ? 1 : -1;
  const value = (r) => r[key] ?? (type === "number" ? -Infinity : "");
  return [...rows].sort((a, b) => sign * (type === "number" ? value(a) - value(b) : String(value(a)).localeCompare(String(value(b)))));
}

for (const table of document.querySelectorAll("table.sortable")) {
  for (const th of table.querySelectorAll("th")) {
    th.addEventListener("click", () => {
      const ascending = th.getAttribute("aria-sort") !== "ascending";
      table.querySelectorAll("th").forEach((other) => other.removeAttribute("aria-sort"));
      th.setAttribute("aria-sort", ascending ? "ascending" : "descending");
      if (results) (table.id === "psm-table" ? showPSMTable : showProteinTable)();
    });
  }
}

const MAX_ROWS = 2000;
let shownPSMs = [];

function showPSMTable() {
  const filter = $("psm-filter").value.trim().toUpperCase();
  const pool = $("show-decoys").checked ? results.psms.filter((p) => !hasQValues() || p.qValue <= qThreshold()) : acceptedPSMs();
  const filtered = pool.filter((p) => !filter || p.sequence.toUpperCase().includes(filter) ||
    p.proteins.toUpperCase().includes(filter) || p.spectrum.toUpperCase().includes(filter));
  shownPSMs = sortRows($("psm-table"), filtered).slice(0, MAX_ROWS);
  $("psm-count").textContent = `${filtered.length.toLocaleString()} PSMs` + (filtered.length > MAX_ROWS ? `, the first ${MAX_ROWS} shown` : "");
  fillRows($("psm-table"), shownPSMs.map((p) => [
    p.index, p.sequence, p.charge, fmt(p.mz, 4), fmt(p.rt, 1), fmt(p.hyperscore, 3), p.qValue === null ? "" : p.qValue.toExponential(2),
    fmt(p.precursorErrorPPM, 2), fmt(p.deltaMass, 3), p.proteins,
  ]), {
    classes: [null, "sequence"],
    onClick: (i, tr) => {
      document.querySelectorAll("#psm-table tr.selected").forEach((r) => r.classList.remove("selected"));
      tr.classList.add("selected");
      showPSM(shownPSMs[i]);
    },
  });
  for (const [i, tr] of [...$("psm-table").tBodies[0].rows].entries()) {
    if (shownPSMs[i].decoy) tr.cells[1].classList.add("decoy");
  }
}

$("psm-filter").addEventListener("input", () => results && showPSMTable());
$("show-decoys").addEventListener("change", () => results && showPSMTable());

//------------------------------------------------------------------------------
// One PSM: sequence with fragment coverage, annotated spectrum, fragment errors

let currentPSM = null;

const selectedIons = () => [...document.querySelectorAll(".ion-type:checked")].map((c) => c.value);

async function showPSM(psm) {
  currentPSM = psm;
  const detail = await call("psm", { index: psm.index, ions: selectedIons(), losses: $("losses").checked });
  if (currentPSM !== psm) return; // another one was chosen meanwhile
  $("psm-title").textContent = `${detail.sequence}, ${detail.charge}+ · spectrum ${psm.spectrum} · RT ${fmt(psm.rt, 1)} s · ` +
    `score ${fmt(psm.hyperscore, 3)}` + (psm.qValue !== null ? ` · q ${psm.qValue.toExponential(2)}` : "");
  showSequence(detail);
  mirrorSpectrum($("spectrum-plot"), detail.peaks ?? { mz: [], intensity: [] }, detail.annotation);
  const unit = detail.annotation?.unit ?? "";
  scatter($("fragment-error-plot"), (detail.annotation?.matches ?? []).map((m) => ({
    x: m.mz, y: m.error, cls: ionType(m.ion), title: `${m.ion}: ${m.error.toFixed(2)} ${unit}`,
  })), { xLabel: "m/z", yLabel: `error (${unit})` });
  const meta = Object.entries(detail.meta).filter(([k]) => k !== "fragment_annotation").sort(([a], [b]) => a.localeCompare(b));
  $("psm-meta").replaceChildren(...meta.map(([k, v]) => {
    const tr = document.createElement("tr");
    tr.innerHTML = "<td></td><td></td>";
    tr.cells[0].textContent = k;
    tr.cells[1].textContent = typeof v === "number" ? +v.toPrecision(6) : Array.isArray(v) ? v.join(", ") : String(v);
    return tr;
  }));
  $("psm-detail").hidden = false;
}

for (const c of document.querySelectorAll(".ion-type, #losses")) {
  c.addEventListener("change", () => currentPSM && showPSM(currentPSM));
}

/// The residues, with a bracket for each cleavage site a prefix (a, b, c) or suffix (x, y, z) ion covers.
function showSequence(detail) {
  const n = detail.residues.length;
  const prefix = new Set(), suffix = new Set();
  for (const m of detail.annotation?.matches ?? []) {
    const match = /^([abcxyz])[.']*(\d+)/.exec(m.ion);
    if (!match) continue;
    const k = Number(match[2]);
    if ("abc".includes(match[1])) prefix.add(k); // cleavage after residue k
    else suffix.add(n - k);
  }
  const view = $("sequence-view");
  view.replaceChildren(...detail.residues.map((r, i) => {
    const div = document.createElement("div");
    div.className = "residue";
    div.textContent = r.code;
    if (r.modification) {
      div.classList.add("modified");
      div.title = r.modification;
    }
    if (prefix.has(i + 1) && i + 1 < n) div.classList.add("b");
    if (suffix.has(i + 1) && i + 1 < n) div.classList.add("y");
    const index = document.createElement("span");
    index.className = "index";
    index.textContent = i + 1;
    div.append(index);
    return div;
  }));
  const coverage = document.createElement("div");
  coverage.className = "coverage";
  const sites = new Set([...prefix, ...suffix].filter((k) => k >= 1 && k < n));
  coverage.textContent = `${sites.size} of ${n - 1} cleavage sites covered · prefix ions at ${prefix.size}, suffix ions at ${suffix.size}` +
    (detail.nTermModification ? ` · N-term: ${detail.nTermModification}` : "") +
    (detail.cTermModification ? ` · C-term: ${detail.cTermModification}` : "");
  view.append(coverage);
}

//------------------------------------------------------------------------------
// Downloads

function download(blob, name) {
  const url = URL.createObjectURL(blob);
  Object.assign(document.createElement("a"), { href: url, download: name }).click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
}

const stem = () => (files.spectra?.name ?? "results").replace(/\.mzML$/i, "");

$("reload").addEventListener("click", () => location.reload());

$("download-idxml").addEventListener("click", async () => {
  download(new Blob([await call("idXML")], { type: "application/xml" }), `${stem()}.idXML`);
});

$("download-tsv").addEventListener("click", () => {
  const columns = ["index", "spectrum", "rt", "mz", "charge", "sequence", "hyperscore", "qValue", "decoy", "precursorErrorPPM", "deltaMass", "proteins"];
  const lines = [columns.join("\t"), ...results.psms.map((p) => columns.map((c) => p[c] ?? "").join("\t"))];
  download(new Blob([lines.join("\n") + "\n"], { type: "text/tab-separated-values" }), `${stem()}_psms.tsv`);
});

//------------------------------------------------------------------------------

try {
  const info = await loadOpenMS();
  buildForm(info.parameters);
} catch (e) {
  status("OpenMS could not be loaded.");
  $("error").textContent = e.message;
  $("error").hidden = false;
}
