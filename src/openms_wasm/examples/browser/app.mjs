// Demo web application of the OpenMS WebAssembly module (see ../../README.md).
// Served from the package directory: the module is ../../openms.mjs.

import createOpenMS from "../../openms.mjs";

const $ = (id) => document.getElementById(id);
const status = (text) => { $("status").textContent = text; };
/// Lets the browser repaint (e.g. a status message) before a long synchronous computation.
const repaint = () => new Promise((resolve) => setTimeout(resolve, 20));

status("Loading OpenMS…");
const OpenMS = await createOpenMS({ printErr: (line) => console.warn(line) });
OpenMS.setLogEnabled("info", false);
$("version").textContent = OpenMS.getVersion();
status(`OpenMS ${OpenMS.getVersion()} is ready.`);

//------------------------------------------------------------------------------
// Tabs

for (const tab of document.querySelectorAll(".tabs button")) {
  tab.addEventListener("click", () => {
    for (const other of document.querySelectorAll(".tabs button")) {
      const selected = other === tab;
      other.setAttribute("aria-selected", String(selected));
      $(other.dataset.tab).hidden = !selected;
    }
  });
}

//------------------------------------------------------------------------------
// Plots (SVG)

const SVG = "http://www.w3.org/2000/svg";

function el(name, attributes = {}, text) {
  const e = document.createElementNS(SVG, name);
  for (const [k, v] of Object.entries(attributes)) e.setAttribute(k, v);
  if (text !== undefined) e.textContent = text;
  return e;
}

function niceTicks(min, max, count = 5) {
  const span = max - min || 1;
  const step = 10 ** Math.floor(Math.log10(span / count));
  const factor = [1, 2, 5, 10].find((f) => span / (step * f) <= count) ?? 10;
  const s = step * factor;
  const ticks = [];
  for (let t = Math.ceil(min / s) * s; t <= max; t += s) ticks.push(t);
  return ticks;
}

/// Axes of a plot in @p svg; returns the mapping of data to pixel coordinates.
function axes(svg, xMin, xMax, yMax, xLabel) {
  const [, , width, height] = svg.getAttribute("viewBox").split(" ").map(Number);
  const m = { left: 48, right: 12, top: 12, bottom: 32 };
  const x = (v) => m.left + ((v - xMin) / (xMax - xMin || 1)) * (width - m.left - m.right);
  const y = (v) => height - m.bottom - (v / (yMax || 1)) * (height - m.top - m.bottom);
  svg.replaceChildren();
  svg.append(el("line", { class: "axis", x1: m.left, x2: width - m.right, y1: y(0), y2: y(0) }));
  for (const t of niceTicks(xMin, xMax, Math.max(2, Math.floor(width / 110)))) {
    svg.append(el("text", { x: x(t), y: height - m.bottom + 16, "text-anchor": "middle" }, +t.toFixed(4)));
  }
  svg.append(el("text", { x: width - m.right, y: height - 4, "text-anchor": "end" }, xLabel));
  svg.append(el("text", { x: m.left - 6, y: y(yMax) + 4, "text-anchor": "end" }, yMax.toExponential(1)));
  return { x, y };
}

function stickPlot(svg, mz, intensity, { labels = [], classes = [] } = {}) {
  if (mz.length === 0) { svg.replaceChildren(); return; }
  const pad = (mz[mz.length - 1] - mz[0]) * 0.05 || 1;
  const yMax = Math.max(...intensity) * (labels.length ? 1.15 : 1.02);
  const { x, y } = axes(svg, mz[0] - pad, mz[mz.length - 1] + pad, yMax, "m/z");
  let lastLabelEnd = -Infinity; // labels are left out where they would overlap the previous one
  for (let i = 0; i < mz.length; ++i) {
    svg.append(el("line", { class: `stick ${classes[i] ?? ""}`, x1: x(mz[i]), x2: x(mz[i]), y1: y(0), y2: y(intensity[i]) }));
    const halfWidth = (labels[i]?.length ?? 0) * 3.5;
    if (labels[i] && x(mz[i]) - halfWidth > lastLabelEnd + 2) {
      svg.append(el("text", { class: `label-${classes[i]}`, x: x(mz[i]), y: y(intensity[i]) - 4, "text-anchor": "middle" }, labels[i]));
      lastLabelEnd = x(mz[i]) + halfWidth;
    }
  }
}

function linePlot(svg, xs, ys, xLabel) {
  if (xs.length === 0) { svg.replaceChildren(); return; }
  const { x, y } = axes(svg, xs[0], xs[xs.length - 1], Math.max(...ys) * 1.02, xLabel);
  const points = Array.from(xs, (v, i) => `${x(v).toFixed(1)},${y(ys[i]).toFixed(1)}`).join(" ");
  svg.append(el("polyline", { class: "line", points }));
}

//------------------------------------------------------------------------------
// Peptide

function showPeptide(event) {
  event?.preventDefault();
  const sequence = $("peptide-sequence").value.trim();
  const charge = Number($("peptide-charge").value) || 1;
  $("peptide-error").hidden = true;
  let peptide;
  try {
    peptide = OpenMS.AASequence.fromString(sequence);
  } catch (e) {
    $("peptide-error").textContent = e.message;
    $("peptide-error").hidden = false;
    return;
  }
  const formula = peptide.getFormula();
  const rows = [
    ["Sequence", peptide.toString()],
    ["Length", peptide.size()],
    ["Formula", formula.toString()],
    ["Monoisotopic mass", peptide.getMonoWeight().toFixed(5)],
    ["Average mass", peptide.getAverageWeight().toFixed(4)],
    ...[1, 2, 3, 4].map((z) => [`m/z ${z}+`, peptide.getMZ(z).toFixed(5)]),
  ];
  $("peptide-masses").replaceChildren(...rows.map(([k, v]) => {
    const tr = document.createElement("tr");
    tr.innerHTML = "<td></td><td></td>";
    tr.cells[0].textContent = k;
    tr.cells[1].textContent = v;
    return tr;
  }));

  // isotope distribution of the ion
  const generator = new OpenMS.CoarseIsotopePatternGenerator(6);
  const ionFormula = peptide.getFormula(OpenMS.ResidueType.Full, charge);
  const isotopes = ionFormula.getIsotopeDistribution(generator);
  const { mz, intensity } = isotopes.getPeaks();
  stickPlot($("peptide-isotopes"), mz.map((m) => m / charge), intensity);
  [generator, ionFormula, isotopes].forEach((o) => o.delete());

  // b and y ions, annotated
  const tsg = new OpenMS.TheoreticalSpectrumGenerator();
  tsg.setParameterValues({ add_metainfo: true, add_first_prefix_ion: true });
  const spectrum = new OpenMS.MSSpectrum();
  tsg.getSpectrum(spectrum, peptide, 1, Math.min(charge, 2));
  const peaks = spectrum.getPeaks();
  const names = spectrum.getStringDataArrays().find((a) => a.name === "IonNames")?.data ?? [];
  stickPlot($("peptide-fragments"), peaks.mz, peaks.intensity, {
    labels: names.map((n) => n.replace(/\+$/, "")),
    classes: names.map((n) => (n.startsWith("y") ? "y" : "b")),
  });
  [tsg, spectrum, formula, peptide].forEach((o) => o.delete());
}

$("peptide-form").addEventListener("submit", showPeptide);
showPeptide();

//------------------------------------------------------------------------------
// Digestion

for (const name of OpenMS.ProteaseDB.getAllNames()) {
  $("digest-enzyme").append(new Option(name, name, false, name === "Trypsin"));
}

function digest(event) {
  event?.preventDefault();
  const protein = $("digest-protein").value.replace(/\s+/g, "").toUpperCase();
  const digestion = new OpenMS.ProteaseDigestion();
  try {
    digestion.setEnzyme($("digest-enzyme").value);
    digestion.setMissedCleavages(Number($("digest-missed").value));
    const peptides = digestion.digestUnmodified(protein, Number($("digest-min").value), Number($("digest-max").value));
    const rows = peptides.map((p, i) => {
      const seq = OpenMS.AASequence.fromString(p);
      const row = [i + 1, p, p.length, seq.getMonoWeight().toFixed(4), seq.getMZ(2).toFixed(4), seq.getMZ(3).toFixed(4)];
      seq.delete();
      return row;
    });
    fillTable($("digest-table"), rows);
    $("digest-summary").textContent = `${peptides.length} peptides from ${protein.length} residues.`;
  } catch (e) {
    $("digest-summary").textContent = e.message;
  } finally {
    digestion.delete();
  }
}

$("digest-form").addEventListener("submit", digest);
digest();

function fillTable(table, rows, onClick) {
  const body = table.tBodies[0];
  body.replaceChildren(...rows.map((cells, index) => {
    const tr = document.createElement("tr");
    for (const c of cells) {
      const td = document.createElement("td");
      td.textContent = c ?? "";
      tr.append(td);
    }
    if (onClick) tr.addEventListener("click", () => onClick(index, tr));
    return tr;
  }));
  table.hidden = false;
}

//------------------------------------------------------------------------------
// mzML viewer

let experiment = null;
let fileName = "";
const MAX_ROWS = 5000;

async function openFile(file) {
  status(`Reading ${file.name}…`);
  await repaint();
  const path = `/tmp/${file.name}`;
  OpenMS.FS.writeFile(path, new Uint8Array(await file.arrayBuffer()));
  experiment?.delete();
  experiment = new OpenMS.MSExperiment();
  const started = performance.now();
  try {
    const handler = new OpenMS.FileHandler();
    handler.loadExperiment(path, experiment);
    handler.delete();
  } catch (e) {
    status(`Could not read ${file.name}: ${e.message}`);
    return;
  } finally {
    OpenMS.FS.unlink(path);
  }
  fileName = file.name.replace(/\.(mzML|mzXML|mgf)(\.gz)?$/i, "");
  status(`Loaded ${file.name} in ${((performance.now() - started) / 1000).toFixed(1)} s.`);
  showExperiment();
}

function showExperiment() {
  experiment.updateRanges();
  const summary = experiment.getSpectrumSummary();
  $("mzml-summary").textContent = experiment.getNrSpectra() > 0
    ? `${experiment.getNrSpectra()} spectra (MS levels ${experiment.getMSLevels().join(", ")}), ` +
      `${experiment.getNrChromatograms()} chromatograms, ${experiment.getNrPeaks()} peaks; ` +
      `RT ${experiment.getMinRT().toFixed(1)}–${experiment.getMaxRT().toFixed(1)} s, ` +
      `m/z ${experiment.getMinMZ().toFixed(2)}–${experiment.getMaxMZ().toFixed(2)}.` +
      (summary.length > MAX_ROWS ? ` The table shows the first ${MAX_ROWS} spectra.` : "")
    : `${experiment.getNrChromatograms()} chromatograms, no spectra.`;
  $("mzml-actions").hidden = false;

  const tic = experiment.calculateTIC(0, 1);
  const { rt, intensity } = tic.getPeaks();
  tic.delete();
  $("tic-title").hidden = $("mzml-tic").hidden = rt.length === 0;
  linePlot($("mzml-tic"), rt, intensity, "RT (s)");

  fillTable($("mzml-table"), summary.slice(0, MAX_ROWS).map((s) => [
    s.index, s.nativeID, s.rt.toFixed(2), s.msLevel, s.size, s.basePeakMZ?.toFixed(4),
    s.precursorMZ?.toFixed(4), s.precursorCharge || "",
  ]), (index, tr) => {
    document.querySelectorAll("#mzml-table tr.selected").forEach((r) => r.classList.remove("selected"));
    tr.classList.add("selected");
    showSpectrum(index);
  });
  if (summary.length > 0) showSpectrum(summary.findIndex((s) => s.size > 0));
}

function showSpectrum(index) {
  if (index < 0) return;
  const spectrum = experiment.getSpectrum(index);
  const { mz, intensity } = spectrum.getPeaks();
  $("spectrum-title").textContent = `Spectrum ${index} (MS${spectrum.getMSLevel()}, RT ${spectrum.getRT().toFixed(2)} s, ${spectrum.getType().toLowerCase()})`;
  $("spectrum-title").hidden = $("mzml-spectrum").hidden = false;
  if (spectrum.getType() === "Profile" || mz.length > 2000) {
    linePlot($("mzml-spectrum"), mz, intensity, "m/z");
  } else {
    stickPlot($("mzml-spectrum"), mz, intensity);
  }
  spectrum.delete();
}

async function busy(button, message, work) {
  const buttons = document.querySelectorAll("#mzml-actions button");
  buttons.forEach((b) => { b.disabled = true; });
  status(message);
  await repaint();
  const started = performance.now();
  try {
    const result = work();
    status(`${result} (${((performance.now() - started) / 1000).toFixed(1)} s)`);
  } catch (e) {
    status(`${button.textContent} failed: ${e.message}`);
  } finally {
    buttons.forEach((b) => { b.disabled = false; });
  }
}

$("action-centroid").addEventListener("click", (event) => busy(event.target, "Centroiding…", () => {
  const picker = new OpenMS.PeakPickerHiRes();
  const centroided = new OpenMS.MSExperiment();
  picker.pickExperiment(experiment, centroided, false);
  picker.delete();
  experiment.delete();
  experiment = centroided;
  showExperiment();
  return `Centroided ${experiment.getNrSpectra()} spectra.`;
}));

$("action-features").addEventListener("click", (event) => busy(event.target, "Detecting features…", () => {
  const traces = new OpenMS.MassTraceList();
  const split = new OpenMS.MassTraceList();
  const features = new OpenMS.FeatureMap();
  const mtd = new OpenMS.MassTraceDetection();
  const epd = new OpenMS.ElutionPeakDetection();
  const ffm = new OpenMS.FeatureFindingMetabo();
  try {
    mtd.run(experiment, traces);
    epd.detectPeaks(traces, split);
    ffm.run(split, features);
    features.sortByIntensity(true);
    const t = features.getFeatureTable();
    fillTable($("features-table"), Array.from(t.rt, (rt, i) => [
      i + 1, rt.toFixed(2), t.mz[i].toFixed(5), t.charge[i], t.intensity[i].toExponential(3), t.quality[i].toFixed(3),
    ]).slice(0, MAX_ROWS));
    $("features-title").hidden = false;
    return `${features.size()} features from ${split.size()} mass traces.`;
  } finally {
    [traces, split, features, mtd, epd, ffm].forEach((o) => o.delete());
  }
}));

$("action-download").addEventListener("click", () => {
  const text = new OpenMS.MzMLFile().storeBuffer(experiment);
  const url = URL.createObjectURL(new Blob([text], { type: "application/xml" }));
  const a = Object.assign(document.createElement("a"), { href: url, download: `${fileName || "experiment"}.mzML` });
  a.click();
  URL.revokeObjectURL(url);
});

const drop = $("drop-zone");
$("mzml-file").addEventListener("change", (event) => event.target.files[0] && openFile(event.target.files[0]));
drop.addEventListener("dragover", (event) => { event.preventDefault(); drop.classList.add("over"); });
drop.addEventListener("dragleave", () => drop.classList.remove("over"));
drop.addEventListener("drop", (event) => {
  event.preventDefault();
  drop.classList.remove("over");
  if (event.dataTransfer.files[0]) openFile(event.dataTransfer.files[0]);
});
