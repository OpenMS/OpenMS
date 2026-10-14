// Searches a spectrum file against a protein database with ProSE, the OpenMS search engine, in
// WebAssembly, and prints what the search found, how long it took and how much memory it used.
//
//   node examples/node/prose.mjs <spectra.mzML> <database.fasta> [options]
//     --threads <n>           threads of the threaded build (default: OMP_NUM_THREADS or all cores)
//     --set <name>=<value>    a ProSE parameter, e.g. --set "fragment:mass_tolerance=0.02"; lists with commas:
//                             --set "modifications:fixed=Carbamidomethyl (C)"
//     --out <file.idXML>      write the identifications
//     --no-rescore            keep ProSE's scores (default: rescore the PSMs with Percolator)
//     --json                  print the summary as JSON
//
// Run from the package directory (or adjust the import to "@openms/openms-wasm").

import fs from "node:fs";
import path from "node:path";
import createOpenMS from "../../openms.mjs";

const args = process.argv.slice(2);
const positional = [];
const settings = {};
let threads, out, json = false, rescore = true;
for (let i = 0; i < args.length; ++i) {
  const a = args[i];
  if (a === "--threads") threads = Number(args[++i]);
  else if (a === "--out") out = args[++i];
  else if (a === "--json") json = true;
  else if (a === "--no-rescore") rescore = false;
  else if (a === "--set") {
    const [name, ...rest] = args[++i].split("=");
    settings[name] = rest.join("=");
  } else positional.push(a);
}
if (positional.length !== 2) {
  console.error("usage: node prose.mjs <spectra.mzML> <database.fasta> [--threads n] [--set name=value]... [--out file.idXML] [--json]");
  process.exit(1);
}

const started = performance.now();
const OpenMS = await createOpenMS({ numThreads: threads, print: () => {}, printErr: () => {} });
OpenMS.setLogEnabled("info", false);
OpenMS.setLogEnabled("warn", false);

// the files' directories, visible to OpenMS without copying (NODEFS)
const mount = (file, at) => {
  OpenMS.FS.mkdir(at);
  OpenMS.FS.mount(OpenMS.FS.filesystems.NODEFS, { root: path.dirname(path.resolve(file)) }, at);
  return `${at}/${path.basename(file)}`;
};
const spectra = mount(positional[0], "/spectra");
const database = mount(positional[1], "/database");

const prose = new OpenMS.ProSEAlgorithm();
const proteins = new OpenMS.ProteinIdentificationList();
const peptides = new OpenMS.PeptideIdentificationList();

// parameters given as strings: convert them to the type of the parameter
const current = prose.getParameters();
const entries = new Map(current.getEntries().map((e) => [e.name, e]));
current.delete();
const values = {};
for (const [name, text] of Object.entries(settings)) {
  const entry = entries.get(name);
  if (!entry) throw new Error(`Unknown ProSE parameter '${name}'`);
  const list = text === "" ? [] : text.split(",").map((s) => s.trim());
  values[name] = { "int": Number(text), "double": Number(text), "bool": text === "true", "string[]": list,
                   "int[]": list.map(Number), "double[]": list.map(Number) }[entry.type] ?? text;
}
prose.setParameterValues(values);

let lastLabel = "";
if (!json && process.stderr.isTTY) {
  prose.setProgressCallback((p) => {
    if (p.event === "start" && p.label) lastLabel = p.label;
    const percent = p.end > p.begin ? Math.round(((p.value - p.begin) / (p.end - p.begin)) * 100) : 0;
    process.stderr.write(`\r\x1b[K${lastLabel} ${p.event === "progress" ? `${percent}%` : ""}`);
  });
}

const searchStarted = performance.now();
const summary = prose.search(spectra, database, proteins, peptides, { rescore });
const searchSeconds = (performance.now() - searchStarted) / 1000;
if (!json && process.stderr.isTTY) process.stderr.write("\r\x1b[K");

if (out) {
  const file = new OpenMS.IdXMLFile();
  file.store("/tmp/out.idXML", proteins, peptides);
  file.delete();
  fs.writeFileSync(out, OpenMS.FS.readFile("/tmp/out.idXML"));
  OpenMS.FS.unlink("/tmp/out.idXML");
}

const memory = OpenMS.getMemoryInfo();
const threading = OpenMS.getThreadingInfo();
const report = {
  status: summary.status,
  threads: threading.numThreads,
  threaded: threading.threaded,
  ms2Spectra: summary.run.ms2Spectra,
  psms: peptides.size(),
  targetPSMs: summary.run.targetPSMs,
  uniquePeptides: summary.run.uniquePeptides,
  uniqueProteins: summary.run.uniqueProteins,
  indexedPeptides: summary.search.indexedPeptides,
  indexedFragments: summary.search.indexedFragments,
  secondsIndexBuild: summary.search.secondsIndexBuild,
  secondsSearch: summary.run.secondsSearch,
  secondsTotal: searchSeconds,
  secondsWithStartup: (performance.now() - started) / 1000,
  heapMB: memory.heapSize / 2 ** 20,
  openSearch: summary.openSearch,
  percolator: summary.percolator,
};
for (const o of [prose, proteins, peptides]) o.delete();

if (json) {
  console.log(JSON.stringify({ ...report, summary }, null, 2));
} else {
  console.log(`${path.basename(positional[0])} against ${path.basename(positional[1])}: ${report.status}`);
  console.log(`  ${report.ms2Spectra} MS2 spectra, ${report.psms} PSMs (${report.targetPSMs} target), ` +
              `${report.uniquePeptides} peptides, ${report.uniqueProteins} proteins`);
  console.log(`  index: ${report.indexedPeptides} peptides, ${report.indexedFragments} fragments, built in ${report.secondsIndexBuild.toFixed(1)} s`);
  console.log(`  search: ${report.secondsTotal.toFixed(1)} s on ${report.threads} thread(s) (${report.threaded ? "threaded" : "single-threaded"} build); ` +
              `memory ${report.heapMB.toFixed(0)} MB`);
}
