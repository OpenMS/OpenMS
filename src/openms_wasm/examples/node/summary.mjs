// Prints a summary of a peak file (mzML, mzXML, MGF) with the OpenMS WebAssembly module.
//
//   node examples/node/summary.mjs <file> [--centroid]
//
// Run from the package directory (or adjust the import to "@openms/openms-wasm").

import path from "node:path";
import createOpenMS from "../../openms.mjs";

const [file, ...flags] = process.argv.slice(2);
if (!file) {
  console.error("usage: node summary.mjs <mzML|mzXML|mgf file> [--centroid]");
  process.exit(1);
}

const OpenMS = await createOpenMS();
OpenMS.setLogEnabled("info", false);

// Make the file's directory visible to OpenMS (Node.js: no copy into memory needed).
const directory = path.dirname(path.resolve(file));
OpenMS.FS.mkdir("/input");
OpenMS.FS.mount(OpenMS.FS.filesystems.NODEFS, { root: directory }, "/input");

// Objects created from JavaScript live in WebAssembly memory: delete() them when done.
const handler = new OpenMS.FileHandler();
const experiment = new OpenMS.MSExperiment();
handler.loadExperiment(`/input/${path.basename(file)}`, experiment);
handler.delete();
experiment.updateRanges();

console.log(`${path.basename(file)}: ${experiment.getNrSpectra()} spectra, ${experiment.getNrChromatograms()} chromatograms, ` +
            `${experiment.getNrPeaks()} peaks`);
console.log(`MS levels ${experiment.getMSLevels().join(", ")}; RT ${experiment.getMinRT().toFixed(1)}-${experiment.getMaxRT().toFixed(1)} s; ` +
            `m/z ${experiment.getMinMZ().toFixed(3)}-${experiment.getMaxMZ().toFixed(3)}`);

if (flags.includes("--centroid")) {
  const picker = new OpenMS.PeakPickerHiRes();
  const centroided = new OpenMS.MSExperiment();
  picker.pickExperiment(experiment, centroided);
  console.log(`centroided: ${centroided.getNrPeaks()} peaks`);
  picker.delete();
  centroided.delete();
}

const rows = experiment.getSpectrumSummary();
console.table(rows.slice(0, 10).map((s) => ({
  rt: s.rt.toFixed(2), msLevel: s.msLevel, peaks: s.size,
  basePeakMZ: s.basePeakMZ?.toFixed(4), precursorMZ: s.precursorMZ?.toFixed(4), charge: s.precursorCharge,
})));
experiment.delete();
