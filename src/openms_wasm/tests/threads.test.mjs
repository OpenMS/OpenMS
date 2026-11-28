// Threads: the multi-threaded module (OPENMS_WASM_THREADS, @openms/openms-wasm-mt) computes with
// Module.numThreads threads; the single-threaded one reports that it has none.

import { test } from "node:test";
import assert from "node:assert/strict";
import { pathToFileURL } from "node:url";
import { loadOpenMS, modulePath, testDataDir, toppDataDir, using } from "./helpers.mjs";

const OpenMS = await loadOpenMS();
const info = OpenMS.getThreadingInfo();
const withoutThreads = !info.threaded && "built without threads";

/// A further instance of the module, with the test data mounted like loadOpenMS() does.
async function createInstance(options) {
  const { default: createOpenMS } = await import(pathToFileURL(modulePath).href);
  const instance = await createOpenMS({ print: () => {}, printErr: () => {}, ...options });
  instance.mountHostDirectory("/data", testDataDir);
  instance.mountHostDirectory("/topp", toppDataDir);
  return instance;
}

/// Features of FeatureFindingMetabo (OpenMP-parallel, like the mzML decoding before it) on a larger input.
function findFeatures(m) {
  return using([new m.MSExperiment(), new m.MassTraceList(), new m.MassTraceList(), new m.FeatureMap(),
                new m.MassTraceDetection(), new m.ElutionPeakDetection(), new m.FeatureFindingMetabo()],
    (input, traces, split, features, mtd, epd, ffm) => {
      using([new m.MzMLFile()], (f) => f.load("/topp/FeatureFinderMetabo_3_input.mzML", input));
      mtd.run(input, traces);
      epd.detectPeaks(traces, split);
      ffm.run(split, features);
      return features.getFeatureTable();
    });
}

test("getThreadingInfo", () => {
  if (info.threaded) {
    assert.equal(info.available, true, "Node.js can start threads");
    assert.ok(info.numThreads >= 1);
  } else {
    assert.deepEqual(info, { threaded: false, available: false, numThreads: 1 });
  }
});

test("numThreads sets the number of threads", { skip: withoutThreads }, async () => {
  for (const n of [1, 3]) {
    const m = await createInstance({ numThreads: n });
    assert.equal(m.getThreadingInfo().numThreads, n);
  }
  await assert.rejects(createInstance({ numThreads: 0 }), /numThreads must be a positive integer/);
});

/// Peaks of all spectra of an mzML file (the binary data is decoded in parallel).
function loadPeaks(m) {
  return using([new m.MSExperiment(), new m.MzMLFile()], (exp, file) => {
    file.load("/topp/InternalCalibration_1_BSA1.mzML", exp);
    const peaks = [];
    for (let i = 0; i < exp.getNrSpectra(); ++i) {
      using([exp.getSpectrum(i)], (s) => peaks.push(s.getPeaks()));
    }
    return peaks;
  });
}

test("results do not depend on the number of threads", { skip: withoutThreads }, async () => {
  const one = await createInstance({ numThreads: 1 });
  const four = await createInstance({ numThreads: 4 });

  const single = findFeatures(one);
  const multi = findFeatures(four);
  assert.equal(single.mz.length, 13);
  for (const column of ["rt", "mz", "intensity", "charge", "quality", "width"]) {
    assert.deepEqual(multi[column], single[column], column);
  }

  const peaksOne = loadPeaks(one);
  const peaksFour = loadPeaks(four);
  assert.equal(peaksOne.length, 1684);
  assert.deepEqual(peaksFour, peaksOne);
});
