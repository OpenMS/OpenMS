// The search engine ProSE in WebAssembly: ProSEAlgorithm with the settings and data of the TOPP
// test TOPP_ProSE_1 (src/tests/topp), which has to give the same PSMs; its progress reports,
// statistics, peak annotations and open search; the search of several files pooled as the ProSE
// tool's -out_merged (TOPP_ProSE_13_merged), with its results as typed columns (tables()) and the
// peaks of the identified spectra (spectra()); and annotateSpectrum() and annotatePeaks() for spectrum viewers.

import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import os from "node:os";
import path from "node:path";
import { pathToFileURL } from "node:url";
import { loadOpenMS, modulePath, toppDataDir, using, assertClose } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

/// ProSE_1.ini
const prose1 = {
  "enzyme": "Trypsin",
  "decoys": "ignore",
  "precursor:mass_tolerance_lower": 5,
  "precursor:mass_tolerance_upper": 5,
  "precursor:mass_tolerance_unit": "ppm",
  "precursor:min_charge": 2,
  "precursor:max_charge": 5,
  "fragment:mass_tolerance": 0.1,
  "fragment:mass_tolerance_unit": "Da",
  "modifications:fixed": ["Carbamidomethyl (C)"],
  "modifications:variable": ["Oxidation (M)"],
  "modifications:variable_max_per_peptide": 2,
  "peptide:min_size": 7,
  "peptide:max_size": 40,
  "peptide:missed_cleavages": 1,
  "report:top_hits": 1,
};

/// Runs ProSE on the data of TOPP_ProSE_1 with @p parameters; returns the summary and the PSMs (plain objects).
function search(parameters, onProgress) {
  return using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues(parameters);
      if (onProgress) prose.setProgressCallback(onProgress);
      const summary = prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides);
      const annotations = [];
      for (let i = 0; i < peptides.size(); ++i) {
        using([peptides.get(i)], (id) => {
          for (const hit of id.getHits()) {
            annotations.push(hit.getPeakAnnotations());
            hit.delete();
          }
        });
      }
      return { summary, proteins: proteins.toArray(), psms: peptides.toArray(), annotations };
    });
}

test("ProSEAlgorithm reproduces TOPP_ProSE_1", () => {
  const { summary, proteins, psms, annotations } = search(prose1);
  assert.equal(summary.status, "EXECUTION_OK");
  assert.equal(summary.openSearch, false);
  assert.equal(summary.modificationAnalysis, null);
  assert.equal(proteins.length, 1);

  // ProSE_1_out.idXML
  const expected = {
    "DFASSGGYVLHLHR": 22.967891913444792,
    "IALSRPNVEVVALNDPFITNDYAAYM(Oxidation)FK": 45.369950862989675,
    "RPGADSDIGGFGGLFDLAQAGFR": 42.332335629931798,
  };
  const hits = psms.flatMap((id) => id.hits);
  assert.deepEqual(hits.map((h) => h.sequence).sort(), Object.keys(expected).sort());
  for (const hit of hits) {
    assertClose(assert, hit.score, expected[hit.sequence], 1e-6 * expected[hit.sequence], hit.sequence);
    assert.equal(hit.charge, 3);
  }
  // every PSM refers to its spectrum
  for (const id of psms) assert.ok(id.spectrumReference.length > 0);

  // the fragment ions ProSE matched (annotate:PSM includes fragment annotations by default)
  assert.equal(annotations.length, 3);
  for (const a of annotations) {
    assert.ok(a.length > 0, "fragment annotations");
    assert.ok(a.every((p) => /^[abcxyz]/.test(p.annotation) && p.charge >= 1 && p.mz > 0));
  }

  // statistics
  assert.ok(summary.run.ms2Spectra > 0);
  assert.equal(summary.run.matchedSpectra, 3);
  assert.equal(summary.run.targetPSMs, 3);
  assert.ok(summary.search.indexedPeptides > 0);
  assert.equal(summary.search.enzyme, "Trypsin");
  assert.ok(summary.search.secondsTotal >= 0);
});

test("ProSEAlgorithm.annotatePsmQValues: q-values without filtering", () => {
  using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues({ ...prose1, "decoys": "generate" });
      const summary = prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides);
      assert.equal(summary.haveDecoys, true);
      const before = peptides.toArray().filter((id) => id.hits.length > 0);
      assert.ok(before.length > 0);
      assert.ok(before.every((id) => id.scoreType !== "q-value"));
      prose.annotatePsmQValues(peptides);
      const after = peptides.toArray().filter((id) => id.hits.length > 0);
      assert.equal(after.length, before.length, "nothing is filtered");
      for (const id of after) {
        assert.equal(id.scoreType, "q-value");
        const q = id.hits[0].score;
        assert.ok(q >= 0 && q <= 1, `q-value ${q}`);
        assert.ok(Object.keys(id.hits[0].metaValues).some((k) => /hyperscore\)_score$/.test(k)), "the search score is kept");
      }
    });
});

test("ProSEAlgorithm: Percolator rescoring, in-process", () => {
  // too few PSMs: skipped, with the reason; the search engine scores stay
  const small = using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues({ ...prose1, "decoys": "generate" });
      return prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides, { rescore: true });
    });
  assert.equal(small.percolator.applied, false);
  assert.match(small.percolator.skipped, /100/);

  // on by default (the parameter 'rescore'); { rescore: false } turns it off, as does the parameter
  for (const [parameters, options, expected] of [[prose1, undefined, "object"], [prose1, { rescore: false }, "null"],
                                                 [{ ...prose1, "rescore": false }, undefined, "null"]]) {
    const summary = using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
      (prose, proteins, peptides) => {
        prose.setParameterValues(parameters);
        return options === undefined
          ? prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides)
          : prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides, options);
      });
    assert.equal(summary.percolator === null ? "null" : typeof summary.percolator, expected,
      JSON.stringify({ rescore: parameters.rescore, options }));
  }

  // the data and settings of TOPP_ProSE_11_rescore (ProSE_6_fdr.ini)
  const settings = {
    "enzyme": "Trypsin", "decoys": "generate",
    "precursor:mass_tolerance_lower": 20, "precursor:mass_tolerance_upper": 20, "precursor:mass_tolerance_unit": "ppm",
    "fragment:mass_tolerance": 10, "fragment:mass_tolerance_unit": "ppm",
    "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
    "FDR:PSM": 0.05, "FDR:protein": 0.05,
  };
  // { rescore: false } sets the parameter 'rescore' for the call: ProSE's own scores (q-values of its FDR:PSM)
  using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues(settings);
      const summary = prose.search("/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML", "/topp/THIRDPARTY/DatabaseSuitability_database.fasta",
        proteins, peptides, { rescore: false });
      assert.equal(summary.status, "EXECUTION_OK");
      assert.equal(summary.percolator, null);
      assert.equal(proteins.toArray()[0].searchEngine, "ProSE", "not rescored");
      assert.ok(peptides.toArray().some((id) => id.hits.length > 0));
      assert.equal(prose.getParameterValues()["rescore"], true, "the parameter is as set");
    });
  // rescored, q-values as scores, filtered at FDR:PSM
  using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues(settings);
      const summary = prose.search("/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML", "/topp/THIRDPARTY/DatabaseSuitability_database.fasta",
        proteins, peptides, { rescore: true });
      assert.equal(summary.status, "EXECUTION_OK");
      assert.equal(summary.percolator.applied, true, summary.percolator.skipped ?? "");
      const ids = peptides.toArray().filter((id) => id.hits.length > 0);
      assert.ok(ids.length > 0);
      // the statistics describe the PSMs returned (after rescoring and FDR), not the search before
      const decoys = ids.filter((id) => id.hits[0].isDecoy).length;
      assert.equal(summary.run.matchedSpectra, ids.length);
      assert.equal(summary.run.targetPSMs, ids.length - decoys);
      assert.equal(summary.run.decoyPSMs, decoys);
      assert.equal(summary.run.fdrApplied, true);
      assert.ok(summary.run.achievedPsmFDR <= 0.05);
      for (const id of ids) {
        assert.equal(id.scoreType, "q-value");
        assert.ok(id.hits[0].score <= 0.05, `q-value ${id.hits[0].score} passes FDR:PSM`);
      }
      assert.equal(proteins.toArray()[0].searchEngine, "Percolator");
      // the parameters are as set
      assert.equal(prose.getParameterValues()["FDR:PSM"], 0.05);
    });
});

// A module with ONNX Runtime (WITH_ONNX): ProSE's PeptDeep features, with the models of the package
// (models/ next to openms.mjs). The PeptDeep class tests check the predictions themselves.
test("ProSEAlgorithm: PeptDeep features for Percolator (ONNX Runtime)", (t) => {
  const models = path.join(path.dirname(modulePath), "models");
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  assert.ok(fs.existsSync(path.join(models, "peptdeep_ms2_dynamic.onnx")), `the models are in ${models}`);
  OpenMS.mountHostDirectory("/models", models);
  try {
    using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
      (prose, proteins, peptides) => {
        prose.setParameterValues({
          "enzyme": "Trypsin", "decoys": "generate",
          "precursor:mass_tolerance_lower": 20, "precursor:mass_tolerance_upper": 20, "precursor:mass_tolerance_unit": "ppm",
          "fragment:mass_tolerance": 10, "fragment:mass_tolerance_unit": "ppm",
          "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
          "FDR:PSM": 0.01, "FDR:protein": 0,
          "peptdeep:enable": true,
          "peptdeep:ms2_model": "/models/peptdeep_ms2_dynamic.onnx", "peptdeep:rt_model": "/models/peptdeep_rt_dynamic.onnx",
        });
        const summary = prose.search("/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML", "/topp/THIRDPARTY/DatabaseSuitability_database.fasta",
          proteins, peptides, { rescore: true });
        assert.equal(summary.status, "EXECUTION_OK");
        assert.equal(summary.percolator.applied, true, summary.percolator.skipped ?? "");
        const ids = peptides.toArray().filter((id) => id.hits.length > 0);
        assert.ok(ids.length > 500, `${ids.length} PSMs`);
        const features = ["ms2_cosine", "ms2_spectral_angle", "ms2_pearson", "ms2_frac_pred_found", "rt_abs_error"];
        for (const id of ids) {
          const meta = id.hits[0].metaValues;
          for (const f of features) assert.equal(typeof meta[f], "number", `${f} of ${id.hits[0].sequence}`);
          assert.ok(meta.ms2_cosine >= -1 && meta.ms2_cosine <= 1);
        }
        // the predictions agree with the spectra of most target PSMs
        const cosines = ids.filter((id) => !id.hits[0].isDecoy).map((id) => id.hits[0].metaValues.ms2_cosine).sort((a, b) => a - b);
        assert.ok(cosines[Math.floor(cosines.length / 2)] > 0.5, `median cosine ${cosines[Math.floor(cosines.length / 2)]}`);
      });
  } finally {
    OpenMS.unmount("/models");
  }
});

/// The progress step of searchFiles() that reads the spectra of the PSMs (keepSpectra), after the steps of the files
const readingSpectra = "Reading the spectra of the PSMs...";

/// Checks progress @p reports of a search of the files with the steps @p fileSteps (their labels, in input order),
/// followed by the outermost steps @p laterSteps outside the files:
/// - the outermost steps (depth 0) are those of the files, each with the file's position (file), then @p laterSteps,
///   with file -1;
/// - every report of a step carries its label, range and depth; its values lie in the range and never decrease;
/// - every report carries the position of the file whose step encloses it.
function assertProgressSteps(reports, fileSteps, laterSteps = []) {
  assert.deepEqual(reports.filter((r) => r.event === "start" && r.depth === 0).map((r) => [r.label, r.file]),
    [...fileSteps.map((label, i) => [label, i]), ...laterSteps.map((label) => [label, -1])]);
  const open = [];
  let file = -1;
  for (const r of reports) {
    assert.ok(["start", "progress", "end"].includes(r.event), r.event);
    if (r.event === "start") {
      assert.equal(r.depth, open.length, `depth of ${r.label}`);
      if (open.length === 0) file = fileSteps.indexOf(r.label, file + 1);
      open.push({ label: r.label, begin: r.begin, end: r.end, value: r.value });
    }
    const step = open.at(-1);
    assert.ok(step !== undefined, `${r.event} report outside a step: ${JSON.stringify(r)}`);
    const where = `${r.event} of "${step.label}": ${JSON.stringify(r)}`;
    assert.equal(r.label, step.label, where);
    assert.deepEqual([r.begin, r.end], [step.begin, step.end], where);
    assert.ok(r.begin <= r.value && r.value <= r.end, `value in range, ${where}`);
    assert.ok(r.value >= step.value, `value does not decrease, ${where}`);
    step.value = r.value;
    assert.equal(r.file, file, `file, ${where}`);
    if (r.event === "end") {
      open.pop();
      assert.equal(r.depth, open.length, `depth, ${where}`);
    }
  }
  assert.equal(open.length, 0, "every step ended");
}

test("ProSEAlgorithm reports its progress", () => {
  const reports = [];
  const { summary } = search(prose1, (p) => reports.push(p));
  assert.equal(summary.status, "EXECUTION_OK");
  const started = reports.filter((r) => r.event === "start").map((r) => r.label);
  assert.ok(started.includes("Building fragment index..."), JSON.stringify(started));
  assert.ok(reports.some((r) => r.event === "end"));
  assertProgressSteps(reports, ["Searching SimpleSearchEngine_1.mzML (1 of 1)..."]);
});

// A module with ONNX Runtime: PeptDeep without its model throws inside a step of the file's step, which the search
// leaves open. The next search of the algorithm has to report its steps from the outermost level again.
test("ProSEAlgorithm: a search that throws leaves the progress of the next search intact", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      let reports = [];
      prose.setProgressCallback((r) => reports.push(r));
      const failing = { ...prose1, "peptdeep:enable": true, "peptdeep:ms2_model": "/no-such-model/ms2.onnx" };
      const steps = ["Searching SimpleSearchEngine_1.mzML (1 of 1)..."];
      for (const [run, laterSteps] of [
        // searchFiles() then reads the spectra of the PSMs
        [() => using([prose.searchFiles(["/topp/SimpleSearchEngine_1.mzML"], "/topp/SimpleSearchEngine_1.fasta", { rescore: false })], () => {}),
         [readingSpectra]],
        [() => prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides, { rescore: false }),
         []],
      ]) {
        prose.setParameterValues(failing);
        reports = [];
        assert.throws(run, /no-such-model/);
        assert.ok(reports.some((r) => r.event === "start" && r.label === "Adding PeptDeep rescoring features..."),
          "the search threw inside a step of the file's step");
        prose.setParameterValues({ "peptdeep:enable": false });
        reports = [];
        run();
        assertProgressSteps(reports, steps, laterSteps);
      }
    });
});

test("ProSEAlgorithm: progress callback has to be a function", () => {
  using([new OpenMS.ProSEAlgorithm()], (prose) => {
    assert.throws(() => prose.setProgressCallback(42), /function/);
    prose.setProgressCallback(null); // removes it
  });
});

test("ProSEAlgorithm: open search gives the modification analysis", () => {
  const { summary } = search({
    ...prose1,
    "precursor:mass_tolerance_lower": 500,
    "precursor:mass_tolerance_upper": 500,
    "precursor:mass_tolerance_unit": "Da",
  });
  assert.equal(summary.status, "EXECUTION_OK");
  assert.equal(summary.openSearch, true);
  const analysis = summary.modificationAnalysis;
  assert.ok(analysis !== null);
  assert.ok(Array.isArray(analysis.deltaMasses.entries));
  assert.ok(Array.isArray(analysis.ptms.entries));
  assert.ok(analysis.deltaMasses.totalPSMs >= 0);
});

//------------------------------------------------------------------------------
// Several files: searchFiles() and ProSEResultSet

/// ProSE_6_fdr.ini: the settings of TOPP_ProSE_13_merged (with its FDR:protein 0.05)
const prose6 = {
  "enzyme": "Trypsin",
  "decoys": "generate",
  "decoy_prefix": "DECOY_",
  "low_memory": false,
  "precursor:mass_tolerance_lower": 20,
  "precursor:mass_tolerance_upper": 20,
  "precursor:mass_tolerance_unit": "ppm",
  "precursor:min_charge": 2,
  "precursor:max_charge": 5,
  "precursor:isotope_error_min": -1,
  "precursor:isotope_error_max": 1,
  "fragment:mass_tolerance": 10,
  "fragment:mass_tolerance_unit": "ppm",
  "fragment:deisotope": "auto",
  "fragment:min_mz": 150,
  "fragment:max_mz": 2000,
  "fragment:min_ion_index": 2,
  "fragment:min_matched_ions": 5,
  "fragment:max_charge": 2,
  "modifications:fixed": ["Carbamidomethyl (C)"],
  "modifications:variable": ["Oxidation (M)"],
  "modifications:variable_max_per_peptide": 2,
  "annotate:PSM": ["ALL"],
  "peptide:min_size": 7,
  "peptide:max_size": 40,
  "peptide:missed_cleavages": 1,
  "peptide:enzyme_specificity": "full",
  "peptide:motif": "",
  "peptide:min_mass": 100,
  "peptide:max_mass": 9000,
  "report:top_hits": 1,
  "FDR:PSM": 0.05,
  "FDR:protein": 0.05,
  "scoring:max_candidates_per_spectrum": 50,
  "ions:add_y_ions": true,
  "ions:add_b_ions": true,
  "ions:add_a_ions": false,
  "ions:add_c_ions": false,
  "ions:add_x_ions": false,
  "ions:add_z_ions": false,
  "calibration:enabled": false,
  "calibration:subset_ratio": 0.1,
  "calibration:min_psms": 50,
  "database:chunk_size": 0,
};
const suitabilityDatabase = "/topp/THIRDPARTY/DatabaseSuitability_database.fasta";
const suitabilitySpectra = path.join(toppDataDir, "THIRDPARTY/DatabaseSuitability_in_spec.mzML");

/// Copies @p source (default: DatabaseSuitability_in_spec.mzML) to @p names (paths relative to a new temporary
/// directory, e.g. "0/run.mzML"), mounted at /searchfiles of module @p m while @p fn runs; @p fn gets the paths of
/// the copies there, and the directory on the host.
function withSpectraCopies(m, names, fn, source = suitabilitySpectra) {
  const dir = fs.mkdtempSync(path.join(os.tmpdir(), "openms-wasm-prose-"));
  try {
    for (const name of names) {
      fs.mkdirSync(path.dirname(path.join(dir, name)), { recursive: true });
      fs.copyFileSync(source, path.join(dir, name));
    }
    m.mountHostDirectory("/searchfiles", dir);
    try {
      return fn(names.map((name) => `/searchfiles/${name}`), dir);
    } finally {
      m.unmount("/searchfiles");
    }
  } finally {
    fs.rmSync(dir, { recursive: true, force: true });
  }
}

/// ProSEAlgorithm.searchFiles() of module @p m; returns the ProSEResultSet, which the caller deletes.
function searchFiles(m, files, options, { parameters = prose6, database = suitabilityDatabase, onProgress } = {}) {
  return using([new m.ProSEAlgorithm()], (prose) => {
    prose.setParameterValues(parameters);
    if (onProgress) prose.setProgressCallback(onProgress);
    return options === undefined ? prose.searchFiles(files, database) : prose.searchFiles(files, database, options);
  });
}

/// The identifications of result set @p set, stored as idXML at @p file and loaded again (plain objects).
function storedIdentifications(set, file) {
  set.storeIdXML(file);
  return using([new OpenMS.IdXMLFile(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (idxml, proteins, peptides) => {
      idxml.load(file, proteins, peptides);
      return { proteins: proteins.toArray(), peptides: peptides.toArray() };
    });
}

/// The counts of FileInfo's identification summary (FileInfo -out_tsv)
function fileInfoCounts({ proteins, peptides }) {
  const matched = peptides.filter((id) => id.hits.length > 0);
  return {
    runs: proteins.length,
    proteinHits: proteins.reduce((n, run) => n + run.hits.length, 0),
    matchedSpectra: matched.length,
    peptideHits: matched.reduce((n, id) => n + id.hits.length, 0),
    distinctPeptides: new Set(matched.flatMap((id) => id.hits.map((hit) => hit.sequence))).size,
  };
}

/// The counts of a FileInfo -out_tsv file of the TOPP tests
function fileInfoReference(file) {
  const lines = fs.readFileSync(path.join(toppDataDir, file), "utf8").split("\n");
  const value = (key) => {
    const line = lines.find((l) => l.startsWith(key));
    assert.ok(line, `${key} in ${file}`);
    return Number(line.split("\t").at(-1));
  };
  return {
    runs: value("general: num. of runs"),
    proteinHits: value("general: num. of protein hits"),
    matchedSpectra: value("general: num. of matched spectra"),
    peptideHits: value("general: num. of peptide hits"),
    distinctPeptides: value("general: num. of non-redundant peptide hits"),
  };
}

/// FuzzyDiff's comparison of the TOPP tests: within @p absdiff, or within a ratio of @p ratio
function assertFuzzy(actual, expected, absdiff, ratio, message) {
  const withinRatio = actual > 0 && expected > 0 && Math.max(actual, expected) / Math.min(actual, expected) <= ratio;
  assert.ok(Math.abs(actual - expected) <= absdiff || withinRatio,
    `${message}: expected ${expected} (± ${absdiff} or a ratio of ${ratio}), got ${actual}`);
}

/// The columns of the three tables of ProSEResultSet.tables() with one entry per row, and their types
const tableColumns = {
  psms: { file: Uint16Array, nativeId: Array, rt: Float64Array, mz: Float64Array, charge: Int8Array, peptide: Uint32Array,
          score: Float64Array, qValue: Float64Array, pep: Float64Array, precursorErrorPpm: Float64Array,
          isotopeError: Int8Array, deltaMass: Float64Array, decoy: Uint8Array },
  peptides: { openms: Array, proforma: Array, sequence: Array, decoy: Uint8Array },
  proteins: { accession: Array, description: Array, sequence: Array, score: Float64Array, qValue: Float64Array,
              decoy: Uint8Array, group: Int32Array },
};

/// Checks the shape of tables @p t (ProSEResultSet.tables()): the types and lengths of the columns, the offsets of the
/// modifications and protein references of the peptides, and that every reference points into its table.
function assertTableShapes(t) {
  assert.ok(Array.isArray(t.files));
  for (const [table, columns] of Object.entries(tableColumns)) {
    const n = t[table].count;
    assert.ok(Number.isInteger(n) && n >= 0, `${table}.count ${n}`);
    for (const [column, type] of Object.entries(columns)) {
      assert.ok(t[table][column] instanceof type, `${table}.${column} is a ${type.name}`);
      assert.equal(t[table][column].length, n, `${table}.${column} has ${table}.count entries`);
    }
  }
  const { psms, peptides, proteins } = t;
  assert.ok(Array.isArray(psms.featureNames) && psms.features instanceof Float64Array);
  assert.equal(psms.features.length, psms.count * psms.featureNames.length, "features: count × featureNames.length");
  assert.ok(psms.peptide.every((k) => k < peptides.count), "every PSM refers to a peptide");
  assert.ok(psms.file.every((f) => f < t.files.length), "every PSM refers to a file");
  const ranges = (offsets, columns) => {
    assert.ok(peptides[offsets] instanceof Uint32Array);
    assert.equal(peptides[offsets].length, peptides.count + 1, `${offsets}: count + 1 entries`);
    assert.equal(peptides[offsets][0], 0);
    for (let k = 0; k < peptides.count; ++k) assert.ok(peptides[offsets][k] <= peptides[offsets][k + 1], offsets);
    const n = peptides[offsets][peptides.count];
    for (const [column, type] of Object.entries(columns)) {
      assert.ok(peptides[column] instanceof type, `peptides.${column} is a ${type.name}`);
      assert.equal(peptides[column].length, n, `peptides.${column} has ${offsets}[count] entries`);
    }
  };
  ranges("modOffsets", { modPosition: Int16Array, modName: Array, modUnimod: Int32Array, modMass: Float64Array });
  ranges("proteinOffsets", { protein: Uint32Array, start: Uint32Array, end: Uint32Array });
  assert.ok(peptides.protein.every((p) => p < proteins.count), "every peptide reference points to a protein");
  for (let k = 0; k < peptides.count; ++k) {
    for (let i = peptides.modOffsets[k]; i < peptides.modOffsets[k + 1]; ++i) {
      assert.ok(peptides.modPosition[i] >= -1 && peptides.modPosition[i] <= peptides.sequence[k].length, "modPosition");
    }
  }
}

/// @p actual is @p expected (NaN for a missing number), up to the precision of idXML
function assertSame(actual, expected, message) {
  if (typeof expected !== "number" || Number.isNaN(expected)) {
    assert.ok(Number.isNaN(actual), `${message}: expected NaN, got ${actual}`);
  } else {
    assertClose(assert, actual, expected, 1e-9 * Math.max(1, Math.abs(expected)), message);
  }
}

/// The text of the idXML file @p file of the host
function idXMLText(file) {
  return fs.readFileSync(file, "utf8");
}

/// The Percolator features of the search engine that the run of idXML @p text records (its 'extra_features')
function extraFeatures(text) {
  const match = /name="extra_features" value="([^"]*)"/.exec(text);
  assert.ok(match, "extra_features in the idXML");
  return match[1].split(",");
}

/// The indistinguishable-protein groups of the run of idXML @p text, as accessions: one array per group, in order
function indistinguishableGroups(text) {
  const accessions = new Map([...text.matchAll(/<ProteinHit id="([^"]*)" accession="([^"]*)"/g)].map((m) => [m[1], m[2]]));
  const groups = [...text.matchAll(/name="indistinguishable_proteins_(\d+)" value="([^"]*)"/g)]
    .map((m) => [Number(m[1]), m[2].split(",").slice(1).map((id) => accessions.get(id))]);
  return groups.sort((a, b) => a[0] - b[0]).map(([, group]) => group);
}

/// The search of TOPP_ProSE_13_merged (two copies of one file, a.mzML and b.mzML; Percolator, the decoys removed),
/// run once for the tests that look at it: its summary, the stored identifications and the progress reports.
let pooledSearch;
function pooled() {
  pooledSearch ??= withSpectraCopies(OpenMS, ["a.mzML", "b.mzML"], (files) => {
    const reports = [];
    const set = searchFiles(OpenMS, files, { removeDecoys: true }, { onProgress: (r) => reports.push(r) });
    try {
      return { files, reports, summary: set.summary(), stored: storedIdentifications(set, "/searchfiles/merged.idXML") };
    } finally {
      set.delete();
    }
  });
  return pooledSearch;
}

test("searchFiles pools two files like ProSE -out_merged", () => {
  const { files, summary, stored } = pooled();
  assert.equal(summary.files.length, 2);
  assert.deepEqual(summary.files.map((f) => f.path), files);
  for (const f of summary.files) {
    assert.equal(f.status, "EXECUTION_OK");
    assert.equal(f.ms2Spectra, summary.files[0].ms2Spectra);
    assert.ok(f.ms2Spectra > 0 && f.targetPSMs > 0 && f.decoyPSMs >= 0, JSON.stringify(f));
    assert.equal(f.percolator.applied, true, f.percolator.skipped ?? "");
    assert.equal(f.percolator.skipped, null);
    assert.ok(f.percolator.warning === null || typeof f.percolator.warning === "string");
    assert.equal(f.peptdeepNce, null, "no PeptDeep rescoring");
    assert.ok(f.secondsSearch >= 0);
  }
  assert.equal(summary.proteinFdrApplied, true);
  assert.equal(summary.psmFdr, 0.05);
  assert.equal(summary.proteinFdr, 0.05);
  assert.equal(summary.decoyString, "DECOY_");
  assert.equal(summary.openSearch, false);
  assert.equal(summary.modificationAnalysis, null);

  // the pooled run, as FileInfo counts it, against the reference of TOPP_ProSE_13_merged with its tolerances
  const counts = fileInfoCounts(stored);
  const expected = fileInfoReference("ProSE_13_merged_info.tsv");
  assert.equal(counts.runs, expected.runs);
  assertFuzzy(counts.proteinHits, expected.proteinHits, 1, 1.01, "protein hits");
  assertFuzzy(counts.matchedSpectra, expected.matchedSpectra, 6, 1.01, "matched spectra");
  assertFuzzy(counts.peptideHits, expected.peptideHits, 6, 1.01, "peptide hits");
  assertFuzzy(counts.distinctPeptides, expected.distinctPeptides, 6, 1.01, "distinct peptides");
  // TOPP_ProSE_13_merged_structure: Percolator's run with protein q-values, and no decoys left
  assert.equal(stored.proteins[0].searchEngine, "Percolator");
  assert.equal(stored.proteins[0].scoreType, "q-value");
  assert.ok(stored.proteins[0].hits.every((h) => !h.accession.startsWith("DECOY_")));
  assert.ok(stored.peptides.every((id) => id.hits.every((h) => !h.isDecoy)));
  // the PSMs of both files
  const fileIndices = new Set(stored.peptides.filter((id) => id.hits.length > 0).map((id) => id.metaValues.id_merge_index));
  assert.deepEqual([...fileIndices].sort(), [0, 1]);
});

test("searchFiles reports one progress step per file", () => {
  const { reports, summary } = pooled();
  // then the step that reads the spectra of the PSMs (keepSpectra), over the MS2 spectra of the files
  assertProgressSteps(reports, ["Searching a.mzML (1 of 2)...", "Searching b.mzML (2 of 2)..."], [readingSpectra]);
  const reading = reports.filter((r) => r.label === readingSpectra);
  assert.deepEqual([reading[0].event, reading.at(-1).event], ["start", "end"]);
  for (const r of reading) {
    assert.deepEqual([r.depth, r.begin, r.end], [0, 0, summary.files[0].ms2Spectra + summary.files[1].ms2Spectra]);
  }
  // the steps inside a file report their own range (e.g. one per spectrum), not that of the file's step
  assert.ok(reports.some((r) => r.event !== "start" && r.end > 1), "a step with a range beyond [0, 1]");
});

test("searchFiles keeps files with the same name apart", () => {
  withSpectraCopies(OpenMS, ["0/run.mzML", "1/run.mzML"], (files) => {
    const reports = [];
    const set = searchFiles(OpenMS, files, { rescore: false }, { onProgress: (r) => reports.push(r) });
    try {
      const summary = set.summary();
      assert.deepEqual(summary.files.map((f) => f.path), files);
      assert.notEqual(summary.files[0].path, summary.files[1].path);
      for (const f of summary.files) {
        assert.equal(f.status, "EXECUTION_OK");
        assert.equal(f.percolator, null, "not rescored");
      }
      // the progress tells them apart by position
      assert.deepEqual(reports.filter((r) => r.event === "start" && r.depth === 0).map((r) => [r.label, r.file]),
        [["Searching run.mzML (1 of 2)...", 0], ["Searching run.mzML (2 of 2)...", 1], [readingSpectra, -1]]);
      // PSMs of both files; the decoys are kept by default
      const { proteins, peptides } = storedIdentifications(set, "/searchfiles/merged.idXML");
      const ids = peptides.filter((id) => id.hits.length > 0);
      assert.deepEqual([...new Set(ids.map((id) => id.metaValues.id_merge_index))].sort(), [0, 1]);
      assert.equal(summary.proteinFdrApplied, true);
      assert.ok(ids.some((id) => id.hits[0].isDecoy), "decoy PSMs");
      assert.ok(proteins[0].hits.some((h) => h.accession.startsWith("DECOY_")), "decoy proteins");
      // the tables: the files in input order, and the file of each PSM by its position
      const t = set.tables();
      assertTableShapes(t);
      assert.deepEqual(t.files, files);
      assert.equal(t.psms.count, ids.length);
      assert.deepEqual(Array.from(t.psms.file), ids.map((id) => id.metaValues.id_merge_index));
    } finally {
      set.delete();
    }
  });
});

test("searchFiles: one file, open search", () => {
  const set = searchFiles(OpenMS, ["/topp/SimpleSearchEngine_1.mzML"], undefined, {
    parameters: { ...prose1, "precursor:mass_tolerance_lower": 500, "precursor:mass_tolerance_upper": 500,
                  "precursor:mass_tolerance_unit": "Da" },
    database: "/topp/SimpleSearchEngine_1.fasta",
  });
  try {
    const summary = set.summary();
    assert.equal(summary.files.length, 1);
    assert.equal(summary.files[0].status, "EXECUTION_OK");
    assert.equal(summary.openSearch, true);
    assert.ok(summary.modificationAnalysis !== null, "the modification analysis of the one file");
    assert.ok(Array.isArray(summary.modificationAnalysis.deltaMasses.entries));
    assert.equal(summary.proteinFdrApplied, false, "decoys: ignore");
    // the mass shift of every PSM; no q-values and no protein FDR without decoys
    const t = set.tables();
    assertTableShapes(t);
    assert.ok(t.psms.count > 0);
    assert.ok(t.psms.deltaMass.every(Number.isFinite), `delta masses ${t.psms.deltaMass}`);
    // the shift is a value of the PSM: the ProForma of its peptide never has it as unlocalized mass ("[+Δ]?SEQ")
    for (let i = 0; i < t.psms.count; ++i) {
      const proforma = t.peptides.proforma[t.psms.peptide[i]];
      assert.ok(!/^\[[^\]]*\]\?/.test(proforma), `no unlocalized shift in ${proforma}`);
      assert.equal(proforma.replace(/\[[^\]]*\]-?/g, ""), t.peptides.sequence[t.psms.peptide[i]], proforma);
    }
    assert.ok(t.psms.qValue.every(Number.isNaN) && t.psms.pep.every(Number.isNaN));
    assert.ok(t.proteins.count > 0 && t.proteins.qValue.every(Number.isNaN));
  } finally {
    set.delete();
  }
  assert.throws(() => searchFiles(OpenMS, [], undefined), /at least one/);
  assert.throws(() => searchFiles(OpenMS, ["/topp/no-such-file.mzML"], undefined), /no-such-file/);
  // the PSMs tell their files apart by path
  assert.throws(() => searchFiles(OpenMS, ["/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.mzML"], undefined),
    /SimpleSearchEngine_1\.mzML.*twice/);
});

test("searchFiles: open search of two files, with the modification analysis of all PSMs", () => {
  // TOPP_ProSE_12_open (ProSE_6_fdr.ini, an open precursor window); its PSMs have mass shifts, unlike those of the
  // open search of TOPP_ProSE_1's data
  const open = { ...prose6, "rescore": false, "precursor:mass_tolerance_lower": 100, "precursor:mass_tolerance_upper": 300,
                 "precursor:mass_tolerance_unit": "Da", "FDR:PSM": 0.0005, "FDR:protein": 0 };
  const analysis = (files) => {
    const set = searchFiles(OpenMS, files, undefined, { parameters: open });
    try {
      const summary = set.summary();
      assert.equal(summary.files.length, files.length);
      assert.ok(summary.files.every((f) => f.status === "EXECUTION_OK"));
      assert.equal(summary.openSearch, true);
      assert.ok(summary.modificationAnalysis !== null);
      return summary.modificationAnalysis.deltaMasses;
    } finally {
      set.delete();
    }
  };
  withSpectraCopies(OpenMS, ["a.mzML", "b.mzML"], (files) => {
    const one = analysis(files.slice(0, 1));
    const both = analysis(files);
    assert.ok(both.entries.length > 0, "delta masses");
    // the analysis of the PSMs of both files (copies of one)
    assert.ok(one.totalPSMs > 0);
    assert.equal(both.totalPSMs, 2 * one.totalPSMs);
    assert.deepEqual(both.entries.map((e) => [e.deltaMass, e.count]), one.entries.map((e) => [e.deltaMass, 2 * e.count]));
  });
});

test("searchFiles: files without PSMs give empty results", () => {
  // the spectra of TOPP_ProSE_1 match no peptide of the DatabaseSuitability database
  const reports = [];
  const set = searchFiles(OpenMS, ["/topp/SimpleSearchEngine_1.mzML"], undefined, { onProgress: (r) => reports.push(r) });
  try {
    // no spectra to read
    assertProgressSteps(reports, ["Searching SimpleSearchEngine_1.mzML (1 of 1)..."]);
    const summary = set.summary();
    assert.equal(summary.files[0].status, "EXECUTION_OK");
    assert.equal(summary.files[0].targetPSMs, 0);
    assert.equal(summary.files[0].percolator.applied, false);
    assert.match(summary.files[0].percolator.skipped, /100/);
    assert.equal(summary.proteinFdrApplied, false);
    withSpectraCopies(OpenMS, [], () => {
      const { proteins, peptides } = storedIdentifications(set, "/searchfiles/empty.idXML");
      assert.equal(peptides.length, 0);
      assert.ok(proteins.every((run) => run.hits.length === 0));
    });
    const t = set.tables();
    assertTableShapes(t);
    assert.deepEqual(t.files, ["/topp/SimpleSearchEngine_1.mzML"]);
    assert.deepEqual([t.psms.count, t.peptides.count, t.proteins.count], [0, 0, 0]);
    assert.deepEqual(t.psms.featureNames, ["ln(hyperscore)"]);
    assert.equal(assertSpectraShapes(set.spectra(), 0), 0);
  } finally {
    set.delete();
  }
});

test("searchFiles: files that Percolator could not all rescore are not pooled", () => {
  withSpectraCopies(OpenMS, ["a.mzML"], ([a]) => {
    using([new OpenMS.ProSEAlgorithm()], (prose) => {
      prose.setParameterValues({ ...prose6, "rescore": false });
      // too few PSMs in SimpleSearchEngine_1.mzML to rescore (TOPP_ProSE_11_rescore_mixed)
      assert.throws(() => prose.searchFiles([a, "/topp/SimpleSearchEngine_1.mzML"], suitabilityDatabase, { rescore: true }),
        /rescored only some.*SimpleSearchEngine_1\.mzML.*rescore: false/s);
      assert.equal(prose.getParameterValues()["rescore"], false, "the parameter is restored");
    });
  });
});

test("tables() is consistent with the result set", () => {
  // one file, rescored with Percolator; the decoys are kept (default)
  withSpectraCopies(OpenMS, ["a.mzML"], (files, dir) => {
    const set = searchFiles(OpenMS, files, undefined);
    try {
      assert.equal(set.summary().files[0].percolator.applied, true);
      const t = set.tables();
      assertTableShapes(t);
      const { psms, peptides, proteins } = t;
      const { proteins: runs, peptides: ids } = storedIdentifications(set, "/searchfiles/merged.idXML");
      const text = idXMLText(path.join(dir, "merged.idXML"));
      assert.deepEqual(t.files, files);

      // One row per identified spectrum (its best hit), in the order of the stored identifications
      const rows = ids.filter((id) => id.hits.length > 0);
      assert.equal(psms.count, rows.length);
      assert.ok(psms.count > 500, `${psms.count} PSMs`);
      const names = psms.featureNames;
      assert.ok(names.includes("ln(hyperscore)"), JSON.stringify(names));
      // the Percolator features of the search engine ('score' is its ln(hyperscore))
      for (const name of extraFeatures(text)) {
        if (name !== "score") assert.ok(names.includes(name), `feature ${name} in ${JSON.stringify(names)}`);
      }
      // the standard .pin features Percolator trains on (precursor charges 2 to 5), which are not meta values of the
      // stored hits: Percolator computes them
      const pin = ["CalcMass", "mass", "peplen", "charge2", "charge3", "charge4", "charge5", "enzN", "enzC", "enzInt", "dm", "absdm"];
      for (const name of pin) assert.ok(names.includes(name), `feature ${name} in ${JSON.stringify(names)}`);
      const feature = (i, name) => psms.features[i * names.length + names.indexOf(name)];
      rows.forEach((id, i) => {
        const hit = id.hits[0];
        const where = `PSM ${i} (${id.spectrumReference})`;
        assert.equal(psms.file[i], id.metaValues.id_merge_index, where);
        assert.equal(psms.nativeId[i], id.spectrumReference, where);
        assertSame(psms.rt[i], id.rt ?? NaN, `rt of ${where}`);
        assertSame(psms.mz[i], id.mz ?? NaN, `mz of ${where}`);
        assert.equal(psms.charge[i], hit.charge, where);
        assert.equal(peptides.openms[psms.peptide[i]], hit.sequence, where);
        assert.equal(psms.decoy[i], hit.isDecoy ? 1 : 0, where);
        // Percolator: q-values as scores, its SVM score and PEP as meta values
        assert.equal(id.scoreType, "q-value");
        assertSame(psms.qValue[i], hit.score, `qValue of ${where}`);
        assertSame(psms.score[i], hit.metaValues["MS:1001492"], `score of ${where}`);
        assertSame(psms.pep[i], hit.metaValues["MS:1001493"], `pep of ${where}`);
        assertSame(psms.precursorErrorPpm[i], hit.metaValues.precursor_mz_error_ppm, `precursorErrorPpm of ${where}`);
        assert.equal(psms.isotopeError[i], hit.metaValues.isotope_error, where);
        assert.ok(Number.isNaN(psms.deltaMass[i]), `no delta mass in a closed search, ${where}`);
        names.forEach((name, j) => {
          if (!pin.includes(name)) assertSame(psms.features[i * names.length + j], hit.metaValues[name], `${name} of ${where}`);
        });
        assert.ok(Number.isFinite(feature(i, "ln(hyperscore)")), `ln(hyperscore) of ${where}`);
        for (const name of pin) assert.ok(Number.isFinite(feature(i, name)), `${name} of ${where}: ${feature(i, name)}`);
        const residues = peptides.sequence[psms.peptide[i]];
        assert.equal(feature(i, "peplen"), residues.length, where);
        for (const z of [2, 3, 4, 5]) assert.equal(feature(i, `charge${z}`), hit.charge === z ? 1 : 0, where);
        assertSame(feature(i, "mass"), id.mz, `mass of ${where}`);
        assertSame(feature(i, "absdm"), Math.abs(feature(i, "dm")), `absdm of ${where}`);
        for (const name of ["enzN", "enzC"]) assert.ok([0, 1].includes(feature(i, name)), `${name} of ${where}`);
        assert.ok(Number.isInteger(feature(i, "enzInt")) && feature(i, "enzInt") >= 0, where);
      });
      // tryptic peptides: enzymatic termini, mostly without internal cleavage sites
      assert.ok(rows.filter((_, i) => feature(i, "enzN") === 1 && feature(i, "enzC") === 1).length > 0.9 * rows.length);
      assert.ok(psms.decoy.includes(1) && psms.decoy.includes(0), "target and decoy PSMs");

      // Peptides: unique by their OpenMS string; ProForma is the sequence with its modifications in brackets
      assert.equal(new Set(peptides.openms).size, peptides.count);
      assert.deepEqual(new Set(peptides.openms), new Set(rows.map((id) => id.hits[0].sequence)));
      const residues = (proforma) => proforma.replace(/^\[[^\]]*\]-/, "").replace(/-\[[^\]]*\]$/, "").replace(/\[[^\]]*\]/g, "");
      const accessionsOf = new Map(rows.map((id) => [id.hits[0].sequence, id.hits[0].proteinAccessions]));
      for (let k = 0; k < peptides.count; ++k) {
        const where = `peptide ${peptides.openms[k]}`;
        assert.equal(residues(peptides.proforma[k]), peptides.sequence[k], where);
        assert.equal(peptides.modOffsets[k + 1] - peptides.modOffsets[k], (peptides.openms[k].match(/\(/g) ?? []).length, where);
        for (let i = peptides.modOffsets[k]; i < peptides.modOffsets[k + 1]; ++i) {
          // carbamidomethyl-C (fixed) and oxidized M (variable)
          const [residue, name, unimod] = peptides.modName[i] === "Oxidation" ? ["M", "Oxidation", 35] : ["C", "Carbamidomethyl", 4];
          assert.equal(peptides.sequence[k][peptides.modPosition[i]], residue, where);
          assert.equal(peptides.modName[i], name, where);
          assert.equal(peptides.modUnimod[i], unimod, where);
        }
        // every reference is the stretch of the peptide in the protein's sequence, decoys included
        const refs = [];
        for (let r = peptides.proteinOffsets[k]; r < peptides.proteinOffsets[k + 1]; ++r) {
          const p = peptides.protein[r];
          assert.equal(proteins.sequence[p].slice(peptides.start[r], peptides.end[r]), peptides.sequence[k],
            `${where} in ${proteins.accession[p]} at ${peptides.start[r]}-${peptides.end[r]}`);
          refs.push(proteins.accession[p]);
        }
        assert.deepEqual([...new Set(refs)].sort(), [...new Set(accessionsOf.get(peptides.openms[k]))].sort(), where);
        const decoys = rows.filter((id) => id.hits[0].sequence === peptides.openms[k]).map((id) => id.hits[0].isDecoy);
        assert.ok(decoys.every((d) => d === (peptides.decoy[k] === 1)), `decoy flag of ${where}`);
      }
      assert.ok(peptides.proteinOffsets[peptides.count] >= peptides.count, "the peptides have proteins");

      // Proteins: the hits of the pooled run, with the sequences of the FASTA (decoys: as ProSE generated them)
      const hits = runs[0].hits;
      assert.deepEqual(proteins.accession, hits.map((h) => h.accession));
      const fasta = new Map(using([new OpenMS.FASTAFile()], (f) => f.load(suitabilityDatabase)).map((e) => [e.identifier, e]));
      const groups = indistinguishableGroups(text);
      assert.ok(groups.length > 0);
      const groupOf = new Map(groups.flatMap((group, g) => group.map((accession) => [accession, g])));
      // the best score of the PSMs of each protein: those of its peptides (targets and decoys)
      const bestScore = [];
      for (let i = 0; i < psms.count; ++i) {
        const k = psms.peptide[i];
        for (let r = peptides.proteinOffsets[k]; r < peptides.proteinOffsets[k + 1]; ++r) {
          const p = peptides.protein[r];
          bestScore[p] = Math.max(bestScore[p] ?? -Infinity, psms.score[i]);
        }
      }
      hits.forEach((h, p) => {
        const decoy = h.accession.startsWith("DECOY_");
        const target = fasta.get(decoy ? h.accession.slice("DECOY_".length) : h.accession);
        assert.ok(target, h.accession);
        assert.equal(proteins.decoy[p], decoy ? 1 : 0, h.accession);
        assert.equal(proteins.description[p], target.description ?? "", h.accession);
        if (decoy) {
          assert.equal(proteins.sequence[p].length, target.sequence.length, h.accession);
          assert.notEqual(proteins.sequence[p], target.sequence, h.accession);
        } else {
          assert.equal(proteins.sequence[p], target.sequence, h.accession);
        }
        assertSame(proteins.qValue[p], h.score, `protein q-value of ${h.accession}`);
        assert.ok(Number.isFinite(proteins.score[p]), `score of ${h.accession}`);
        assertSame(proteins.score[p], bestScore[p] ?? NaN, `score of ${h.accession}: the best of its PSMs`);
        assert.equal(proteins.group[p], groupOf.get(h.accession) ?? -1, `group of ${h.accession}`);
      });
      assert.ok(proteins.decoy.includes(1), "decoy proteins");
    } finally {
      set.delete();
    }
  });
});

/// Stores the spectra of TOPP_ProSE_1 and, as further spectra, the b and y ions of @p peptides ([OpenMS sequence,
/// precursor charge]) as mzML at @p file.
function storeSpectraWith(file, peptides) {
  using([new OpenMS.MSExperiment(), new OpenMS.MzMLFile(), new OpenMS.TheoreticalSpectrumGenerator()], (experiment, mzml, generator) => {
    mzml.load("/topp/SimpleSearchEngine_1.mzML", experiment);
    for (const [sequence, charge] of peptides) {
      using([OpenMS.AASequence.fromString(sequence), new OpenMS.MSSpectrum()], (peptide, spectrum) => {
        generator.getSpectrum(spectrum, peptide, 1, 2);
        const { mz } = spectrum.getPeaks();
        spectrum.setPeaks(mz, new Float64Array(mz.length).fill(100));
        spectrum.setMSLevel(2);
        spectrum.setRT(5000 + experiment.size());
        spectrum.setNativeID(`spectrum=${experiment.size()}`);
        spectrum.setPrecursors([{ mz: peptide.getMZ(charge), charge }]);
        experiment.addSpectrum(spectrum);
      });
    }
    mzml.store(file, experiment);
  });
}

test("tables() encodes terminal modifications", () => {
  // the data of TOPP_ProSE_1 (with an oxidized M), and spectra of a peptide with an N-terminal and one with a C-terminal
  // modification
  const acetylated = ".(Acetyl)EWATYTC(Carbamidomethyl)VLGFHVYVPVR";
  const amidated = "LVTLAESPR.(Amidated)";
  const oxidized = "IALSRPNVEVVALNDPFITNDYAAYM(Oxidation)FK";
  withSpectraCopies(OpenMS, [], () => {
    storeSpectraWith("/searchfiles/terminal.mzML", [[acetylated, 3], [amidated, 2]]);
    const set = searchFiles(OpenMS, ["/searchfiles/terminal.mzML"], undefined, {
      parameters: { ...prose1, "modifications:variable": ["Acetyl (N-term)", "Oxidation (M)", "Amidated (C-term)"] },
      database: "/topp/SimpleSearchEngine_1.fasta",
    });
    try {
      const t = set.tables();
      assertTableShapes(t);
      const { psms, peptides, proteins } = t;
      assert.equal(psms.count, 5);
      const peptide = (openms) => {
        const k = peptides.openms.indexOf(openms);
        assert.ok(k >= 0, `${openms} in ${JSON.stringify(peptides.openms)}`);
        return k;
      };
      const mods = (k) => {
        const m = [];
        for (let i = peptides.modOffsets[k]; i < peptides.modOffsets[k + 1]; ++i) {
          m.push([peptides.modPosition[i], peptides.modName[i], peptides.modUnimod[i], Number(peptides.modMass[i].toFixed(4))]);
        }
        return m;
      };
      const proteinsOf = (k) => {
        const p = [];
        for (let r = peptides.proteinOffsets[k]; r < peptides.proteinOffsets[k + 1]; ++r) {
          assert.equal(proteins.sequence[peptides.protein[r]].slice(peptides.start[r], peptides.end[r]), peptides.sequence[k]);
          p.push(proteins.accession[peptides.protein[r]]);
        }
        return p.sort();
      };

      // N-terminal: position -1, "[Acetyl]-" in ProForma; the fixed modification is written like any other
      let k = peptide(acetylated);
      assert.equal(peptides.sequence[k], "EWATYTCVLGFHVYVPVR");
      assert.equal(peptides.proforma[k], "[Acetyl]-EWATYTC[Carbamidomethyl]VLGFHVYVPVR");
      assert.deepEqual(mods(k), [[-1, "Acetyl", 1, 42.0106], [6, "Carbamidomethyl", 4, 57.0215]]);
      assert.deepEqual(proteinsOf(k), ["BSA2"]);
      // C-terminal: position = length of the sequence, "-[Amidated]"; a peptide of three proteins
      k = peptide(amidated);
      assert.equal(peptides.proforma[k], "LVTLAESPR-[Amidated]");
      assert.deepEqual(mods(k), [[9, "Amidated", 2, -0.984]]);
      assert.deepEqual(proteinsOf(k), ["BSA3", "BSA333", "test"]);
      // a modified residue
      k = peptide(oxidized);
      assert.equal(peptides.proforma[k], "IALSRPNVEVVALNDPFITNDYAAYM[Oxidation]FK");
      assert.deepEqual(mods(k), [[25, "Oxidation", 35, 15.9949]]);
      // unmodified
      k = peptide("DFASSGGYVLHLHR");
      assert.equal(peptides.proforma[k], "DFASSGGYVLHLHR");
      assert.deepEqual(mods(k), []);
      // the OpenMS strings are those of AASequence, which reads them back
      for (const s of peptides.openms) using([OpenMS.AASequence.fromString(s)], (a) => assert.equal(a.toString(), s));

      // Without decoys: ProSE's ln(hyperscore) as score (TOPP_ProSE_1), no q-values, PEPs or protein q-values
      const row = psms.peptide.indexOf(peptide(oxidized));
      assertClose(assert, psms.score[row], 45.369950862989675, 1e-6, "ln(hyperscore)");
      assert.equal(psms.charge[row], 3);
      const hyperscore = psms.featureNames.indexOf("ln(hyperscore)");
      for (let i = 0; i < psms.count; ++i) {
        assert.equal(psms.features[i * psms.featureNames.length + hyperscore], psms.score[i], "the score is the feature");
        assert.ok(Number.isNaN(psms.qValue[i]) && Number.isNaN(psms.pep[i]) && Number.isNaN(psms.deltaMass[i]));
        assert.equal(psms.decoy[i], 0);
      }
      assert.ok(proteins.qValue.every(Number.isNaN));
      assert.ok(proteins.decoy.every((d) => d === 0));
    } finally {
      set.delete();
    }
  });
});

test("tables(): the decoys of an unspecific search", () => {
  // ProSE reverses whole proteins for the decoys of an unspecific search (peptide:enzyme_specificity none), not their
  // peptides. With the spectra of TOPP_ProSE_1, a spectrum of a stretch of the decoy of BSA2 (its reversed sequence
  // RVPVYVHFGLVCTYTAWEKFM...), which only that decoy has.
  const decoyPeptide = "VHFGLVC(Carbamidomethyl)TYTAWEK";
  withSpectraCopies(OpenMS, [], () => {
    storeSpectraWith("/searchfiles/unspecific.mzML", [[decoyPeptide, 2]]);
    const set = searchFiles(OpenMS, ["/searchfiles/unspecific.mzML"], { rescore: false }, {
      parameters: { ...prose1, "decoys": "generate", "peptide:enzyme_specificity": "none", "FDR:PSM": 0, "FDR:protein": 0 },
      database: "/topp/SimpleSearchEngine_1.fasta",
    });
    try {
      const t = set.tables();
      assertTableShapes(t);
      const { psms, peptides, proteins } = t;
      const k = peptides.openms.indexOf(decoyPeptide);
      assert.ok(k >= 0, `${decoyPeptide} in ${JSON.stringify(peptides.openms)}`);
      assert.equal(peptides.decoy[k], 1);
      assert.ok(psms.decoy[psms.peptide.indexOf(k)] === 1);
      // the decoy proteins have the sequences ProSE searched: every peptide is at [start, end) of each of its proteins
      let decoyReferences = 0;
      for (let j = 0; j < peptides.count; ++j) {
        for (let r = peptides.proteinOffsets[j]; r < peptides.proteinOffsets[j + 1]; ++r) {
          const p = peptides.protein[r];
          assert.equal(proteins.sequence[p].slice(peptides.start[r], peptides.end[r]), peptides.sequence[j],
            `${peptides.openms[j]} in ${proteins.accession[p]} at ${peptides.start[r]}-${peptides.end[r]}`);
          decoyReferences += proteins.decoy[p];
        }
      }
      assert.ok(decoyReferences > 0, "references to decoy proteins");
      const decoy = proteins.accession.indexOf("DECOY_BSA2");
      assert.ok(decoy >= 0, JSON.stringify(proteins.accession));
      assert.equal(proteins.decoy[decoy], 1);
      assert.equal(proteins.sequence[decoy], [..."IALSRPNVEVVALNDPFITNDYAAYMFKEWATYTCVLGFHVYVPVR"].reverse().join(""));
    } finally {
      set.delete();
    }
  });
});

/// Checks the shape of the stored spectra @p s (ProSEResultSet.spectra()) of a result set with @p rows PSM rows: the
/// types and lengths of the arrays, the offsets, and that every row refers to a spectrum or is -1. Returns the number
/// of spectra.
function assertSpectraShapes(s, rows) {
  const types = { psmSpectrum: Int32Array, offsets: Uint32Array, mz: Float32Array, intensity: Float32Array,
                  precursorMz: Float64Array, collisionEnergy: Float32Array };
  assert.deepEqual(Object.keys(s).sort(), Object.keys(types).sort());
  for (const [name, type] of Object.entries(types)) assert.ok(s[name] instanceof type, `${name} is a ${type.name}`);
  const n = s.precursorMz.length;
  assert.equal(s.psmSpectrum.length, rows, "psmSpectrum: one entry per PSM row");
  assert.equal(s.collisionEnergy.length, n, "collisionEnergy: one entry per spectrum");
  assert.equal(s.offsets.length, n + 1, "offsets: n + 1 entries");
  assert.equal(s.offsets[0], 0);
  for (let k = 0; k < n; ++k) assert.ok(s.offsets[k] <= s.offsets[k + 1], "offsets");
  assert.equal(s.mz.length, s.offsets[n], "mz: offsets[n] entries");
  assert.equal(s.intensity.length, s.offsets[n], "intensity: offsets[n] entries");
  assert.ok(s.psmSpectrum.every((k) => k >= -1 && k < n), "psmSpectrum: -1 or a spectrum");
  return n;
}

/// [m/z, intensity] pairs (as Float32) sorted by m/z, then intensity
function peakPairs(mz, intensity) {
  return Array.from(mz, (m, j) => [Math.fround(m), Math.fround(intensity[j])]).sort((a, b) => a[0] - b[0] || a[1] - b[1]);
}

/// The peaks of spectrum @p k of the stored spectra @p s, which are sorted by m/z, as peakPairs()
function storedPeaks(s, k) {
  const mz = s.mz.subarray(s.offsets[k], s.offsets[k + 1]);
  for (let j = 1; j < mz.length; ++j) assert.ok(mz[j - 1] <= mz[j], `the peaks of spectrum ${k} are sorted by m/z`);
  return peakPairs(mz, s.intensity.subarray(s.offsets[k], s.offsets[k + 1]));
}

/// The spectra of mzML file @p file, by native ID: their MS level, peaks (peakPairs()), whether the file has them sorted
/// by m/z, and precursors
function fileSpectra(file) {
  return using([new OpenMS.MSExperiment(), new OpenMS.MzMLFile()], (experiment, mzml) => {
    mzml.setOptions({ sortSpectraByMZ: false }); // the peaks in the order of the file
    mzml.load(file, experiment);
    const spectra = new Map();
    for (let i = 0; i < experiment.size(); ++i) {
      using([experiment.getSpectrum(i)], (spectrum) => {
        const { mz, intensity } = spectrum.getPeaks();
        spectra.set(spectrum.getNativeID(),
          { msLevel: spectrum.getMSLevel(), peaks: peakPairs(mz, intensity), sorted: spectrum.isSorted(),
            precursors: spectrum.getPrecursors() });
      });
    }
    return spectra;
  });
}

/// Stores the spectra of mzML file @p source at @p file, with an activation energy of @p energy on every precursor (a
/// marker of the file)
function storeWithActivationEnergy(source, file, energy) {
  using([new OpenMS.MSExperiment(), new OpenMS.MzMLFile()], (experiment, mzml) => {
    mzml.setOptions({ sortSpectraByMZ: false });
    mzml.load(source, experiment);
    for (let i = 0; i < experiment.size(); ++i) {
      using([experiment.getSpectrum(i)], (spectrum) => {
        spectrum.setPrecursors(spectrum.getPrecursors().map((p) => ({ ...p, activationEnergy: energy })));
        experiment.setSpectrum(i, spectrum);
      });
    }
    mzml.store(file, experiment);
  });
}

test("spectra() holds the peaks of the identified spectra", () => {
  // two files with the same name and spectra, the second with an activation energy of 33 eV on every precursor: PSMs
  // of both with the same native IDs, whose spectra are those of their own file
  withSpectraCopies(OpenMS, ["0/run.mzML"], ([first], dir) => {
    fs.mkdirSync(path.join(dir, "1"));
    const files = [first, "/searchfiles/1/run.mzML"];
    storeWithActivationEnergy(first, files[1], 33);
    const set = searchFiles(OpenMS, files, { rescore: false });
    try {
      const { psms } = set.tables();
      const s = set.spectra();
      const n = assertSpectraShapes(s, psms.count);
      const { psmFdr } = set.summary();
      assert.equal(psmFdr, 0.05);

      // every target PSM at 1 % FDR has its spectrum, with peaks
      let targets = 0;
      for (let i = 0; i < psms.count; ++i) {
        if (psms.decoy[i] === 1 || !(psms.qValue[i] <= 0.01)) continue;
        const k = s.psmSpectrum[i];
        assert.ok(k >= 0, `the spectrum of PSM ${i} (${psms.nativeId[i]})`);
        assert.ok(s.offsets[k + 1] > s.offsets[k], `the peaks of PSM ${i}`);
        ++targets;
      }
      assert.ok(targets > 100, `${targets} target PSMs at 1 % FDR`);
      // every PSM that passes FDR:PSM, decoys too (here all: the search filtered them), has its spectrum: the file has
      // every spectrum ProSE identified
      assert.ok(psms.decoy.includes(1), "decoy PSMs");
      assert.ok(psms.qValue.every((q) => q <= psmFdr));
      assert.ok(s.psmSpectrum.every((k) => k >= 0), "every PSM has its spectrum");

      // the spectra of the PSMs' own files, as read: the peaks sorted by m/z, the precursor of the PSM, and the
      // collision energy of the file (the first records none, the second 33 eV)
      assert.deepEqual([...new Set(psms.file)].sort(), [0, 1], "PSMs of both files");
      const expected = files.map(fileSpectra);
      const keyOf = new Map(); // spectrum -> its file and native ID
      for (let i = 0; i < psms.count; ++i) {
        const k = s.psmSpectrum[i];
        const key = `${psms.file[i]} ${psms.nativeId[i]}`;
        assert.ok(!keyOf.has(k) || keyOf.get(k) === key, `spectrum ${k}: ${keyOf.get(k)} and ${key}`);
        keyOf.set(k, key);
        const e = expected[psms.file[i]].get(psms.nativeId[i]);
        assert.ok(e && e.msLevel === 2, key);
        assert.deepEqual(storedPeaks(s, k), e.peaks, `the peaks of ${key}`);
        assertClose(assert, s.precursorMz[k], psms.mz[i], 1e-4, `precursor m/z of ${key}`);
        assert.equal(s.precursorMz[k], e.precursors[0].mz, `precursor m/z of ${key}`);
        if (psms.file[i] === 1) {
          assert.equal(s.collisionEnergy[k], 33, `the collision energy of ${key}`);
        } else {
          assert.ok(Number.isNaN(s.collisionEnergy[k]), `no collision energy for ${key}: ${s.collisionEnergy[k]}`);
        }
      }
      assert.equal(keyOf.size, n, "every spectrum is that of a PSM");
      assert.equal(new Set(keyOf.values()).size, n, "one spectrum per file and native ID");
      // the same native ID in both files: two spectra, with the same peaks, each with the energy of its file
      const inSecond = new Map();
      for (let i = 0; i < psms.count; ++i) if (psms.file[i] === 1) inSecond.set(psms.nativeId[i], s.psmSpectrum[i]);
      let both = 0;
      for (let i = 0; i < psms.count; ++i) {
        if (psms.file[i] !== 0 || !inSecond.has(psms.nativeId[i])) continue;
        const [k0, k1] = [s.psmSpectrum[i], inSecond.get(psms.nativeId[i])];
        assert.notEqual(k0, k1, psms.nativeId[i]);
        assert.deepEqual(storedPeaks(s, k0), storedPeaks(s, k1), psms.nativeId[i]);
        assert.ok(Number.isNaN(s.collisionEnergy[k0]) && s.collisionEnergy[k1] === 33, psms.nativeId[i]);
        ++both;
      }
      assert.ok(both > 100, `${both} native IDs with PSMs in both files`);
    } finally {
      set.delete();
    }
  });
});

/// Stores the spectra of TOPP_ProSE_1 as mzML at @p file, changed for the tests of spectra(): an MS1 spectrum ("scan=1")
/// first, the first MS2 spectrum with a collision energy of 27.5 eV, and the peaks of the second in reverse order. (The
/// native IDs have a '=': otherwise MzMLFile writes "spectrum=<index>" for all.)
function storeSpectraToKeep(file) {
  using([new OpenMS.MSExperiment(), new OpenMS.MSExperiment(), new OpenMS.MzMLFile(), new OpenMS.MSSpectrum()],
    (original, changed, mzml, ms1) => {
      mzml.load("/topp/SimpleSearchEngine_1.mzML", original);
      ms1.setMSLevel(1);
      ms1.setNativeID("scan=1");
      ms1.setPeaks([400, 500], [10, 20]);
      changed.addSpectrum(ms1);
      for (let i = 0; i < original.size(); ++i) {
        using([original.getSpectrum(i)], (spectrum) => {
          if (i === 0) {
            spectrum.setPrecursors(spectrum.getPrecursors().map((p) => ({ mz: p.mz, charge: p.charge, activationEnergy: 27.5 })));
          }
          if (i === 1) {
            const { mz, intensity } = spectrum.getPeaks();
            spectrum.setPeaks(Array.from(mz).reverse(), Array.from(intensity).reverse());
          }
          changed.addSpectrum(spectrum);
        });
      }
      mzml.store(file, changed);
    });
}

/// ProSEAlgorithm.searchFiles() of the data of TOPP_ProSE_1 changed by storeSpectraToKeep() (three PSMs, no q-values)
function searchSpectraToKeep(options, onProgress) {
  storeSpectraToKeep("/searchfiles/keep.mzML");
  return searchFiles(OpenMS, ["/searchfiles/keep.mzML"], options,
    { parameters: prose1, database: "/topp/SimpleSearchEngine_1.fasta", onProgress });
}

test("spectra() marks PSMs without a readable spectrum with -1", () => {
  withSpectraCopies(OpenMS, [], () => {
    const set = searchSpectraToKeep(undefined);
    try {
      const { nativeId } = set.tables().psms;
      assert.equal(nativeId.length, 3);
      const spectra = () => {
        const s = set.spectra();
        return [s.psmSpectrum, assertSpectraShapes(s, 3)];
      };
      let [psmSpectrum, n] = spectra();
      assert.equal(n, 3);
      assert.ok(psmSpectrum.every((k) => k >= 0));

      // a native ID that the file does not have (the test hook re-reads the file)
      set._debugSetNativeId(0, "scan=99999999");
      assert.equal(set.tables().psms.nativeId[0], "scan=99999999");
      [psmSpectrum, n] = spectra();
      assert.equal(psmSpectrum[0], -1);
      assert.ok(psmSpectrum[1] >= 0 && psmSpectrum[2] >= 0 && psmSpectrum[1] !== psmSpectrum[2]);
      assert.equal(n, 2);
      // that of an MS1 spectrum: only MS2 spectra are kept
      assert.equal(fileSpectra("/searchfiles/keep.mzML").get("scan=1").msLevel, 1);
      set._debugSetNativeId(0, "scan=1");
      [psmSpectrum, n] = spectra();
      assert.equal(psmSpectrum[0], -1);
      assert.equal(n, 2);
      // two PSMs of one spectrum share it
      set._debugSetNativeId(0, nativeId[1]);
      [psmSpectrum, n] = spectra();
      assert.ok(psmSpectrum[0] >= 0);
      assert.equal(psmSpectrum[0], psmSpectrum[1]);
      assert.equal(n, 2);

      assert.throws(() => set._debugSetNativeId(3, "x"), /Index/);
      assert.throws(() => set._debugSetNativeId(-1, "x"), /Index/);
    } finally {
      set.delete();
    }
  });
});

test("spectra(): the peaks as read, sorted by m/z, and the collision energy of the precursor", () => {
  withSpectraCopies(OpenMS, [], () => {
    const set = searchSpectraToKeep(undefined);
    try {
      const { psms } = set.tables();
      const s = set.spectra();
      assert.equal(assertSpectraShapes(s, psms.count), 3);
      // no q-values (decoys: ignore): the spectra of all PSMs
      assert.ok(psms.qValue.every(Number.isNaN));
      assert.ok(s.psmSpectrum.every((k) => k >= 0));
      const expected = fileSpectra("/searchfiles/keep.mzML");
      assert.equal(expected.get("spectrum=0").precursors[0].activationEnergy, 27.5, "the file records the energy");
      assert.equal(expected.get("spectrum=1").sorted, false, "the file has the peaks of spectrum=1 in reverse order");
      for (let i = 0; i < psms.count; ++i) {
        const k = s.psmSpectrum[i];
        const e = expected.get(psms.nativeId[i]);
        // the second spectrum of TOPP_ProSE_1 has its peaks in reverse order in the file
        assert.deepEqual(storedPeaks(s, k), e.peaks, psms.nativeId[i]);
        assert.equal(s.precursorMz[k], e.precursors[0].mz);
        if (psms.nativeId[i] === "spectrum=0") {
          assert.equal(s.collisionEnergy[k], 27.5);
        } else {
          assert.ok(Number.isNaN(s.collisionEnergy[k]), `no energy (0 in OpenMS) for ${psms.nativeId[i]}`);
        }
      }
      assert.deepEqual(new Set(psms.nativeId), new Set(["spectrum=0", "spectrum=1", "spectrum=2"]));
    } finally {
      set.delete();
    }
  });
});

test("spectra() without keepSpectra: no peaks", () => {
  withSpectraCopies(OpenMS, [], () => {
    const reports = [];
    const set = searchSpectraToKeep({ keepSpectra: false }, (r) => reports.push(r));
    try {
      // no step that reads the spectra
      assertProgressSteps(reports, ["Searching keep.mzML (1 of 1)..."]);
      const s = set.spectra();
      assert.equal(assertSpectraShapes(s, 3), 0);
      assert.deepEqual(Array.from(s.psmSpectrum), [-1, -1, -1]);
      assert.deepEqual(Array.from(s.offsets), [0]);
    } finally {
      set.delete();
    }
  });
});

/// Stores the spectra of TOPP_ProSE_1 as mzML at @p file, and two more under native IDs it has: spectrum=1 again at
/// RT + 100 s, with an activation energy of 40 eV, and spectrum=2 again at the same RT.
function storeRepeatedNativeIds(file) {
  using([new OpenMS.MSExperiment(), new OpenMS.MzMLFile()], (experiment, mzml) => {
    mzml.load("/topp/SimpleSearchEngine_1.mzML", experiment);
    using([experiment.getSpectrum(1), experiment.getSpectrum(2)], (one, two) => {
      one.setRT(one.getRT() + 100);
      one.setPrecursors(one.getPrecursors().map((p) => ({ ...p, activationEnergy: 40 })));
      experiment.addSpectrum(one);
      experiment.addSpectrum(two);
    });
    mzml.store(file, experiment);
  });
}

test("spectra(): a native ID that a file has more than once", () => {
  withSpectraCopies(OpenMS, [], () => {
    storeRepeatedNativeIds("/searchfiles/repeated.mzML");
    const set = searchFiles(OpenMS, ["/searchfiles/repeated.mzML"], undefined,
      { parameters: prose1, database: "/topp/SimpleSearchEngine_1.fasta" });
    try {
      const { psms } = set.tables();
      let s = set.spectra();
      const rowsOf = (nativeId) => psms.nativeId.flatMap((id, i) => (id === nativeId ? [i] : []));
      // spectrum=1 twice, at different RTs: the PSM of each has the spectrum at its RT
      const ones = rowsOf("spectrum=1").sort((a, b) => psms.rt[a] - psms.rt[b]);
      assert.equal(ones.length, 2, JSON.stringify(psms.nativeId));
      const [early, late] = ones.map((i) => s.psmSpectrum[i]);
      assert.ok(early >= 0 && late >= 0 && early !== late, `spectra ${early}, ${late}`);
      assert.ok(Number.isNaN(s.collisionEnergy[early]));
      assert.equal(s.collisionEnergy[late], 40);
      // spectrum=2 twice, at the same RT: neither is the one at the RT of its PSMs
      const twos = rowsOf("spectrum=2");
      assert.ok(twos.length > 0);
      assert.ok(twos.every((i) => s.psmSpectrum[i] === -1), `spectrum=2: ${twos.map((i) => s.psmSpectrum[i])}`);
      // spectrum=0 once
      const [zero] = rowsOf("spectrum=0");
      assert.ok(s.psmSpectrum[zero] >= 0);
      // only the spectra of PSMs are kept: that of spectrum=0, and the two of spectrum=1
      assert.equal(assertSpectraShapes(s, psms.count), 3);
      assert.equal(new Set(Array.from(s.psmSpectrum).filter((k) => k >= 0)).size, 3);

      // a PSM at an RT of neither spectrum of a native ID that the file has twice: none
      set._debugSetNativeId(zero, "spectrum=1");
      s = set.spectra();
      assert.equal(s.psmSpectrum[zero], -1);
      assert.equal(assertSpectraShapes(s, psms.count), 2);
    } finally {
      set.delete();
    }
  });
});

test("spectra(): the PSMs of a file that cannot be read again have none", () => {
  withSpectraCopies(OpenMS, [], (_, dir) => {
    const files = ["0/run.mzML", "1/run.mzML"];
    for (const f of files) {
      fs.mkdirSync(path.dirname(path.join(dir, f)));
      storeSpectraToKeep(`/searchfiles/${f}`);
    }
    const set = searchFiles(OpenMS, files.map((f) => `/searchfiles/${f}`), undefined,
      { parameters: prose1, database: "/topp/SimpleSearchEngine_1.fasta" });
    try {
      const { psms } = set.tables();
      const before = set.spectra();
      assert.deepEqual([...new Set(psms.file)].sort(), [0, 1]);
      assert.ok(before.psmSpectrum.every((k) => k >= 0));
      // the second file is gone: reading the files again (through the test hook) does not throw
      fs.unlinkSync(path.join(dir, files[1]));
      set._debugSetNativeId(0, psms.nativeId[0]);
      const after = set.spectra();
      assertSpectraShapes(after, psms.count);
      for (let i = 0; i < psms.count; ++i) {
        const [k0, k1] = [before.psmSpectrum[i], after.psmSpectrum[i]];
        if (psms.file[i] === 1) {
          assert.equal(k1, -1, `PSM ${i} of the file that is gone`);
        } else {
          assert.ok(k1 >= 0, `PSM ${i}`);
          assert.deepEqual(storedPeaks(after, k1), storedPeaks(before, k0), `the peaks of PSM ${i}`);
          assert.equal(after.precursorMz[k1], before.precursorMz[k0]);
          assert.ok(Object.is(after.collisionEnergy[k1], before.collisionEnergy[k0]));
        }
      }
    } finally {
      set.delete();
    }
  });
});

test("spectra() with FDR:PSM 0 (no filter): the spectra of all PSMs, whatever their q-values", () => {
  // Percolator gives every PSM a q-value, which FDR:PSM 0 does not filter: these PSMs pass, and keep their spectra
  withSpectraCopies(OpenMS, ["a.mzML"], (files) => {
    const set = searchFiles(OpenMS, files, { rescore: true }, { parameters: { ...prose6, "FDR:PSM": 0, "FDR:protein": 0 } });
    try {
      assert.equal(set.summary().files[0].percolator.applied, true);
      assert.equal(set.summary().psmFdr, 0);
      const { psms } = set.tables();
      assert.ok(psms.qValue.every(Number.isFinite), "q-values");
      assert.ok(psms.qValue.some((q) => q > 0), "q-values above FDR:PSM");
      const s = set.spectra();
      assert.equal(assertSpectraShapes(s, psms.count), psms.count);
      assert.ok(s.psmSpectrum.every((k) => k >= 0), "every PSM has its spectrum");
    } finally {
      set.delete();
    }
  });
});

test("ProSEResultSet.delete() releases the results", async () => {
  // A module of its own, whose memory grows with these searches only. WebAssembly memory never shrinks, and it grows
  // in steps of 20 % (Emscripten's MEMORY_GROWTH_GEOMETRIC_STEP): result sets that stay allocated show as a step after
  // a few searches (three, without delete()). One thread: the search allocates the same way each time. (With several
  // threads, the heap of a threaded search grows from search to search also without result sets, as that of search().)
  const { default: createOpenMS } = await import(pathToFileURL(modulePath).href);
  const m = await createOpenMS({ print: () => {}, printErr: () => {}, numThreads: 1 });
  m.mountHostDirectory("/topp", toppDataDir);
  withSpectraCopies(m, ["a.mzML", "b.mzML"], (files) => {
    const heapAfterSearch = () => {
      const set = searchFiles(m, files, { rescore: false });
      assert.equal(set.summary().files.length, 2);
      set.delete();
      assert.throws(() => set.summary(), /deleted/);
      return m.getMemoryInfo().heapSize;
    };
    const h1 = heapAfterSearch();
    for (let i = 2; i <= 4; ++i) {
      const h = heapAfterSearch();
      assert.ok(h - h1 < 0.2 * h1, `the heap grew from ${h1} to ${h} bytes by search ${i}`);
    }
  });
});

// WebAssembly: ProSE's PeptDeep batches are smaller than natively (500): ONNX Runtime has no memory arena there, and the
// memory of large tensors is not reused.
test("ProSEAlgorithm: peptdeep:batch_size is 128 by default", () => {
  using([new OpenMS.ProSEAlgorithm()], (prose) => {
    using([prose.getDefaults()], (defaults) => assert.equal(defaults.toObject()["peptdeep:batch_size"], 128));
    assert.equal(prose.getParameterValues()["peptdeep:batch_size"], 128);
  });
});

// A module with ONNX Runtime: a PeptDeep-rescored search does not take memory per file beyond the search's own. A module
// of its own, whose memory grows with these searches only. The spectra of DatabaseSuitability_in_spec.mzML twice (the
// copy with other native IDs and later retention times): the sample of confident PSMs on which PeptDeepRescoring selects
// the collision energy then holds most peptides twice. Before ONNXPredictorBase did without ONNX Runtime's memory
// pattern and PeptDeepRescoring predicted each peptide of that sample once, the heap was 1.4 GB after the first search
// and 3.3 GB after the second. The heap of a threaded search still grows from search to search, also without PeptDeep
// (see "ProSEResultSet.delete() releases the results"), but by far less than the bound.
test("ProSEAlgorithm: repeated PeptDeep searches do not grow the heap per file (ONNX Runtime)", async (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  const { default: createOpenMS } = await import(pathToFileURL(modulePath).href);
  const m = await createOpenMS({ print: () => {}, printErr: () => {}, numThreads: 4 });
  m.mountHostDirectory("/topp", toppDataDir);
  m.mountHostDirectory("/models", path.join(path.dirname(modulePath), "models"));
  m.FS.mkdir("/twice");
  using([new m.MzMLFile(), new m.MSExperiment(), new m.MSExperiment()], (mzml, spectra, twice) => {
    mzml.load("/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML", spectra);
    let lastRT = 0;
    for (let i = 0; i < spectra.size(); ++i) using([spectra.getSpectrum(i)], (s) => { lastRT = Math.max(lastRT, s.getRT()); });
    for (let copy = 0; copy < 2; ++copy) {
      for (let i = 0; i < spectra.size(); ++i) {
        using([spectra.getSpectrum(i)], (s) => {
          s.setNativeID(`scan=${copy * spectra.size() + i + 1}`);
          s.setRT(s.getRT() + copy * (lastRT + 1));
          twice.addSpectrum(s);
        });
      }
    }
    mzml.store("/twice/spectra.mzML", twice);
  });
  const heapAfterSearch = () => using([new m.ProSEAlgorithm(), new m.ProteinIdentificationList(), new m.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues({
        "enzyme": "Trypsin", "decoys": "generate",
        "precursor:mass_tolerance_lower": 20, "precursor:mass_tolerance_upper": 20, "precursor:mass_tolerance_unit": "ppm",
        "fragment:mass_tolerance": 10, "fragment:mass_tolerance_unit": "ppm",
        "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
        "FDR:PSM": 0.01, "FDR:protein": 0,
        "peptdeep:enable": true,
        "peptdeep:ms2_model": "/models/peptdeep_ms2_dynamic.onnx", "peptdeep:rt_model": "/models/peptdeep_rt_dynamic.onnx",
      });
      const summary = prose.search("/twice/spectra.mzML", suitabilityDatabase, proteins, peptides, { rescore: false });
      assert.equal(summary.status, "EXECUTION_OK");
      return m.getMemoryInfo().heapSize;
    });
  const h1 = heapAfterSearch();
  const h2 = heapAfterSearch();
  assert.ok(h2 - h1 < 400 * 1024 * 1024, `the heap grew from ${h1} to ${h2} bytes by the second search`);
});

test("annotateSpectrum: b and y ions of a theoretical spectrum", () => {
  using([OpenMS.AASequence.fromString("PEPTIDEK"), new OpenMS.TheoreticalSpectrumGenerator(), new OpenMS.MSSpectrum()],
    (peptide, generator, spectrum) => {
      generator.getSpectrum(spectrum, peptide, 1, 1);
      const { mz, intensity } = spectrum.getPeaks();
      // reverse the peaks: the measured spectrum need not be sorted
      spectrum.setPeaks(Array.from(mz).reverse(), Array.from(intensity).reverse());

      const a = OpenMS.annotateSpectrum(spectrum, peptide, 2, { tolerance: 0.01 });
      assert.equal(a.unit, "Da");
      assert.equal(a.matches.length, mz.length, "every peak is a b or y ion");
      const ions = a.matches.map((m) => m.ion);
      for (const ion of ["b2+", "y1+", "y7+"]) assert.ok(ions.includes(ion), `${ion} in ${ions}`);
      for (const m of a.matches) {
        assertClose(assert, m.error, 0, 1e-6, m.ion);
        assert.equal(m.charge, 1);
      }
      assert.equal(a.theoretical.mz.length, a.theoretical.ion.length);

      // ppm, and only y ions
      const y = OpenMS.annotateSpectrum(spectrum, peptide, 2, { tolerance: 5, unit: "ppm", ions: ["y"] });
      assert.equal(y.unit, "ppm");
      assert.ok(y.matches.length > 0 && y.matches.every((m) => m.ion.startsWith("y")));

      assert.throws(() => OpenMS.annotateSpectrum(spectrum, peptide, 2, { unit: "mmu" }), /unit/);
      assert.throws(() => OpenMS.annotateSpectrum(spectrum, peptide, 2, { ions: ["q"] }), /ion type/);
    });
});

/// The stored peaks of PSM row @p row of the stored spectra @p s (ProSEResultSet.spectra()): views of its arrays
function psmPeaks(s, row) {
  const k = s.psmSpectrum[row];
  assert.ok(k >= 0, `PSM ${row} has a spectrum`);
  return { mz: s.mz.subarray(s.offsets[k], s.offsets[k + 1]), intensity: s.intensity.subarray(s.offsets[k], s.offsets[k + 1]) };
}

/// The ion label for the fields of a PeakMatch or theoretical ion: 'b3' for charge 1, 'y5++' for charge 2, 'y7-H2O'
/// with a loss
function ionLabel({ type, ordinal, charge, loss }) {
  return `${type}${ordinal}${loss ? "-" + loss : ""}${charge > 1 ? "+".repeat(charge) : ""}`;
}

/// The PSMs at 1 % FDR of the search of DatabaseSuitability_in_spec.mzML (without Percolator), every @p stride-th, with
/// the tables and stored spectra of the result set; @p fn gets them.
function withPsmSpectra(stride, fn) {
  withSpectraCopies(OpenMS, ["run.mzML"], (files) => {
    const set = searchFiles(OpenMS, files, { rescore: false });
    try {
      const { psms, peptides } = set.tables();
      const s = set.spectra();
      const rows = [];
      for (let i = 0; i < psms.count; i += stride) {
        if (psms.decoy[i] === 0 && psms.qValue[i] <= 0.01 && s.psmSpectrum[i] >= 0) rows.push(i);
      }
      fn({ psms, peptides, s, rows });
    } finally {
      set.delete();
    }
  });
}

const ppm10 = { tolerance: 10, unit: "ppm" };

test("annotatePeaks matches b and y ions of a stored spectrum", () => {
  withPsmSpectra(5, ({ psms, peptides, s, rows }) => {
    assert.ok(rows.length > 30, `${rows.length} PSMs`);
    const types = new Set();
    let matched = 0;
    for (const row of rows) {
      const peptide = peptides.openms[psms.peptide[row]];
      const n = peptides.sequence[psms.peptide[row]].length;
      const { mz, intensity } = psmPeaks(s, row);
      const r = OpenMS.annotatePeaks(mz, intensity, peptide, psms.charge[row], ppm10);
      const what = `PSM ${row} ${peptide}/${psms.charge[row]}`;
      assert.deepEqual(Object.keys(r).sort(), ["matches", "theoretical"]);

      // the theoretical ions: b1..b(n-1) and y1..y(n-1) at charge 1, unique, sorted by m/z
      const labels = new Set(r.theoretical.map((t) => t.ion));
      assert.equal(labels.size, r.theoretical.length, `${what}: unique ions`);
      for (let i = 1; i < n; ++i) {
        for (const type of ["b", "y"]) {
          const t = r.theoretical.find((t) => t.ion === `${type}${i}`);
          assert.ok(t, `${what}: ${type}${i} in theoretical`);
          assert.deepEqual(Object.keys(t).sort(), ["charge", "ion", "mz", "ordinal", "type"]);
          assert.deepEqual([t.type, t.ordinal, t.charge], [type, i, 1]);
        }
      }
      for (const t of r.theoretical) {
        assert.equal(ionLabel({ ...t, loss: "" }), t.ion, what);
        assert.ok(t.mz > 0 && t.charge >= 1 && t.charge <= Math.max(1, psms.charge[row] - 1), `${what}: ${t.ion}`);
      }
      for (let i = 1; i < r.theoretical.length; ++i) assert.ok(r.theoretical[i - 1].mz <= r.theoretical[i].mz);

      // the matches
      assert.ok(r.matches.length > 0, `${what}: matches`);
      for (const m of r.matches) {
        assert.deepEqual(Object.keys(m).sort(), ["charge", "errorPpm", "ion", "loss", "ordinal", "peak", "theoreticalMz", "type"]);
        assert.ok(Number.isInteger(m.peak) && m.peak >= 0 && m.peak < mz.length, `${what}: peak ${m.peak}`);
        assert.ok(m.type === "b" || m.type === "y", `${what}: ${m.ion}`);
        types.add(m.type);
        assert.equal(m.loss, "");
        assert.equal(ionLabel(m), m.ion);
        assert.ok(Number.isInteger(m.ordinal) && m.ordinal >= 1 && m.ordinal < n, `${what}: ${m.ion}`);
        assert.ok(r.theoretical.some((t) => t.ion === m.ion && t.mz === m.theoreticalMz && t.charge === m.charge), `${what}: ${m.ion}`);
        // in ppm, from the peak as passed: (observed - theoretical) / theoretical
        assertClose(assert, m.errorPpm, (mz[m.peak] - m.theoreticalMz) / m.theoreticalMz * 1e6, 1e-6, `${what}: ${m.ion}`);
        assert.ok(Math.abs(m.errorPpm) <= 10, `${what}: ${m.ion} is ${m.errorPpm} ppm off`);
      }
      matched += r.matches.length;
    }
    assert.deepEqual([...types].sort(), ["b", "y"]);
    assert.ok(matched > 5 * rows.length, `${matched} matches in ${rows.length} PSMs`);

    // the spectra of this file are made of b and y ions: the first PSM's peaks 5 ppm above and below theirs are matched
    // within 10 ppm (errorPpm 5 and -5, the m/z being Float32) but not within 4
    const row = rows[0];
    const peptide = peptides.openms[psms.peptide[row]];
    const { mz, intensity } = psmPeaks(s, row);
    for (const shift of [5, -5]) {
      const moved = Float32Array.from(mz, (x) => x * (1 + shift * 1e-6));
      const r = OpenMS.annotatePeaks(moved, intensity, peptide, psms.charge[row], ppm10);
      assert.equal(r.matches.length, mz.length, `${shift} ppm`);
      for (const m of r.matches) assertClose(assert, m.errorPpm, shift, 0.2, `${m.ion} ${shift} ppm`);
      assert.deepEqual(OpenMS.annotatePeaks(moved, intensity, peptide, psms.charge[row], { tolerance: 4, unit: "ppm" }).matches, []);
      // and in Da, ppm all the same: 0.01 Da is more than 5 ppm for ions below 2000
      const da = OpenMS.annotatePeaks(moved, intensity, peptide, psms.charge[row], { tolerance: 0.01 });
      assert.equal(da.matches.length, mz.length, `${shift} ppm in Da`);
      for (const m of da.matches) assertClose(assert, m.errorPpm, shift, 0.2, `${m.ion} ${shift} ppm in Da`);
    }
  });
});

test("annotatePeaks: peak is the index into the arrays as passed, whatever their order", () => {
  withPsmSpectra(25, ({ psms, peptides, s, rows }) => {
    assert.ok(rows.length >= 5, `${rows.length} PSMs`);
    let seed = 20261016; // a fixed shuffle
    const random = () => (seed = (seed * 1103515245 + 12345) % 2147483648) / 2147483648;
    // the matches as [m/z of the peak, its intensity, ion, theoretical m/z, error], found by the index
    const matches = (r, mz, intensity) => r.matches
      .map((m) => [mz[m.peak], intensity[m.peak], m.ion, m.theoreticalMz, m.errorPpm])
      .sort((a, b) => a[0] - b[0] || (a[2] < b[2] ? -1 : 1));
    let found = 0;
    for (const row of rows) {
      const peptide = peptides.openms[psms.peptide[row]];
      const charge = psms.charge[row];
      const { mz, intensity } = psmPeaks(s, row);
      // no peaks of equal m/z: which of them an ion takes depends on their order (see the test of equal m/z below)
      assert.equal(new Set(mz).size, mz.length, `${peptide}: peaks of different m/z`);
      const sorted = OpenMS.annotatePeaks(mz, intensity, peptide, charge, ppm10);
      assert.ok(sorted.matches.length > 0);
      for (const m of sorted.matches) {
        assert.ok(Math.abs(mz[m.peak] - m.theoreticalMz) / m.theoreticalMz * 1e6 <= 10, `${peptide}: ${m.ion} is peak ${m.peak}`);
      }

      const order = Array.from(mz, (_, i) => i); // order[j]: the position of the j-th peak of the copy in the original
      for (let i = order.length - 1; i > 0; --i) {
        const j = Math.floor(random() * (i + 1));
        [order[i], order[j]] = [order[j], order[i]];
      }
      assert.ok(order.some((o, j) => o !== j), "shuffled");
      const mz2 = Float32Array.from(order, (o) => mz[o]);
      const intensity2 = Float32Array.from(order, (o) => intensity[o]);
      const shuffled = OpenMS.annotatePeaks(mz2, intensity2, peptide, charge, ppm10);
      for (const m of shuffled.matches) {
        assert.ok(Math.abs(mz2[m.peak] - m.theoreticalMz) / m.theoreticalMz * 1e6 <= 10, `${peptide}: ${m.ion} is peak ${m.peak}`);
      }
      // the same ions on the same peaks, with the same errors, but at the indices of the copy
      assert.deepEqual(matches(shuffled, mz2, intensity2), matches(sorted, mz, intensity), peptide);
      assert.deepEqual(shuffled.matches.map((m) => order[m.peak]).sort((a, b) => a - b),
                       sorted.matches.map((m) => m.peak).sort((a, b) => a - b), peptide);
      assert.deepEqual(shuffled.theoretical, sorted.theoretical);

      // the same matching as annotateSpectrum() of a spectrum with these peaks, which names the ions as TheoreticalSpectrumGenerator
      const spectrumMatches = using([OpenMS.AASequence.fromString(peptide), new OpenMS.MSSpectrum()], (aa, spectrum) => {
        spectrum.setPeaks(mz2, intensity2);
        return OpenMS.annotateSpectrum(spectrum, aa, charge, ppm10).matches;
      });
      assert.deepEqual(spectrumMatches.map((m) => [m.mz, m.intensity, m.ion, m.charge, m.theoreticalMZ, m.error]),
                       shuffled.matches.map((m) => [mz2[m.peak], intensity2[m.peak], m.charge > 1 ? m.ion : m.ion + "+", m.charge,
                                                    m.theoreticalMz, m.errorPpm]), `${peptide}: as annotateSpectrum`);

      // JavaScript arrays and the reverse order
      const plain = OpenMS.annotatePeaks(Array.from(mz2), Array.from(intensity2), peptide, charge, ppm10);
      assert.deepEqual(plain, shuffled);
      const reversed = OpenMS.annotatePeaks(Float32Array.from(mz).reverse(), Float32Array.from(intensity).reverse(), peptide, charge, ppm10);
      assert.deepEqual(reversed.matches.map((m) => mz.length - 1 - m.peak).sort((a, b) => a - b),
                       sorted.matches.map((m) => m.peak).sort((a, b) => a - b), peptide);
      found += sorted.matches.length;
    }
    assert.ok(found > 5 * rows.length);
  });
});

/// The theoretical ions of PEPTIDEK (precursor charge 3) for all ion types with losses, which include z1++ and y1-NH3++ of the
/// same m/z, and that m/z. (Called inside the tests: with annotatePeaks() called while this file loads, node aborted at
/// exit in every run of the whole suite on the threaded build.)
const sameMzIons = () => {
  const theoretical = OpenMS.annotatePeaks(new Float32Array(0), new Float32Array(0), "PEPTIDEK", 3,
                                           { ions: ["a", "b", "c", "x", "y", "z"], losses: true }).theoretical;
  const [z1, y1] = ["z1++", "y1-NH3++"].map((ion) => theoretical.find((t) => t.ion === ion));
  return { theoretical, mz: z1.mz, equal: z1.mz === y1.mz };
};
const annotateAll = (mz, intensity, tolerance) => OpenMS.annotatePeaks(Float32Array.from(mz), Float32Array.from(intensity), "PEPTIDEK", 3,
  { ions: ["a", "b", "c", "x", "y", "z"], losses: true, ...tolerance });
const peakIons = (r) => r.matches.map((m) => [m.peak, m.ion]).sort((a, b) => a[0] - b[0] || (a[1] < b[1] ? -1 : 1));

test("annotatePeaks in ppm: each ion takes the closest peak, so ions of the same m/z share it", () => {
  const sameMz = sameMzIons();
  assert.ok(sameMz.equal, "z1++ and y1-NH3++ have the same m/z");
  const ppm = { tolerance: 10, unit: "ppm" };
  // one peak at their m/z: both ions have it
  assert.deepEqual(peakIons(annotateAll([sameMz.mz], [1], ppm)), [[0, "y1-NH3++"], [0, "z1++"]]);
  // two peaks of that m/z, or one more 1 ppm above: still both ions on the first
  assert.deepEqual(peakIons(annotateAll([sameMz.mz, sameMz.mz], [1, 2], ppm)), [[0, "y1-NH3++"], [0, "z1++"]]);
  assert.deepEqual(peakIons(annotateAll([sameMz.mz, sameMz.mz * (1 + 1e-6)], [1, 2], ppm)), [[0, "y1-NH3++"], [0, "z1++"]]);
  // a peak is in several matches of the theoretical spectrum as peaks: all its ions are matched, on fewer peaks
  const peaks = sameMz.theoretical.map((t) => t.mz);
  const r = annotateAll(peaks, peaks.map(() => 1), ppm);
  assert.deepEqual(r.matches.map((m) => m.ion).sort(), sameMz.theoretical.map((t) => t.ion).sort());
  assert.ok(new Set(r.matches.map((m) => m.peak)).size < r.matches.length, "several ions share a peak");
  for (const m of r.matches) assert.ok(Math.abs(peaks[m.peak] - m.theoreticalMz) / m.theoreticalMz * 1e6 <= 10, m.ion);
});

test("annotatePeaks in Da: matching is one to one, in m/z order", () => {
  const sameMz = sameMzIons();
  const da = { tolerance: 0.001 };
  // one peak at the m/z of two ions: only one of them has it
  const one = annotateAll([sameMz.mz], [1], da);
  assert.equal(one.matches.length, 1);
  assert.equal(one.matches[0].peak, 0);
  assert.ok(["z1++", "y1-NH3++"].includes(one.matches[0].ion), one.matches[0].ion);
  // two peaks: one ion each
  const two = annotateAll([sameMz.mz, sameMz.mz], [1, 2], da);
  assert.deepEqual(two.matches.map((m) => m.peak).sort(), [0, 1]);
  assert.deepEqual(two.matches.map((m) => m.ion).sort(), ["y1-NH3++", "z1++"]);
  // the theoretical spectrum as peaks: a peak has at most one ion, in order of m/z
  const peaks = sameMz.theoretical.map((t) => t.mz);
  const r = annotateAll(peaks, peaks.map(() => 1), da);
  assert.equal(new Set(r.matches.map((m) => m.peak)).size, r.matches.length, "a peak has one ion");
  for (let i = 1; i < r.matches.length; ++i) {
    assert.ok(r.matches[i - 1].theoreticalMz <= r.matches[i].theoreticalMz && peaks[r.matches[i - 1].peak] <= peaks[r.matches[i].peak]);
  }
});

test("annotatePeaks: peaks of the same m/z", () => {
  const y4 = Math.fround(OpenMS.annotatePeaks(new Float32Array(0), new Float32Array(0), "PEPTIDEK", 2, {})
                              .theoretical.find((t) => t.ion === "y4").mz);
  const annotate = (mz, intensity, options) => OpenMS.annotatePeaks(Float32Array.from(mz), Float32Array.from(intensity), "PEPTIDEK", 2, options);
  const ppm = { tolerance: 10, unit: "ppm" };
  const da = { tolerance: 0.01 };
  // in Da the alignment takes one of the peaks of equal m/z, which one is not specified (it depends on floating-point
  // details of the alignment and the neighbouring peaks): exactly one match, on one of them, and the same whatever the
  // intensities. @p equal: where the equal peaks are in @p peaks
  const takesOneOf = (peaks, equal, ion, options, what) => {
    const taken = [(i) => i + 1, (i) => peaks.length - i].map((intensity) => {
      const matches = peakIons(annotate(peaks, peaks.map((_, i) => intensity(i)), options));
      assert.equal(matches.length, 1, `${what}: ${JSON.stringify(matches)}`);
      assert.ok(equal.includes(matches[0][0]) && matches[0][1] === ion, `${what}: ${JSON.stringify(matches)}`);
      return matches[0][0];
    });
    assert.equal(taken[0], taken[1], `${what}: the peak does not depend on the intensities`);
  };

  // ppm: the first of the equal peaks, whatever their intensities (the sort of the peaks is stable)
  for (const n of [2, 3, 4]) {
    const peaks = new Array(n).fill(y4);
    for (const intensity of [peaks.map((_, i) => i + 1), peaks.map((_, i) => n - i)]) {
      assert.deepEqual(peakIons(annotate(peaks, intensity, ppm)), [[0, "y4"]], `${n} peaks`);
    }
    takesOneOf(peaks, peaks.map((_, i) => i), "y4", da, `Da, ${n} peaks`);
  }
  // among other peaks, and not next to each other in the arrays
  assert.deepEqual(peakIons(annotate([900, y4, 100, y4], [1, 1, 1, 1], ppm)), [[1, "y4"]]);
  assert.deepEqual(peakIons(annotate([y4, 900, y4, 100, y4], [1, 1, 1, 1, 1], ppm)), [[0, "y4"]]);
  takesOneOf([900, y4, 100, y4], [1, 3], "y4", da, "Da, around 100 and 900");
  // a smaller peak before two equal ones 4 mDa from y3 (where the alignment does not take the second), on either side
  const y3 = OpenMS.annotatePeaks(new Float32Array(0), new Float32Array(0), "PEPTIDEK", 2, {}).theoretical.find((t) => t.ion === "y3").mz;
  for (const offset of [-0.004, 0.004]) {
    const p = Math.fround(y3 + offset);
    takesOneOf([10, p, p], [1, 2], "y3", da, `Da, 10 and two peaks at y3 ${offset}`);
    assert.deepEqual(peakIons(annotate([10, p, p], [1, 1, 1], { tolerance: 20, unit: "ppm" })), [[1, "y3"]], `ppm, 10 and two peaks at y3 ${offset}`);
  }
});

test("annotatePeaks: the names of ions with neutral losses, other ion types and precursor peaks", () => {
  const options = { ions: ["a", "b", "c", "x", "y", "z"], losses: true, precursor: true, tolerance: 0.001 };
  // [peptide, charge, ions that have to be there]
  const peptides = [
    ["PEPTIDEK", 3, ["b3", "b3++", "y7-H2O", "y7-H2O++", "y1-NH3", "a2", "c3++", "x4", "z4++", "b4-H2O",
                     "[M+2H]++", "[M+2H-H2O]++", "[M+2H-NH3]++"]],
    ["M(Oxidation)RESK", 4, ["b1", "b1++", "y4", "y4+++", "y4-H2O+++", "y3-NH3++", "[M+3H]+++", "[M+3H-H2O]+++"]],
    [".(Acetyl)DAHKSEVAHR", 3, ["b2-H2O", "b2-H2O++", "y3-NH3++", "y9", "[M+2H]++"]],
    ["PEPS(Phospho)TIDEK", 3, ["b4", "b4-H3PO4", "b4-H3PO4++", "y6-H3PO4", "y6-H3PO4++", "y5", "[M+2H-H2O]++"]],
    ["PEPTIDER.(Amidated)", 2, ["b7", "y1", "y7-H2O", "[M+H]", "[M+H-NH3]"]],
  ];
  for (const [peptide, charge, expected] of peptides) {
    // every ion as a peak of its own: it has to be matched, with its own name
    const theoretical = OpenMS.annotatePeaks(new Float32Array(0), new Float32Array(0), peptide, charge, options).theoretical;
    const ions = new Set(theoretical.map((t) => t.ion));
    for (const ion of expected) assert.ok(ions.has(ion), `${peptide}: ${ion} in ${[...ions]}`);
    assert.equal(ions.size, theoretical.length, `${peptide}: unique ions`);
    const mz = Float32Array.from(theoretical, (t) => t.mz);
    const { matches } = OpenMS.annotatePeaks(mz, new Float32Array(mz.length).fill(1), peptide, charge, options);
    // (in Da, ions of the same m/z, e.g. z1++ and y1-NH3++, are two peaks of that m/z: each takes one)
    assert.deepEqual(matches.map((m) => m.ion).sort(), [...ions].sort(), `${peptide}: every ion is matched`);

    const byIon = new Map(theoretical.map((t) => [t.ion, t]));
    const kinds = new Set();
    for (const m of matches) {
      const what = `${peptide}: ${m.ion}`;
      const t = byIon.get(m.ion);
      assert.ok(t, what);
      assert.deepEqual([m.type, m.ordinal, m.charge, m.theoreticalMz], [t.type, t.ordinal, t.charge, t.mz], what);
      assert.ok(Math.abs(m.errorPpm) < 2, `${what}: ${m.errorPpm} ppm (float32 m/z)`);
      assert.ok(Number.isInteger(m.charge) && m.charge >= 1 && Number.isInteger(m.ordinal), what);
      if (m.type === "precursor") {
        assert.equal(m.ordinal, 0, what);
        const [, adduct, plus] = m.ion.match(/^\[M\+(\d*)H(?:-[A-Za-z0-9]+)?\](\+*)$/) ?? assert.fail(what);
        assert.equal(adduct === "" ? 1 : Number(adduct), m.charge, what);
        assert.equal(plus.length, m.charge > 1 ? m.charge : 0, what); // one '+' per charge, none for charge 1
        assert.equal(m.ion.includes("-"), m.loss !== "", what);
        kinds.add(m.loss ? "precursor loss" : "precursor");
      } else {
        assert.match(m.type, /^[a-z]$/, what);
        assert.ok(m.ordinal >= 1, what);
        assert.equal(ionLabel(m), m.ion, what);
        kinds.add(m.loss ? `${m.type} loss` : m.type);
      }
      if (m.loss) {
        // the loss is a formula, and the ion has the m/z of the ion without the loss less the formula's mass
        const parent = m.ion.replace(`-${m.loss}`, "");
        const p = byIon.get(parent);
        assert.ok(p, `${what}: the ion without the loss, ${parent}`);
        const mass = using([OpenMS.EmpiricalFormula.fromString(m.loss)], (f) => f.getMonoWeight());
        assertClose(assert, p.mz - m.theoreticalMz, mass / m.charge, 1e-6, `${what}: loses ${m.loss}`);
      }
    }
    for (const kind of ["a", "b", "c", "x", "y", "z", "b loss", "y loss", "precursor", "precursor loss"]) {
      assert.ok(kinds.has(kind), `${peptide}: matches of ${kind}`);
    }
    // the loss of a modified residue is named by its formula, conventionally spelled for the common ones
    if (peptide === "M(Oxidation)RESK") {
      assert.ok(matches.some((m) => m.loss === "CH4OS"), "the loss of oxidized M");
    }
    if (peptide === "PEPS(Phospho)TIDEK") {
      const phospho = matches.filter((m) => m.loss === "H3PO4");
      assert.ok(phospho.length > 10, `${phospho.length} ions that lost H3PO4`);
      assert.ok(!matches.some((m) => m.loss === "H3O4P" || m.ion.includes("H3O4P")), "not in the formula's alphabetical order");
      assert.ok(phospho.some((m) => m.ion === "y6-H3PO4" && m.type === "y" && m.ordinal === 6 && m.charge === 1));
    }
  }
});

test("annotatePeaks: options, empty and invalid input", () => {
  const peptide = "PEPTIDEK";
  const none = new Float32Array(0);
  const empty = OpenMS.annotatePeaks(none, none, peptide, 2, { tolerance: 0.01 });
  assert.deepEqual(empty.matches, []);
  assert.deepEqual(empty.theoretical.map((t) => t.ion).sort(), ["b1", "b2", "b3", "b4", "b5", "b6", "b7",
                                                              "y1", "y2", "y3", "y4", "y5", "y6", "y7"].sort());

  // the defaults of annotateSpectrum: b and y ions of the charges 1 to charge - 1, without losses and precursor
  for (const options of [undefined, null, {}]) {
    const ions = OpenMS.annotatePeaks(none, none, peptide, 3, options).theoretical.map((t) => t.ion);
    assert.equal(ions.length, 28, JSON.stringify(options));
    assert.ok(ions.includes("y5++") && ions.includes("b2++") && !ions.some((i) => i.includes("+++")));
  }
  assert.equal(OpenMS.annotatePeaks(none, none, peptide, 1, null).theoretical.length, 14, "charge 1: at least charge 1");
  assert.ok(OpenMS.annotatePeaks(none, none, peptide, 2, { maxCharge: 3 }).theoretical.some((t) => t.ion === "y5+++"));
  assert.deepEqual([...new Set(OpenMS.annotatePeaks(none, none, peptide, 2, { ions: ["y", "a"] }).theoretical.map((t) => t.type))].sort(),
                   ["a", "y"]);

  // a y4 ion 5 mDa off, and a peak 1 Da away from the b2 ion
  const y4 = empty.theoretical.find((t) => t.ion === "y4").mz;
  const b2 = empty.theoretical.find((t) => t.ion === "b2").mz;
  const mz = Float32Array.of(y4 + 0.005, b2 + 1);
  const intensity = Float32Array.of(10, 20);
  const da = OpenMS.annotatePeaks(mz, intensity, peptide, 2, { tolerance: 0.01 });
  assert.deepEqual(da.matches.map((m) => [m.peak, m.ion]), [[0, "y4"]]);
  // ppm whatever the unit of the tolerance
  assertClose(assert, da.matches[0].errorPpm, 0.005 / y4 * 1e6, 0.5, "5 mDa as ppm");
  assert.ok(da.matches[0].errorPpm > 0);
  assert.deepEqual(OpenMS.annotatePeaks(mz, intensity, peptide, 2, { tolerance: 20, unit: "ppm" }).matches.map((m) => m.ion), ["y4"]);
  assertClose(assert, OpenMS.annotatePeaks(mz, intensity, peptide, 2, { tolerance: 20, unit: "ppm" }).matches[0].errorPpm,
              da.matches[0].errorPpm, 1e-9, "the same error in either unit");
  assert.deepEqual(OpenMS.annotatePeaks(mz, intensity, peptide, 2, { tolerance: 5, unit: "ppm" }).matches, [], "5 mDa is more than 5 ppm");
  assert.deepEqual(OpenMS.annotatePeaks(mz, intensity, peptide, 2, { tolerance: 1.5 }).matches.map((m) => m.ion).sort(), ["b2", "y4"]);

  const rejects = (what, f, pattern) => assert.throws(f, pattern, what);
  rejects("charge 0", () => OpenMS.annotatePeaks(mz, intensity, peptide, 0, null), /InvalidValue/);
  rejects("charge -2", () => OpenMS.annotatePeaks(mz, intensity, peptide, -2, null), /InvalidValue/);
  rejects("empty peptide", () => OpenMS.annotatePeaks(mz, intensity, "", 2, null), /InvalidValue/);
  // (AASequence.fromString's own errors: a syntax error is a ParseError, an unknown residue an InvalidValue)
  rejects("a peptide with a digit", () => OpenMS.annotatePeaks(mz, intensity, "PEPT1DE", 2, null), /Parse Error.*PEPT1DE/);
  rejects("a modification that is not closed", () => OpenMS.annotatePeaks(mz, intensity, "PEPTIDE(", 2, null), /Parse Error/);
  rejects("an unknown residue", () => OpenMS.annotatePeaks(mz, intensity, "PEPTXDE", 2, null), /InvalidValue/);
  rejects("unit", () => OpenMS.annotatePeaks(mz, intensity, peptide, 2, { unit: "mmu" }), /unit/);
  rejects("ion type", () => OpenMS.annotatePeaks(mz, intensity, peptide, 2, { ions: ["q"] }), /ion type/);
  rejects("lengths", () => OpenMS.annotatePeaks(mz, Float32Array.of(1), peptide, 2, null), /IllegalArgument/);
  rejects("NaN", () => OpenMS.annotatePeaks(Float32Array.of(100, NaN), Float32Array.of(1, 1), peptide, 2, null), /InvalidValue/);
});

/// Runs @p fn with the PeptDeep models of the package (models/ next to openms.mjs) mounted at /models; the default
/// model of predictMS2() is /models/peptdeep_ms2_dynamic.onnx.
function withModels(fn, mountPoint = "/models") {
  const models = path.join(path.dirname(modulePath), "models");
  assert.ok(fs.existsSync(path.join(models, "peptdeep_ms2_dynamic.onnx")), `the models are in ${models}`);
  OpenMS.mountHostDirectory(mountPoint, models);
  try {
    return fn();
  } finally {
    OpenMS.unmount(mountPoint);
  }
}

const avmPeptide = "AVM(Oxidation)DDFAAFVEK";

/// The `ion` label of a predicted ion: as annotatePeaks names it
const predictedLabel = (i) => `${i.type}${i.ordinal}${i.charge > 1 ? "+".repeat(i.charge) : ""}`;

// A module with ONNX Runtime (WITH_ONNX): the b and y intensities that PeptDeep predicts, for the mirror plot of a PSM
test("predictMS2 predicts b and y intensities", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  withModels(() => {
    const ions = OpenMS.predictMS2(avmPeptide, 2, 30, "QE");
    assert.ok(Array.isArray(ions) && ions.length > 0, "ions are predicted");
    const length = 12; // residues of AVMDDFAAFVEK: ordinals 1 to 11
    for (const ion of ions) {
      assert.deepEqual(Object.keys(ion).sort(), ["charge", "intensity", "ion", "mz", "ordinal", "type"]);
      assert.ok(ion.type === "b" || ion.type === "y", ion.type);
      assert.ok(ion.charge === 1 || ion.charge === 2, `charge ${ion.charge}`);
      assert.ok(Number.isInteger(ion.ordinal) && ion.ordinal >= 1 && ion.ordinal <= length - 1, `ordinal ${ion.ordinal}`);
      assert.ok(ion.intensity > 0.01 && ion.intensity <= 1, `${ion.ion}: intensity ${ion.intensity}`);
      assert.equal(ion.ion, predictedLabel(ion), "the ion is named as annotatePeaks names it");
    }
    // relative to the most intense returned b/y ion, which is exactly 1
    assert.equal(Math.max(...ions.map((i) => i.intensity)), 1);
    assert.equal(new Set(ions.map((i) => i.ion)).size, ions.length, "every ion once");
    // by ascending m/z
    for (let k = 1; k < ions.length; ++k) assert.ok(ions[k - 1].mz <= ions[k].mz, "sorted by m/z");
    // y1 of the C-terminal K
    const y1 = ions.find((i) => i.ion === "y1");
    assert.ok(y1, "y1 is predicted");
    assertClose(assert, y1.mz, 147.1128, 1e-3, "m/z of y1");
    // both b and y ions, of charge 1 and 2
    assert.deepEqual([...new Set(ions.map((i) => i.type))].sort(), ["b", "y"]);
    assert.deepEqual([...new Set(ions.map((i) => i.charge))].sort(), [1, 2]);
    const doubly = ions.filter((i) => i.charge === 2);
    assert.ok(doubly.length > 0 && doubly.every((i) => i.ion.endsWith("++") && !i.ion.endsWith("+++")), "doubly charged ions are named with '++'");
    assert.ok(ions.filter((i) => i.charge === 1).every((i) => !i.ion.includes("+")), "singly charged ions without '+'");

    // m/z from the theoretical spectrum of annotatePeaks (the modification of the peptide string is in it)
    const none = new Float32Array(0);
    const theoretical = new Map(OpenMS.annotatePeaks(none, none, avmPeptide, 3, { ions: ["b", "y"] }).theoretical.map((i) => [i.ion, i.mz]));
    for (const ion of ions) assertClose(assert, ion.mz, theoretical.get(ion.ion), 1e-9, `m/z of ${ion.ion}`);
    // and from the residues, as PeptDeepLibraryPredictor computes them: b(n) is the first n residues, y(n) the last n
    using([OpenMS.AASequence.fromString(avmPeptide)], (sequence) => {
      assert.equal(sequence.size(), length);
      for (const ion of ions) {
        const part = ion.type === "b" ? sequence.getPrefix(ion.ordinal) : sequence.getSuffix(ion.ordinal);
        const mz = using([part], (p) => p.getMZ(ion.charge, ion.type === "b" ? OpenMS.ResidueType.BIon : OpenMS.ResidueType.YIon));
        assertClose(assert, ion.mz, mz, 1e-6, `m/z of ${ion.ion} from its residues`);
      }
    });

    // what the model of the package (peptdeep_ms2_dynamic.onnx) predicts for this peptide: y9 is the base peak, then b2
    // (a regression pin of the model's output, to be updated with the model; the orientation itself is pinned by the
    // proline effect in the next test)
    const strongest = [...ions].sort((a, b) => b.intensity - a.intensity);
    assert.equal(strongest[0].ion, "y9");
    assert.equal(strongest[1].ion, "b2");
    assertClose(assert, strongest[1].intensity, 0.89, 0.05, "intensity of b2");

    // the model is the default one, and the answer is the same for an explicit path
    assert.deepEqual(OpenMS.predictMS2(avmPeptide, 2, 30, "QE", "/models/peptdeep_ms2_dynamic.onnx"), ions);
  });
});

// The channels of the model are read with the y ions numbered from the C-terminus. That follows from chemistry, not
// from the code: the proline effect, cleavage N-terminal to a proline dominates the spectrum of a peptide, which makes
// the y ion that begins with the proline the strongest y ion.
test("predictMS2: the y ions are numbered from the C-terminus (proline effect)", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  withModels(() => {
    // proline effect: cleavage N-terminal to P dominates; a reversed y numbering would name these y5, y5, y2
    for (const [pep, y] of [["AGLQFPVGR", 4], ["VEADIPGHGQEVLIR", 10], ["ELPEHTVK", 6]]) {
      const ys = OpenMS.predictMS2(pep, 2, 30, "QE").filter((i) => i.type === "y" && i.charge === 1);
      assert.equal(ys.reduce((a, b) => (b.intensity > a.intensity ? b : a)).ordinal, y, pep);
      assert.equal(pep[pep.length - y], "P", `${pep}: y${y} begins with the proline`);
    }
  });
});

test("predictMS2: charge, instrument and collision energy are inputs of the model", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  withModels(() => {
    const peptide = "LHVNSAR";
    const at = (charge, nce, instrument, sequence = peptide) => OpenMS.predictMS2(sequence, charge, nce, instrument);
    // a precursor of charge 1 has no doubly charged fragments; others can have them (this peptide's, for the charges 2 to 4)
    for (const charge of [2, 3, 4]) {
      assert.ok(at(charge, 30, "QE", avmPeptide).some((i) => i.charge === 2), `charge ${charge}`);
    }
    // (PEPTIDEK and EDPQTFYYAVAVVK: the model gives doubly charged ions of 7 % and 81 % to a precursor of charge 1)
    for (const sequence of [peptide, avmPeptide, "PEPTIDEK", "EDPQTFYYAVAVVK", "SLYVDEGGHNLTQVHKGEE"]) {
      assert.ok(at(1, 30, "QE", sequence).length > 0);
      assert.ok(at(1, 30, "QE", sequence).every((i) => i.charge === 1), `${sequence}, charge 1: no doubly charged fragments`);
    }
    // the intensities are relative to the strongest returned ion: in these (to a precursor of charge 1) the model's
    // own maximum is a doubly charged or a modloss ion, 18 % and 20 % of which is the strongest b or y ion
    for (const sequence of ["SDNKVPFVGWQVRAVYSGK", "NILMDAFYEPMQPLRQCR"]) {
      const ions = at(1, 30, "QE", sequence);
      assert.ok(ions.length > 1);
      assert.equal(Math.max(...ions.map((i) => i.intensity)), 1, sequence);
      assert.ok(ions.every((i) => i.intensity > 0.01 && i.intensity <= 1 && i.charge === 1));
    }
    for (const charge of [1, 2, 3]) {
      const ions = at(charge, 30, "QE");
      assert.equal(Math.max(...ions.map((i) => i.intensity)), 1, `charge ${charge}`);
      assert.ok(ions.every((i) => i.ordinal >= 1 && i.ordinal <= 6));
    }
    // the instrument is peptdeep's, case-insensitive; a name it does not know is the catch-all, as in ML::instrumentIndex
    assert.deepEqual(at(2, 30, "qe"), at(2, 30, "QE"));
    assert.deepEqual(at(2, 30, "TIMSTOF"), at(2, 30, "timsTOF"));
    assert.deepEqual(at(2, 30, "no such instrument"), at(2, 30, ""));
    assert.deepEqual(at(2, 30, "Orbitrap Exploris"), at(2, 30, "no such instrument"));
    const differs = (a, b) => JSON.stringify(a.map((i) => i.intensity)) !== JSON.stringify(b.map((i) => i.intensity));
    assert.ok(differs(at(2, 30, "QE"), at(2, 30, "timsTOF")), "the instrument changes the prediction");
    assert.ok(differs(at(2, 30, "QE"), at(2, 30, "no such instrument")), "so does the catch-all");
    assert.ok(differs(at(2, 20, "QE"), at(2, 40, "QE")), "the collision energy changes the prediction");
    assert.ok(differs(at(2, 30, "QE"), at(3, 30, "QE")), "the charge too");
  });
});

test("predictMS2 loads a model once per path, and a failed load is not kept", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  const model = "/models-late/peptdeep_ms2_dynamic.onnx";
  // not mounted yet: the library's exception for a missing file
  assert.throws(() => OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", model), /FileNotFound.*models-late/);
  assert.throws(() => OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", model), /FileNotFound.*models-late/, "the failure is not kept");
  withModels(() => {
    const expected = OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE");
    withModels(() => {
      const time = () => {
        const start = performance.now();
        const ions = OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", model);
        return { ions, ms: performance.now() - start };
      };
      // now it is there: loaded by the first call, a session of its own (not that of the default path)
      const first = time();
      assert.deepEqual(first.ions, expected, "the same model gives the same prediction");
      const later = Array.from({ length: 7 }, time);
      for (const { ions } of later) assert.deepEqual(ions, expected);
      const median = later.map((r) => r.ms).sort((a, b) => a - b)[3];
      t.diagnostic(`first call ${first.ms.toFixed(0)} ms, median of 7 later calls ${median.toFixed(1)} ms`);
      assert.ok(first.ms > 2 * median, `the first call loads the model (${first.ms} ms), later ones do not (median ${median} ms)`);
    }, "/models-late");
  });
});

test("predictMS2: errors", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  const rejects = (what, f, pattern) => assert.throws(f, pattern, what);
  withModels(() => {
    rejects("a model that is not there", () => OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", "/models/no-such-model.onnx"), /FileNotFound.*no-such-model/);
    rejects("charge 0", () => OpenMS.predictMS2("PEPTIDEK", 0, 30, "QE"), /InvalidValue/);
    rejects("charge -1", () => OpenMS.predictMS2("PEPTIDEK", -1, 30, "QE"), /InvalidValue/);
    rejects("NaN as energy", () => OpenMS.predictMS2("PEPTIDEK", 2, NaN, "QE"), /InvalidValue/);
    rejects("Infinity as energy", () => OpenMS.predictMS2("PEPTIDEK", 2, Infinity, "QE"), /InvalidValue/);
    rejects("empty peptide", () => OpenMS.predictMS2("", 2, 30, "QE"), /InvalidValue/);
    // (AASequence.fromString's own errors: a syntax error is a ParseError, an unknown residue an InvalidValue)
    rejects("a peptide with a digit", () => OpenMS.predictMS2("PEPT1DE", 2, 30, "QE"), /Parse Error.*PEPT1DE/);
    rejects("an unknown residue", () => OpenMS.predictMS2("PEPTXDE", 2, 30, "QE"), /InvalidValue/);
    rejects("an unknown modification", () => OpenMS.predictMS2("PEPT(NoSuchMod)IDE", 2, 30, "QE"), /NoSuchMod/);
    // an omitted model path may be passed as undefined or null: the default model; other types are rejected
    const expected = OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE");
    assert.deepEqual(OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", undefined), expected);
    assert.deepEqual(OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", null), expected);
    rejects("a number as model path", () => OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", 5), /InvalidValue.*number/);
    rejects("an empty model path", () => OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE", ""), /FileNotFound/);
    // a peptide of one residue has no cleavage
    assert.deepEqual(OpenMS.predictMS2("K", 2, 30, "QE"), []);
    // and the errors leave the module as it was
    assert.ok(OpenMS.predictMS2("PEPTIDEK", 2, 30, "QE").length > 0);
  });
});

test("predictMS2 throws without ONNX Runtime", (t) => {
  if (OpenMS.hasONNXRuntime()) return t.skip("the module was built with ONNX Runtime");
  assert.throws(() => OpenMS.predictMS2(avmPeptide, 2, 30, "QE"), /without ONNX Runtime/);
  assert.throws(() => OpenMS.predictMS2(avmPeptide, 2, 30, "QE", "/models/peptdeep_ms2_dynamic.onnx"), /without ONNX Runtime/);
  assert.throws(() => OpenMS.predictMS2(avmPeptide, 2, 30, "QE", undefined), /without ONNX Runtime/);
});

// ProSE's PeptDeep rescoring records the collision energy its predictions used for each file
test("searchFiles reports the PeptDeep NCE", (t) => {
  if (!OpenMS.hasONNXRuntime()) return t.skip("the module was built without ONNX Runtime");
  withModels(() => {
    const parameters = {
      ...prose6,
      "peptdeep:enable": true,
      "peptdeep:ms2_model": "/models/peptdeep_ms2_dynamic.onnx", "peptdeep:rt_model": "/models/peptdeep_rt_dynamic.onnx",
    };
    withSpectraCopies(OpenMS, ["a.mzML", "b.mzML"], (files) => {
      const set = searchFiles(OpenMS, files, undefined, { parameters });
      try {
        const { files: summaries } = set.summary();
        assert.equal(summaries.length, 2);
        for (const f of summaries) {
          assert.equal(f.status, "EXECUTION_OK");
          assert.equal(typeof f.peptdeepNce, "number", `NCE of ${f.path}`);
          assert.ok(f.peptdeepNce >= 20 && f.peptdeepNce <= 40, `NCE ${f.peptdeepNce} of ${f.path}`);
        }
      } finally {
        set.delete();
      }
    });
  });
});

test("getMemoryInfo", () => {
  const { heapSize, heapMax } = OpenMS.getMemoryInfo();
  assert.ok(heapSize > 0 && heapSize <= heapMax);
});
