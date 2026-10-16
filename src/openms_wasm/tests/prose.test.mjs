// The search engine ProSE in WebAssembly: ProSEAlgorithm with the settings and data of the TOPP
// test TOPP_ProSE_1 (src/tests/topp), which has to give the same PSMs; its progress reports,
// statistics, peak annotations and open search; the search of several files pooled as the ProSE
// tool's -out_merged (TOPP_ProSE_13_merged), with its results as typed columns (tables()) and the
// peaks of the identified spectra (spectra()); and annotateSpectrum() for spectrum viewers.

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

/// Checks progress @p reports of a search of the files with the steps @p fileSteps (their labels, in input order):
/// - the outermost steps (depth 0) are those of the files, each with the file's position (file);
/// - every report of a step carries its label, range and depth; its values lie in the range and never decrease;
/// - every report carries the position of the file whose step encloses it.
function assertProgressSteps(reports, fileSteps) {
  assert.deepEqual(reports.filter((r) => r.event === "start" && r.depth === 0).map((r) => [r.label, r.file]),
    fileSteps.map((label, i) => [label, i]));
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
      for (const run of [
        () => using([prose.searchFiles(["/topp/SimpleSearchEngine_1.mzML"], "/topp/SimpleSearchEngine_1.fasta", { rescore: false })], () => {}),
        () => prose.search("/topp/SimpleSearchEngine_1.mzML", "/topp/SimpleSearchEngine_1.fasta", proteins, peptides, { rescore: false }),
      ]) {
        prose.setParameterValues(failing);
        reports = [];
        assert.throws(run, /no-such-model/);
        assert.ok(reports.some((r) => r.event === "start" && r.label === "Adding PeptDeep rescoring features..."),
          "the search threw inside a step of the file's step");
        prose.setParameterValues({ "peptdeep:enable": false });
        reports = [];
        run();
        assertProgressSteps(reports, steps);
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
  const { reports } = pooled();
  assertProgressSteps(reports, ["Searching a.mzML (1 of 2)...", "Searching b.mzML (2 of 2)..."]);
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
        [["Searching run.mzML (1 of 2)...", 0], ["Searching run.mzML (2 of 2)...", 1]]);
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
  const set = searchFiles(OpenMS, ["/topp/SimpleSearchEngine_1.mzML"], undefined);
  try {
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

test("spectra() holds the peaks of the identified spectra", () => {
  // two copies of one file: PSMs of both with the same native IDs, whose spectra are those of their own file
  withSpectraCopies(OpenMS, ["a.mzML", "b.mzML"], (files) => {
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

      // the spectra of the files, as read: the peaks sorted by m/z, the precursor of the PSM, and no collision energy
      // (the file records none)
      const expected = fileSpectra(files[0]); // b.mzML is a copy of a.mzML
      const keyOf = new Map(); // spectrum -> its file and native ID
      for (let i = 0; i < psms.count; ++i) {
        const k = s.psmSpectrum[i];
        const key = `${psms.file[i]} ${psms.nativeId[i]}`;
        assert.ok(!keyOf.has(k) || keyOf.get(k) === key, `spectrum ${k}: ${keyOf.get(k)} and ${key}`);
        keyOf.set(k, key);
        const e = expected.get(psms.nativeId[i]);
        assert.ok(e && e.msLevel === 2, psms.nativeId[i]);
        assert.deepEqual(storedPeaks(s, k), e.peaks, `the peaks of ${key}`);
        assertClose(assert, s.precursorMz[k], psms.mz[i], 1e-4, `precursor m/z of ${key}`);
        assert.equal(s.precursorMz[k], e.precursors[0].mz, `precursor m/z of ${key}`);
        assert.ok(Number.isNaN(s.collisionEnergy[k]), `no collision energy for ${key}: ${s.collisionEnergy[k]}`);
      }
      assert.equal(keyOf.size, n, "every spectrum is that of a PSM");
      assert.equal(new Set(keyOf.values()).size, n, "one spectrum per file and native ID");
      // the same native ID in both files: two spectra (with the same peaks)
      const inB = new Map();
      for (let i = 0; i < psms.count; ++i) if (psms.file[i] === 1) inB.set(psms.nativeId[i], s.psmSpectrum[i]);
      let both = 0;
      for (let i = 0; i < psms.count; ++i) {
        if (psms.file[i] !== 0 || !inB.has(psms.nativeId[i])) continue;
        const [ka, kb] = [s.psmSpectrum[i], inB.get(psms.nativeId[i])];
        assert.notEqual(ka, kb, psms.nativeId[i]);
        assert.deepEqual(storedPeaks(s, ka), storedPeaks(s, kb), psms.nativeId[i]);
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
function searchSpectraToKeep(options) {
  storeSpectraToKeep("/searchfiles/keep.mzML");
  return searchFiles(OpenMS, ["/searchfiles/keep.mzML"], options, { parameters: prose1, database: "/topp/SimpleSearchEngine_1.fasta" });
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
    const set = searchSpectraToKeep({ keepSpectra: false });
    try {
      const s = set.spectra();
      assert.equal(assertSpectraShapes(s, 3), 0);
      assert.deepEqual(Array.from(s.psmSpectrum), [-1, -1, -1]);
      assert.deepEqual(Array.from(s.offsets), [0]);
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

test("getMemoryInfo", () => {
  const { heapSize, heapMax } = OpenMS.getMemoryInfo();
  assert.ok(heapSize > 0 && heapSize <= heapMax);
});
