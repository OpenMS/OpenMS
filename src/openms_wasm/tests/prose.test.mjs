// The search engine ProSE in WebAssembly: ProSEAlgorithm with the settings and data of the TOPP
// test TOPP_ProSE_1 (src/tests/topp), which has to give the same PSMs; its progress reports,
// statistics, peak annotations and open search; and annotateSpectrum() for spectrum viewers.

import { test } from "node:test";
import assert from "node:assert/strict";
import { loadOpenMS, using, assertClose } from "./helpers.mjs";

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

  // the data and settings of TOPP_ProSE_11_rescore (ProSE_6_fdr.ini): rescored, q-values as scores, filtered at FDR:PSM
  using([new OpenMS.ProSEAlgorithm(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (prose, proteins, peptides) => {
      prose.setParameterValues({
        "enzyme": "Trypsin", "decoys": "generate",
        "precursor:mass_tolerance_lower": 20, "precursor:mass_tolerance_upper": 20, "precursor:mass_tolerance_unit": "ppm",
        "fragment:mass_tolerance": 10, "fragment:mass_tolerance_unit": "ppm",
        "modifications:fixed": ["Carbamidomethyl (C)"], "modifications:variable": ["Oxidation (M)"],
        "FDR:PSM": 0.05, "FDR:protein": 0.05,
      });
      const summary = prose.search("/topp/THIRDPARTY/DatabaseSuitability_in_spec.mzML", "/topp/THIRDPARTY/DatabaseSuitability_database.fasta",
        proteins, peptides, { rescore: true });
      assert.equal(summary.status, "EXECUTION_OK");
      assert.equal(summary.percolator.applied, true, summary.percolator.skipped ?? "");
      const ids = peptides.toArray().filter((id) => id.hits.length > 0);
      assert.ok(ids.length > 0);
      for (const id of ids) {
        assert.equal(id.scoreType, "q-value");
        assert.ok(id.hits[0].score <= 0.05, `q-value ${id.hits[0].score} passes FDR:PSM`);
      }
      assert.equal(proteins.toArray()[0].searchEngine, "Percolator");
      // the parameters are as set (the search ran with the FDR filters deferred)
      assert.equal(prose.getParameterValues()["FDR:PSM"], 0.05);
    });
});

test("ProSEAlgorithm reports its progress", () => {
  const reports = [];
  const { summary } = search(prose1, (p) => reports.push(p));
  assert.equal(summary.status, "EXECUTION_OK");
  const started = reports.filter((r) => r.event === "start").map((r) => r.label);
  assert.ok(started.includes("Building fragment index..."), JSON.stringify(started));
  assert.ok(reports.some((r) => r.event === "end"));
  for (const r of reports) {
    assert.ok(["start", "progress", "end"].includes(r.event));
    assert.ok(r.depth >= 0 && r.begin <= r.end);
  }
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
