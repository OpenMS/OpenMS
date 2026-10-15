// Reading and writing files: through the virtual file system and from strings.

import { test } from "node:test";
import assert from "node:assert/strict";
import fs from "node:fs";
import path from "node:path";
import { loadOpenMS, using, assertClose, testDataDir } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

/// MzMLFile_test: MzMLFile_1.mzML has 4 spectra (15, 10, 15 and 0 peaks) and 2 chromatograms
function checkMzMLFile1(exp) {
  assert.equal(exp.getNrSpectra(), 4);
  assert.equal(exp.getNrChromatograms(), 2);
  using([exp.getSpectrum(0), exp.getSpectrum(1)], (s0, s1) => {
    assert.equal(s0.size(), 15);
    assert.equal(s1.size(), 10);
    assert.equal(s0.getMSLevel(), 1);
    assert.equal(s1.getMSLevel(), 2);
    assertClose(assert, s0.getRT(), 5.1, 1e-6, "RT");
    const { mz, intensity } = s0.getPeaks();
    assert.ok(mz instanceof Float64Array && intensity instanceof Float32Array);
    assertClose(assert, mz[0], 0, 1e-6, "first m/z");
    assertClose(assert, intensity[0], 15, 1e-6, "first intensity");
    const precursors = s1.getPrecursors();
    assert.equal(precursors.length, 2);
    // the selected ion m/z (5.5555); MzMLFile_test compares it to the isolation window target 5.55 with tolerance 0.01
    assertClose(assert, precursors[0].mz, 5.5555, 1e-6, "precursor m/z");
    assert.equal(precursors[0].charge, 2);
    assertClose(assert, precursors[0].intensity, 120053, 1e-3, "precursor intensity");
    assertClose(assert, precursors[0].driftTime, 8.1, 1e-6, "precursor drift time");
    assert.ok(precursors[0].activationMethods.includes("CID"));
    assert.equal(s0.getNativeID(), "index=0");
    assert.equal(s0.getMetaValue("maldi_spot_id"), "M0");
    assert.equal(s0.getType(), "Centroid");
  });
}

test("mountHostDirectory and unmount", () => {
  OpenMS.mountHostDirectory("/mounted", testDataDir);
  assert.ok(OpenMS.FS.stat("/mounted/MzMLFile_1.mzML").size > 1000);
  assert.ok(OpenMS.FS.readdir("/mounted").includes("MzMLFile_1.mzML"));
  OpenMS.unmount("/mounted");
  assert.throws(() => OpenMS.FS.stat("/mounted/MzMLFile_1.mzML"));
  assert.throws(() => OpenMS.unmount("/mounted"), /mounted/);
  assert.throws(() => OpenMS.mountHostDirectory("/no/such/parent/x", testDataDir), /no\/such\/parent/);
});

test("MzMLFile from the (mounted) file system", () => {
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment()], (file, exp) => {
    file.load("/data/MzMLFile_1.mzML", exp);
    checkMzMLFile1(exp);
    const summary = exp.getSpectrumSummary();
    assert.equal(summary.length, 4);
    assert.equal(summary[1].msLevel, 2);
    assert.deepEqual(exp.getMSLevels(), [1, 2]);
  });
});

test("MzMLFile from a string, and round trip through MEMFS", () => {
  const text = fs.readFileSync(path.join(testDataDir, "MzMLFile_1.mzML"), "utf8");
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment()],
    (file, exp, fromFile, fromString) => {
      file.loadBuffer(text, exp);
      checkMzMLFile1(exp);

      file.store("/tmp/roundtrip.mzML", exp);
      assert.ok(OpenMS.FS.stat("/tmp/roundtrip.mzML").size > 1000);
      file.load("/tmp/roundtrip.mzML", fromFile);
      checkMzMLFile1(fromFile);

      file.loadBuffer(file.storeBuffer(exp), fromString);
      checkMzMLFile1(fromString);
    });
});

test("MzMLFile load options", () => {
  using([new OpenMS.MzMLFile(), new OpenMS.MSExperiment()], (file, exp) => {
    file.setOptions({ msLevels: [2] });
    assert.deepEqual(file.getOptions().msLevels, [2]);
    file.load("/data/MzMLFile_1.mzML", exp);
    assert.deepEqual(exp.getMSLevels(), [2]);
  });
});

test("FileHandler", () => {
  assert.equal(OpenMS.FileHandler.getType("/data/MzMLFile_1.mzML"), "mzML");
  assert.equal(OpenMS.FileHandler.getTypeByFileName("x.featureXML"), "featureXML");
  using([new OpenMS.FileHandler(), new OpenMS.MSExperiment(), new OpenMS.FeatureMap()], (handler, exp, features) => {
    handler.loadExperiment("/data/MzMLFile_1.mzML", exp);
    checkMzMLFile1(exp);
    handler.loadFeatures("/data/FeatureXMLFile_1.featureXML", features);
    assert.equal(features.size(), 2);
  });
});

test("FASTAFile", () => {
  using([new OpenMS.FASTAFile()], (file) => {
    const entries = file.load("/data/FASTAFile_test.fasta");
    assert.ok(entries.length > 0);
    assert.ok(entries.every((e) => typeof e.identifier === "string" && e.sequence.length > 0));
    file.store("/tmp/out.fasta", entries);
    assert.deepEqual(file.load("/tmp/out.fasta"), entries);
  });
});

test("IdXMLFile", () => {
  using([new OpenMS.IdXMLFile(), new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()],
    (file, proteins, peptides) => {
      file.load("/data/IdXMLFile_whole.idXML", proteins, peptides);
      // IdXMLFile_test: 2 protein and 3 peptide identification runs
      assert.equal(proteins.size(), 2);
      assert.equal(peptides.size(), 3);
      const ids = peptides.toArray();
      assert.equal(ids[0].hits.length, 2);
      assert.equal(ids[0].hits[0].sequence, "PEPTIDER");
      assertClose(assert, ids[0].hits[0].score, 0.9, 1e-9, "score");
      assert.equal(ids[0].hits[0].charge, 1);
      const runs = proteins.toArray();
      assert.equal(runs[0].searchEngine, "Mascot");
      assert.equal(runs[0].scoreType, "MOWSE");
      assert.equal(runs[0].metaValues.name, "ProteinIdentification");
      using([peptides.get(0)], (id) => {
        assert.equal(id.getHitCount(), 2);
        const hits = id.getHits();
        try {
          assert.equal(hits[0].getSequenceString(), "PEPTIDER");
        } finally {
          hits.forEach((h) => h.delete());
        }
      });
      file.store("/tmp/out.idXML", proteins, peptides);
      using([new OpenMS.ProteinIdentificationList(), new OpenMS.PeptideIdentificationList()], (p2, pep2) => {
        file.load("/tmp/out.idXML", p2, pep2);
        assert.deepEqual(pep2.toArray().map((i) => i.hits.map((h) => h.sequence)),
                         ids.map((i) => i.hits.map((h) => h.sequence)));
      });
    });
});

test("FeatureXMLFile", () => {
  using([new OpenMS.FeatureXMLFile(), new OpenMS.FeatureMap()], (file, features) => {
    file.load("/data/FeatureXMLFile_1.featureXML", features);
    assert.equal(features.size(), 2);
    const table = features.getFeatureTable();
    assert.equal(table.rt.length, 2);
    using([features.get(0)], (f) => {
      assertClose(assert, f.getRT(), table.rt[0], 1e-9, "RT");
      assertClose(assert, f.getMZ(), table.mz[0], 1e-9, "m/z");
    });
  });
});

test("MSExperiment.getSpectrumSummary: m/z range of peaks not sorted by m/z", () => {
  using([new OpenMS.MSExperiment(), new OpenMS.MSSpectrum()], (exp, spectrum) => {
    spectrum.setPeaks([500, 100, 300], [1, 2, 3]);
    exp.addSpectrum(spectrum);
    const summary = exp.getSpectrumSummary();
    assert.equal(summary.length, 1);
    assert.equal(summary[0].minMZ, 100);
    assert.equal(summary[0].maxMZ, 500);
  });
});

test("MSExperiment.get2DPeakData: MS level 0 is all levels", () => {
  using([new OpenMS.MSExperiment(), new OpenMS.MSSpectrum(), new OpenMS.MSSpectrum()], (exp, ms1, ms2) => {
    ms1.setRT(1);
    ms1.setMSLevel(1);
    ms1.setPeaks([100], [10]);
    ms2.setRT(2);
    ms2.setMSLevel(2);
    ms2.setPeaks([200], [20]);
    exp.addSpectrum(ms1);
    exp.addSpectrum(ms2);
    // (minRT, maxRT, minMZ, maxMZ, msLevel)
    const all = exp.get2DPeakData(0, 10, 0, 1000, 0);
    assert.deepEqual(Array.from(all.rt), [1, 2]);
    assert.deepEqual(Array.from(all.mz), [100, 200]);
    const ms2Only = exp.get2DPeakData(0, 10, 0, 1000, 2);
    assert.deepEqual(Array.from(ms2Only.rt), [2]);
    assert.deepEqual(Array.from(ms2Only.mz), [200]);
  });
});
