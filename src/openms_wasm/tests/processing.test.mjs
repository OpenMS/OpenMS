// Parameters, signal processing and feature detection. The expected results are those of
// the OpenMS class tests (PeakPickerHiRes_test, FeatureFindingMetabo_test), i.e. the
// WebAssembly build has to reproduce what the native library computes.

import { test } from "node:test";
import assert from "node:assert/strict";
import { loadOpenMS, using, assertClose } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

function loadMzML(path, exp) {
  using([new OpenMS.MzMLFile()], (file) => file.load(path, exp));
}

/// TEST_REAL_SIMILAR's default: relative 1e-5 or absolute 1e-5
function assertSimilar(actual, expected, message) {
  const ok = Math.abs(actual - expected) <= 1e-5 || Math.abs(actual - expected) <= 1e-5 * Math.abs(expected);
  assert.ok(ok, `${message}: expected ${expected}, got ${actual}`);
}

test("Param", () => {
  using([new OpenMS.Param()], (p) => {
    p.setValue("int", 3);
    p.setValue("double", 1.5);
    p.setValue("string", "abc");
    p.setValue("list", [1.5, 2.5]);
    p.setValue("section:flag", true, "a boolean");
    assert.equal(p.size(), 5);
    assert.deepEqual(p.keys(), ["int", "double", "string", "list", "section:flag"]);
    assert.equal(p.getValue("int"), 3);
    assert.equal(p.getValue("double"), 1.5);
    assert.equal(p.getValue("string"), "abc");
    assert.deepEqual(p.getValue("list"), [1.5, 2.5]);
    assert.equal(p.getValue("section:flag"), "true");
    assert.equal(p.getDescription("section:flag"), "a boolean");
    // a value is converted to the type of the entry it replaces
    p.setValue("double", 2);
    assert.equal(p.getEntry("double").type, "double");
    assert.throws(() => p.setValue("int", 2.5), /Integer parameter value expected/);
    assert.throws(() => p.getValue("missing"), /ElementNotFound/);
    assert.throws(() => p.update({ missing: 1 }), /Unknown parameter/);
  });
});

test("DefaultParamHandler: parameters of an algorithm", () => {
  using([new OpenMS.PeakPickerHiRes()], (picker) => {
    assert.equal(picker.getName(), "PeakPickerHiRes");
    using([picker.getDefaults()], (defaults) => {
      const sn = defaults.getEntry("signal_to_noise");
      assert.equal(sn.type, "double");
      assert.ok(sn.description.length > 0);
      assert.ok(defaults.getEntries().length > 5);
    });
    picker.setParameterValues({ signal_to_noise: 4 });
    assert.equal(picker.getParameterValues().signal_to_noise, 4);
    // restrictions are checked when the parameters are set
    assert.throws(() => picker.setParameterValues({ signal_to_noise: -1 }), /InvalidParameter|Error/);
  });
});

test("PeakPickerHiRes on a profile spectrum", () => {
  using([new OpenMS.PeakPickerHiRes(), new OpenMS.MSSpectrum(), new OpenMS.MSSpectrum()], (picker, input, output) => {
    picker.setParameterValues({ signal_to_noise: 0 });
    input.setPeaks([100.0, 100.01, 100.02, 100.03, 100.04], [200, 250, 450, 250, 200]);
    picker.pick(input, output);
    assert.equal(output.size(), 1);
    assertSimilar(output.getPeak(0).mz, 100.02, "m/z");
    assertSimilar(output.getPeak(0).intensity, 450, "intensity");
  });
});

for (const sn of [1, 4]) {
  test(`PeakPickerHiRes on Orbitrap data (S/N ${sn}) reproduces the native result`, () => {
    using([new OpenMS.PeakPickerHiRes(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment(), new OpenMS.MSExperiment()],
      (picker, input, picked, expected) => {
        loadMzML("/data/PeakPickerHiRes_orbitrap.mzML", input);
        loadMzML(`/data/PeakPickerHiRes_orbitrap_sn${sn}_out.mzML`, expected);
        picker.setParameterValues({ signal_to_noise: sn });
        picker.pickExperiment(input, picked);
        assert.equal(picked.size(), expected.size());
        for (let i = 0; i < picked.size(); ++i) {
          using([picked.getSpectrum(i), expected.getSpectrum(i)], (p, e) => {
            const a = p.getPeaks(), b = e.getPeaks();
            assert.ok(a.mz.length > 0);
            assert.equal(a.mz.length, b.mz.length, `peaks in spectrum ${i}`);
            for (let k = 0; k < a.mz.length; ++k) {
              assertSimilar(a.mz[k], b.mz[k], `m/z of peak ${k} in spectrum ${i}`);
              assertSimilar(a.intensity[k], b.intensity[k], `intensity of peak ${k} in spectrum ${i}`);
            }
          });
        }
      });
  });
}

test("smoothing and filtering", () => {
  using([new OpenMS.MSSpectrum(), new OpenMS.GaussFilter(), new OpenMS.NLargest(2), new OpenMS.Normalizer()],
    (spectrum, gauss, nlargest, normalizer) => {
      const mz = Array.from({ length: 50 }, (_, i) => 500 + i * 0.01);
      const intensity = mz.map((m) => 1000 * Math.exp(-((m - 500.25) ** 2) / 0.002));
      spectrum.setPeaks(mz, intensity);
      gauss.setParameterValues({ gaussian_width: 0.05 });
      gauss.filter(spectrum);
      assert.equal(spectrum.size(), 50);
      const smoothed = spectrum.getPeaks().intensity;
      assert.ok(Math.max(...smoothed) < 1000, "smoothing lowers the apex");
      nlargest.filterPeakSpectrum(spectrum);
      assert.equal(spectrum.size(), 2);
      normalizer.filterPeakSpectrum(spectrum);
      assertClose(assert, Math.max(...spectrum.getPeaks().intensity), 1, 1e-6, "normalized to max 1");
    });
});

test("SpectrumAlignment", () => {
  using([new OpenMS.SpectrumAlignment(), new OpenMS.MSSpectrum(), new OpenMS.MSSpectrum()], (aligner, a, b) => {
    a.setPeaks([100, 200, 300], [1, 1, 1]);
    b.setPeaks([100.1, 250, 299.95], [1, 1, 1]);
    aligner.setParameterValues({ tolerance: 0.2 });
    assert.deepEqual(aligner.getSpectrumAlignment(a, b), [[0, 0], [2, 2]]);
  });
});

test("FeatureFindingMetabo reproduces the native result", () => {
  using([new OpenMS.MSExperiment(), new OpenMS.MassTraceList(), new OpenMS.MassTraceList(), new OpenMS.FeatureMap(),
         new OpenMS.MassTraceDetection(), new OpenMS.ElutionPeakDetection(), new OpenMS.FeatureFindingMetabo()],
    (input, traces, split, features, mtd, epd, ffm) => {
      loadMzML("/data/FeatureFindingMetabo_input1.mzML", input);
      mtd.run(input, traces);
      assert.ok(traces.size() > 0);
      epd.detectPeaks(traces, split);
      assert.ok(split.size() > 0); // elution peak detection splits and filters the traces
      const t = split.toArray()[0];
      assert.ok(t.mz > 0 && t.rt > 0 && t.intensity > 0);

      // FeatureFindingMetabo_test: 83 features with 13C isotope scoring, 81 with the default
      ffm.setParameterValues({ mz_scoring_13C: true });
      ffm.run(split, features);
      assert.equal(features.size(), 83);
      ffm.setParameterValues({ mz_scoring_13C: false, report_convex_hulls: true });
      ffm.run(split, features);
      assert.equal(features.size(), 81);
      const table = features.getFeatureTable();
      assert.equal(table.mz.length, 81);
      using([features.get(0)], (f) => assert.ok(f.getMassTraceBoundingBoxes().length >= 1));
    });
});
