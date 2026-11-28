// Sequences, formulas, isotopes, modifications, digestion, theoretical spectra.
// Reference values are those of the corresponding OpenMS class tests.

import { test } from "node:test";
import assert from "node:assert/strict";
import { loadOpenMS, using, assertClose } from "./helpers.mjs";

const OpenMS = await loadOpenMS();

test("version and shared data", () => {
  assert.match(OpenMS.getVersion(), /^\d+\.\d+\.\d+/);
  assert.equal(OpenMS.getOpenMSDataPath(), "/share/OpenMS");
});

test("AASequence masses", () => {
  using([OpenMS.AASequence.fromString("DFPIANGER")], (seq) => {
    assert.equal(seq.size(), 9);
    assert.equal(seq.toString(), "DFPIANGER");
    assert.equal(seq.isModified(), false);
    // AASequence_test: getMonoWeight of "DFPIANGER" (full, charge 1) = 1018.4952
    assertClose(assert, seq.getMonoWeight(OpenMS.ResidueType.Full, 1), 1018.4952, 1e-3, "mono weight +1");
    assertClose(assert, seq.getMonoWeight(), 1017.4879, 1e-3, "neutral mono weight");
    assertClose(assert, seq.getMZ(2), (1017.4879 + 2 * 1.007276) / 2, 1e-3, "m/z 2+");
    using([seq.getFormula()], (f) => assert.equal(f.toString(), "C44H67N13O15"));
    using([seq.getPrefix(3), seq.getSuffix(2)], (prefix, suffix) => {
      assert.equal(prefix.toString(), "DFP");
      assert.equal(suffix.toString(), "ER");
    });
    assert.equal(seq.getResidue(0).oneLetterCode, "D");
    assert.deepEqual(Object.keys(seq.getAAFrequencies()).sort(), ["A", "D", "E", "F", "G", "I", "N", "P", "R"]);
  });
});

test("AASequence modifications", () => {
  using([OpenMS.AASequence.fromString("PEPTM(Oxidation)IDEK")], (seq) => {
    assert.equal(seq.isModified(), true);
    assert.equal(seq.toUnmodifiedString(), "PEPTMIDEK");
    assert.equal(seq.getResidue(4).modificationName, "Oxidation");
    using([OpenMS.AASequence.fromString("PEPTMIDEK")], (unmodified) => {
      assertClose(assert, seq.getMonoWeight() - unmodified.getMonoWeight(), 15.9949, 1e-3, "oxidation mass shift");
    });
    seq.setNTerminalModification("Acetyl");
    assert.equal(seq.hasNTerminalModification(), true);
    assert.equal(seq.toString(), ".(Acetyl)PEPTM(Oxidation)IDEK");
  });
});

test("C++ exceptions become JavaScript errors", () => {
  // e.g. "Parse Error: Cannot convert string to peptide modification: missing ']' in: PEPT[X"
  assert.throws(() => OpenMS.AASequence.fromString("PEPT[X"), (e) => e instanceof Error && /Parse ?Error/.test(e.message) && e.message.includes("PEPT[X"));
  assert.throws(() => OpenMS.EmpiricalFormula.fromString("C6H12O6Xx"), /Error/);
  assert.throws(() => OpenMS.ModificationsDB.getModification("NoSuchModification"), /InvalidValue: the value 'NoSuchModification'/);
});

test("EmpiricalFormula and isotope distributions", () => {
  using([OpenMS.EmpiricalFormula.fromString("C6H12O6")], (glucose) => {
    assertClose(assert, glucose.getMonoWeight(), 180.06339, 1e-4, "glucose mono weight");
    assertClose(assert, glucose.getAverageWeight(), 180.156, 1e-2, "glucose average weight");
    assert.deepEqual(glucose.toObject(), { C: 6, H: 12, O: 6 });
    assert.equal(glucose.getNumberOfAtoms(), 24);

    using([new OpenMS.CoarseIsotopePatternGenerator(3)], (generator) => {
      using([glucose.getIsotopeDistribution(generator)], (dist) => {
        assert.equal(dist.size(), 3);
        const { mz, intensity } = dist.getPeaks();
        assert.ok(mz instanceof Float64Array);
        assertClose(assert, mz[0], 180.06339, 1e-3, "monoisotopic peak");
        assert.ok(intensity[0] > intensity[1] && intensity[1] > intensity[2]);
        assertClose(assert, intensity.reduce((a, b) => a + b, 0), 1, 1e-6, "normalized");
      });
    });

    using([new OpenMS.FineIsotopePatternGenerator(1e-3)], (generator) => {
      using([generator.run(glucose)], (dist) => {
        assert.ok(dist.size() > 3, "fine structure resolves more isotopologues");
        assertClose(assert, dist.getMostAbundant().mz, 180.06339, 1e-3, "most abundant");
      });
    });
  });
});

test("ModificationsDB and ResidueDB", () => {
  const ox = OpenMS.ModificationsDB.getModification("Oxidation", "M", "any");
  assert.equal(ox.fullId, "Oxidation (M)");
  assertClose(assert, ox.diffMonoMass, 15.994915, 1e-5, "oxidation");
  assert.equal(ox.unimodAccession, "UniMod:35");
  assert.ok(OpenMS.ModificationsDB.getAllSearchModifications().includes("Phospho (S)"));
  assert.ok(OpenMS.ModificationsDB.searchModificationsByDiffMonoMass(79.966, 0.01, "S").includes("Phospho (S)"));

  const lys = OpenMS.ResidueDB.getResidue("K");
  assert.equal(lys.threeLetterCode, "Lys");
  assertClose(assert, lys.monoWeight, 146.1055, 1e-3, "lysine");
  assert.equal(OpenMS.ResidueDB.getResidueNames().length, 20);
});

test("ProteaseDigestion", () => {
  assert.ok(OpenMS.ProteaseDB.getAllNames().includes("Trypsin"));
  using([new OpenMS.ProteaseDigestion(), OpenMS.AASequence.fromString("ACDKDEFKRGHIKPLMNR")], (digestion, protein) => {
    digestion.setEnzyme("Trypsin");
    assert.equal(digestion.getEnzymeName(), "Trypsin");
    digestion.setMissedCleavages(0);
    // Trypsin cleaves after K/R, not before P (KP)
    assert.deepEqual(digestion.digestUnmodified("ACDKDEFKRGHIKPLMNR", 1, 0), ["ACDK", "DEFK", "R", "GHIKPLMNR"]);
    const peptides = digestion.digest(protein, 2, 0);
    try {
      assert.deepEqual(peptides.map((p) => p.toString()), ["ACDK", "DEFK", "GHIKPLMNR"]);
    } finally {
      peptides.forEach((p) => p.delete());
    }
    digestion.setMissedCleavages(1);
    assert.equal(digestion.digestUnmodified("ACDKDEFKRGHIKPLMNR", 1, 0).length, 7);
  });
});

test("TheoreticalSpectrumGenerator", () => {
  using([new OpenMS.TheoreticalSpectrumGenerator(), OpenMS.AASequence.fromString("IFSQVGK"), new OpenMS.MSSpectrum()],
    (generator, peptide, spectrum) => {
      // TheoreticalSpectrumGenerator_test: b and y ions of IFSQVGK, charge 1 (no b1 by default)
      generator.getSpectrum(spectrum, peptide, 1, 1);
      const expected = [147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362];
      const { mz } = spectrum.getPeaks();
      assert.equal(mz.length, expected.length);
      expected.forEach((e, i) => assertClose(assert, mz[i], e, 1e-3, `peak ${i}`));
      assert.equal(spectrum.getMSLevel(), 2);
      assert.equal(spectrum.getType(), "Centroid");
      assertClose(assert, spectrum.getPrecursors()[0].mz, peptide.getMZ(2, OpenMS.ResidueType.Full), 1e-6, "precursor m/z");

      spectrum.clear(true);
      generator.getSpectrum(spectrum, peptide, 1, 2);
      assert.equal(spectrum.size(), 22);

      generator.setParameterValues({ add_first_prefix_ion: true, add_metainfo: true });
      assert.equal(generator.getParameterValues().add_first_prefix_ion, true);
      spectrum.clear(true);
      generator.getSpectrum(spectrum, peptide, 1, 1);
      assert.equal(spectrum.size(), 12);
      assertClose(assert, spectrum.getPeak(0).mz, 114.091, 1e-3, "b1");
      const names = spectrum.getStringDataArrays().find((a) => a.name === "IonNames");
      assert.ok(names.data.includes("y1+"));
    });
});
