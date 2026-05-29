"""
Tests for newly bound CHEMISTRY methods in the nanobind build:

- EmpiricalFormula.getNumberOf(Element)
- EmpiricalFormula.getLightestIsotopeWeight()
- Residue.getHydrophobicity(HydrophobicityScaleMethod) + the enum
- TheoreticalSpectrumGenerator.getPrefixAndSuffixIonsMZ(peptide, charge)
- ModificationsDB.searchModification(ResidueModification)
"""
import math
import numpy as np
import pyopenms as poms


def test_empirical_formula_get_number_of():
    ef = poms.EmpiricalFormula("C6H12O6")
    db = poms.ElementDB.getInstance()
    carbon = db.getElement("Carbon")
    hydrogen = db.getElement("Hydrogen")
    oxygen = db.getElement("Oxygen")
    assert ef.getNumberOf(carbon) == 6
    assert ef.getNumberOf(hydrogen) == 12
    assert ef.getNumberOf(oxygen) == 6


def test_empirical_formula_lightest_isotope_weight():
    ef = poms.EmpiricalFormula("H2O")
    w = ef.getLightestIsotopeWeight()
    assert isinstance(w, float)
    assert w > 0.0
    # the lightest isotope weight must not exceed the monoisotopic weight
    assert w <= ef.getMonoWeight() + 1e-6


def test_hydrophobicity_scale_method_enum():
    # enum is exposed at module scope and convertible to int
    assert int(poms.HydrophobicityScaleMethod.KYTE_DOOLITTLE) == 0
    assert int(poms.HydrophobicityScaleMethod.EISENBERG) == 1
    assert int(poms.HydrophobicityScaleMethod.EISENBERG_CONSENSUS) == 6


def test_residue_get_hydrophobicity():
    # Leucine on the Kyte-Doolittle scale is the documented value 3.8
    leu = poms.ResidueDB.getInstance().getResidue("L")
    val = leu.getHydrophobicity(poms.HydrophobicityScaleMethod.KYTE_DOOLITTLE)
    assert math.isclose(val, 3.8, abs_tol=1e-6)


def test_tsg_prefix_and_suffix_ions_mz():
    tsg = poms.TheoreticalSpectrumGenerator()
    seq = poms.AASequence.fromString("PEPTIDE")
    out = tsg.getPrefixAndSuffixIonsMZ(seq, 1)
    assert isinstance(out, np.ndarray)
    assert out.ndim == 1
    assert out.size > 0
    assert np.all(out > 0.0)


def test_modificationsdb_search_modification():
    db = poms.ModificationsDB.getInstance()
    mod = db.getModification("Oxidation")
    assert mod is not None
    found = db.searchModification(mod)
    assert found is not None
    assert found.getId() == mod.getId()
