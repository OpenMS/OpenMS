"""
Tests for newly bound CHEMISTRY methods in the nanobind build:

- EmpiricalFormula.getNumberOf(Element)
- EmpiricalFormula.getLightestIsotopeWeight()
- Residue.getHydrophobicity(HydrophobicityScaleMethod) + the enum
- TheoreticalSpectrumGenerator.getPrefixAndSuffixIonsMZ(peptide, charge)
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


def test_proteomics_pka_scale_enum():
    # enum is exposed at module scope and convertible to int
    assert int(poms.ProteomicsPkaScale.LEHNINGER) == 0
    assert int(poms.ProteomicsPkaScale.EMBOSS) == 1
    assert int(poms.ProteomicsPkaScale.SILLERO) == 2
    assert int(poms.ProteomicsPkaScale.BJELLQVIST) == 3


def test_isoelectric_point_compute_charge():
    seq = poms.AASequence.fromString("PEPTIDE")
    # At low pH, peptide should be positively charged
    charge_low = poms.IsoelectricPoint.computeCharge(seq, 1.0)
    assert charge_low > 0.0
    # At high pH, peptide should be negatively charged
    charge_high = poms.IsoelectricPoint.computeCharge(seq, 14.0)
    assert charge_high < 0.0


def test_isoelectric_point_compute_pi():
    seq = poms.AASequence.fromString("PEPTIDE")
    pi = poms.IsoelectricPoint.computePI(seq)
    # pI should be between 0 and 14
    assert 0.0 < pi < 14.0
    # At the pI, the charge should be approximately zero
    charge_at_pi = poms.IsoelectricPoint.computeCharge(seq, pi)
    assert math.isclose(charge_at_pi, 0.0, abs_tol=0.01)


def test_isoelectric_point_alanine():
    ala = poms.AASequence.fromString("A")
    pi = poms.IsoelectricPoint.computePI(ala)
    # Lehninger: (9.69 + 2.34)/2 = 6.015
    assert math.isclose(pi, 6.015, abs_tol=0.01)


def test_isoelectric_point_emboss_scale():
    seq = poms.AASequence.fromString("ACDEFGHIKLMNPQRSTVWY")
    pi_leh = poms.IsoelectricPoint.computePI(seq, poms.ProteomicsPkaScale.LEHNINGER)
    pi_emb = poms.IsoelectricPoint.computePI(seq, poms.ProteomicsPkaScale.EMBOSS)
    # Different scales should yield different pI values
    assert not math.isclose(pi_leh, pi_emb, abs_tol=0.001)


def test_isoelectric_point_bjellqvist_scale():
    # Bjellqvist: pI(A) uses N-term pKa 7.59 (A-specific), C-term pKa 3.55 => pI ~= 5.57
    ala = poms.AASequence.fromString("A")
    pi_bjell = poms.IsoelectricPoint.computePI(ala, poms.ProteomicsPkaScale.BJELLQVIST)
    assert math.isclose(pi_bjell, 5.57, abs_tol=0.01)
    # Bjellqvist pI for A must differ from Lehninger (6.015)
    pi_leh = poms.IsoelectricPoint.computePI(ala, poms.ProteomicsPkaScale.LEHNINGER)
    assert not math.isclose(pi_bjell, pi_leh, abs_tol=0.1)
    # Proline N-term pKa 8.36 => pI(P) ~= (8.36 + 3.55) / 2 = 5.955
    pro = poms.AASequence.fromString("P")
    pi_pro = poms.IsoelectricPoint.computePI(pro, poms.ProteomicsPkaScale.BJELLQVIST)
    assert math.isclose(pi_pro, 5.955, abs_tol=0.01)
