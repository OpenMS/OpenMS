"""
Test hash functions for various OpenMS data types.
Tests that:
1. Equal objects have equal hashes
2. Objects can be used as dictionary keys
3. Objects can be used in sets
"""
import unittest

from pyopenms import (
    Peak1D, Peak2D, ChromatogramPeak, MobilityPeak1D,
    AASequence, EmpiricalFormula, PeptideHit, PeptideEvidence,
    FeatureHandle, DateTime, Adduct
)


class TestPeak1DHash(unittest.TestCase):
    """Test hash function for Peak1D"""

    def test_equal_peaks_have_equal_hash(self):
        p1 = Peak1D()
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        p2 = Peak1D()
        p2.setMZ(100.5)
        p2.setIntensity(1000.0)

        self.assertEqual(hash(p1), hash(p2))

    def test_different_peaks_have_different_hash(self):
        p1 = Peak1D()
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        p2 = Peak1D()
        p2.setMZ(200.5)
        p2.setIntensity(1000.0)

        self.assertNotEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p1 = Peak1D()
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        d = {p1: "value1"}
        self.assertEqual(d[p1], "value1")

        # Test that equal peak retrieves same value
        p2 = Peak1D()
        p2.setMZ(100.5)
        p2.setIntensity(1000.0)
        self.assertEqual(d[p2], "value1")

    def test_can_use_in_set(self):
        p1 = Peak1D()
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        p2 = Peak1D()
        p2.setMZ(100.5)
        p2.setIntensity(1000.0)

        p3 = Peak1D()
        p3.setMZ(200.5)
        p3.setIntensity(1000.0)

        s = {p1, p2, p3}
        self.assertEqual(len(s), 2)  # p1 and p2 are equal


class TestPeak2DHash(unittest.TestCase):
    """Test hash function for Peak2D"""

    def test_equal_peaks_have_equal_hash(self):
        p1 = Peak2D()
        p1.setRT(10.5)
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        p2 = Peak2D()
        p2.setRT(10.5)
        p2.setMZ(100.5)
        p2.setIntensity(1000.0)

        self.assertEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p1 = Peak2D()
        p1.setRT(10.5)
        p1.setMZ(100.5)
        p1.setIntensity(1000.0)

        d = {p1: 42}
        self.assertEqual(d[p1], 42)


class TestChromatogramPeakHash(unittest.TestCase):
    """Test hash function for ChromatogramPeak"""

    def test_equal_peaks_have_equal_hash(self):
        p1 = ChromatogramPeak()
        p1.setRT(10.5)
        p1.setIntensity(1000.0)

        p2 = ChromatogramPeak()
        p2.setRT(10.5)
        p2.setIntensity(1000.0)

        self.assertEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p1 = ChromatogramPeak()
        p1.setRT(10.5)
        p1.setIntensity(1000.0)

        d = {p1: "test"}
        self.assertEqual(d[p1], "test")


class TestMobilityPeak1DHash(unittest.TestCase):
    """Test hash function for MobilityPeak1D"""

    def test_equal_peaks_have_equal_hash(self):
        p1 = MobilityPeak1D()
        p1.setMobility(1.5)
        p1.setIntensity(1000.0)

        p2 = MobilityPeak1D()
        p2.setMobility(1.5)
        p2.setIntensity(1000.0)

        self.assertEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p1 = MobilityPeak1D()
        p1.setMobility(1.5)
        p1.setIntensity(1000.0)

        d = {p1: 99}
        self.assertEqual(d[p1], 99)


class TestAASequenceHash(unittest.TestCase):
    """Test hash function for AASequence"""

    def test_equal_sequences_have_equal_hash(self):
        seq1 = AASequence.fromString("PEPTIDE")
        seq2 = AASequence.fromString("PEPTIDE")

        self.assertEqual(hash(seq1), hash(seq2))

    def test_different_sequences_have_different_hash(self):
        seq1 = AASequence.fromString("PEPTIDE")
        seq2 = AASequence.fromString("PROTEIN")

        self.assertNotEqual(hash(seq1), hash(seq2))

    def test_can_use_in_dict(self):
        seq = AASequence.fromString("PEPTIDE")
        d = {seq: "peptide_value"}

        seq2 = AASequence.fromString("PEPTIDE")
        self.assertEqual(d[seq2], "peptide_value")

    def test_can_use_in_set(self):
        seqs = {
            AASequence.fromString("PEPTIDE"),
            AASequence.fromString("PEPTIDE"),  # duplicate
            AASequence.fromString("PROTEIN")
        }
        self.assertEqual(len(seqs), 2)


class TestEmpiricalFormulaHash(unittest.TestCase):
    """Test hash function for EmpiricalFormula"""

    def test_equal_formulas_have_equal_hash(self):
        f1 = EmpiricalFormula("H2O")
        f2 = EmpiricalFormula("H2O")

        self.assertEqual(hash(f1), hash(f2))

    def test_different_formulas_have_different_hash(self):
        f1 = EmpiricalFormula("H2O")
        f2 = EmpiricalFormula("CO2")

        self.assertNotEqual(hash(f1), hash(f2))

    def test_can_use_in_dict(self):
        f = EmpiricalFormula("H2O")
        d = {f: "water"}

        f2 = EmpiricalFormula("H2O")
        self.assertEqual(d[f2], "water")

    def test_isotope_formulas_have_different_hash(self):
        """Test that C and (13)C have different hashes"""
        f1 = EmpiricalFormula("C")
        f2 = EmpiricalFormula("(13)C")

        self.assertNotEqual(hash(f1), hash(f2))


class TestPeptideEvidenceHash(unittest.TestCase):
    """Test hash function for PeptideEvidence"""

    def test_equal_evidences_have_equal_hash(self):
        pe1 = PeptideEvidence()
        pe1.setProteinAccession("sp|P12345|TEST")
        pe1.setStart(10)
        pe1.setEnd(20)

        pe2 = PeptideEvidence()
        pe2.setProteinAccession("sp|P12345|TEST")
        pe2.setStart(10)
        pe2.setEnd(20)

        self.assertEqual(hash(pe1), hash(pe2))

    def test_can_use_in_dict(self):
        pe = PeptideEvidence()
        pe.setProteinAccession("sp|P12345|TEST")
        pe.setStart(10)
        pe.setEnd(20)

        d = {pe: "evidence"}
        self.assertEqual(d[pe], "evidence")


class TestDateTimeHash(unittest.TestCase):
    """Test hash function for DateTime"""

    def test_equal_datetimes_have_equal_hash(self):
        dt1 = DateTime()
        dt1.set("2024-12-25 10:30:00")

        dt2 = DateTime()
        dt2.set("2024-12-25 10:30:00")

        self.assertEqual(hash(dt1), hash(dt2))

    def test_can_use_in_dict(self):
        dt = DateTime()
        dt.set("2024-12-25 10:30:00")

        d = {dt: "christmas"}
        self.assertEqual(d[dt], "christmas")


class TestAdductHash(unittest.TestCase):
    """Test hash function for Adduct"""

    def test_equal_adducts_have_equal_hash(self):
        a1 = Adduct(1, 1, 22.989769, "Na", -5.0, 0.0, "")
        a2 = Adduct(1, 1, 22.989769, "Na", -5.0, 0.0, "")

        self.assertEqual(hash(a1), hash(a2))

    def test_can_use_in_dict(self):
        a = Adduct(1, 1, 22.989769, "Na", -5.0, 0.0, "")

        d = {a: "sodium"}
        self.assertEqual(d[a], "sodium")


if __name__ == "__main__":
    unittest.main()
