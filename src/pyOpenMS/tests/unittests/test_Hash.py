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
    FeatureHandle, DateTime, Adduct, DataValue, MetaInfoInterface,
    CVTerm, Software, Precursor, ProteinHit, PeptideIdentification,
    Product, SpectrumSettings, ChromatogramSettings
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


class TestDataValueHash(unittest.TestCase):
    """Test hash function for DataValue - critical epsilon test"""

    def test_equal_int_values_have_equal_hash(self):
        dv1 = DataValue(42)
        dv2 = DataValue(42)
        self.assertEqual(hash(dv1), hash(dv2))

    def test_equal_string_values_have_equal_hash(self):
        dv1 = DataValue("test")
        dv2 = DataValue("test")
        self.assertEqual(hash(dv1), hash(dv2))

    def test_double_epsilon_compatibility(self):
        """Critical test: doubles within epsilon (1e-6) must have equal hashes"""
        dv1 = DataValue(1.0000001)
        dv2 = DataValue(1.0000002)  # Within epsilon
        # These are equal per operator== so must have equal hash
        self.assertEqual(hash(dv1), hash(dv2))

    def test_different_doubles_have_different_hash(self):
        dv1 = DataValue(1.0)
        dv2 = DataValue(2.0)
        self.assertNotEqual(hash(dv1), hash(dv2))

    def test_can_use_in_dict(self):
        dv = DataValue(42)
        d = {dv: "answer"}
        self.assertEqual(d[dv], "answer")


class TestMetaInfoInterfaceHash(unittest.TestCase):
    """Test hash function for MetaInfoInterface - critical order-independence test"""

    def test_equal_meta_have_equal_hash(self):
        m1 = MetaInfoInterface()
        m1.setMetaValue("key1", "value1")

        m2 = MetaInfoInterface()
        m2.setMetaValue("key1", "value1")

        self.assertEqual(hash(m1), hash(m2))

    def test_order_independent_hash(self):
        """Critical test: insertion order should not affect hash"""
        m1 = MetaInfoInterface()
        m1.setMetaValue("name", "test")
        m1.setMetaValue("score", 1.5)

        m2 = MetaInfoInterface()
        m2.setMetaValue("score", 1.5)  # Different order
        m2.setMetaValue("name", "test")

        self.assertEqual(hash(m1), hash(m2))

    def test_can_use_in_dict(self):
        m = MetaInfoInterface()
        m.setMetaValue("key", "value")
        d = {m: "meta"}
        self.assertEqual(d[m], "meta")


class TestPeptideHitHash(unittest.TestCase):
    """Test hash function for PeptideHit"""

    def test_equal_hits_have_equal_hash(self):
        ph1 = PeptideHit()
        ph1.setSequence(AASequence.fromString("PEPTIDE"))
        ph1.setScore(10.5)
        ph1.setCharge(2)

        ph2 = PeptideHit()
        ph2.setSequence(AASequence.fromString("PEPTIDE"))
        ph2.setScore(10.5)
        ph2.setCharge(2)

        self.assertEqual(hash(ph1), hash(ph2))

    def test_different_hits_have_different_hash(self):
        ph1 = PeptideHit()
        ph1.setSequence(AASequence.fromString("PEPTIDE"))

        ph2 = PeptideHit()
        ph2.setSequence(AASequence.fromString("PROTEIN"))

        self.assertNotEqual(hash(ph1), hash(ph2))

    def test_can_use_in_set(self):
        hits = set()
        for seq in ["PEPTIDE", "PROTEIN", "PEPTIDE"]:  # PEPTIDE twice
            ph = PeptideHit()
            ph.setSequence(AASequence.fromString(seq))
            hits.add(ph)
        self.assertEqual(len(hits), 2)


class TestCVTermHash(unittest.TestCase):
    """Test hash function for CVTerm"""

    def test_equal_terms_have_equal_hash(self):
        cv1 = CVTerm()
        cv1.setAccession("MS:1000001")
        cv1.setName("sample number")

        cv2 = CVTerm()
        cv2.setAccession("MS:1000001")
        cv2.setName("sample number")

        self.assertEqual(hash(cv1), hash(cv2))

    def test_can_use_in_dict(self):
        cv = CVTerm()
        cv.setAccession("MS:1000001")
        d = {cv: "term"}
        self.assertEqual(d[cv], "term")


class TestSoftwareHash(unittest.TestCase):
    """Test hash function for Software"""

    def test_equal_software_have_equal_hash(self):
        sw1 = Software()
        sw1.setName("OpenMS")
        sw1.setVersion("3.0")

        sw2 = Software()
        sw2.setName("OpenMS")
        sw2.setVersion("3.0")

        self.assertEqual(hash(sw1), hash(sw2))

    def test_can_use_in_dict(self):
        sw = Software()
        sw.setName("OpenMS")
        d = {sw: "software"}
        self.assertEqual(d[sw], "software")


class TestPrecursorHash(unittest.TestCase):
    """Test hash function for Precursor"""

    def test_equal_precursors_have_equal_hash(self):
        p1 = Precursor()
        p1.setMZ(500.0)
        p1.setCharge(2)

        p2 = Precursor()
        p2.setMZ(500.0)
        p2.setCharge(2)

        self.assertEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p = Precursor()
        p.setMZ(500.0)
        d = {p: "precursor"}
        self.assertEqual(d[p], "precursor")


class TestProteinHitHash(unittest.TestCase):
    """Test hash function for ProteinHit"""

    def test_equal_hits_have_equal_hash(self):
        ph1 = ProteinHit()
        ph1.setAccession("sp|P12345|TEST")
        ph1.setScore(100.0)

        ph2 = ProteinHit()
        ph2.setAccession("sp|P12345|TEST")
        ph2.setScore(100.0)

        self.assertEqual(hash(ph1), hash(ph2))

    def test_can_use_in_set(self):
        hits = set()
        for acc in ["P12345", "P67890", "P12345"]:
            ph = ProteinHit()
            ph.setAccession(acc)
            hits.add(ph)
        self.assertEqual(len(hits), 2)


class TestPeptideIdentificationHash(unittest.TestCase):
    """Test hash function for PeptideIdentification"""

    def test_equal_ids_have_equal_hash(self):
        pi1 = PeptideIdentification()
        pi1.setIdentifier("test_run")
        pi1.setScoreType("XCorr")

        pi2 = PeptideIdentification()
        pi2.setIdentifier("test_run")
        pi2.setScoreType("XCorr")

        self.assertEqual(hash(pi1), hash(pi2))

    def test_can_use_in_dict(self):
        pi = PeptideIdentification()
        pi.setIdentifier("test")
        d = {pi: "identification"}
        self.assertEqual(d[pi], "identification")


class TestProductHash(unittest.TestCase):
    """Test hash function for Product"""

    def test_equal_products_have_equal_hash(self):
        p1 = Product()
        p1.setMZ(500.0)

        p2 = Product()
        p2.setMZ(500.0)

        self.assertEqual(hash(p1), hash(p2))

    def test_can_use_in_dict(self):
        p = Product()
        p.setMZ(500.0)
        d = {p: "product"}
        self.assertEqual(d[p], "product")


if __name__ == "__main__":
    unittest.main()
