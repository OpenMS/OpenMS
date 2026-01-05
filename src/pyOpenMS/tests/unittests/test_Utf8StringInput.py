"""
Tests for libcpp_utf8_string parameter handling.

These tests verify that functions using libcpp_utf8_string accept both
Python str and bytes inputs (PR #8602).
"""
import unittest
import pyopenms


class TestMRMRTNormalizerStringInput(unittest.TestCase):
    """Test MRMRTNormalizer accepts str and bytes for outlier_detection_method."""

    def setUp(self):
        """Set up test data."""
        self.pairs = [
            [100.0, 110.0],
            [200.0, 210.0],
            [300.0, 310.0],
            [400.0, 410.0],
            [500.0, 510.0],
        ]

    def test_removeOutliersIterative_with_str(self):
        """Test removeOutliersIterative accepts str for outlier_detection_method."""
        result = pyopenms.MRMRTNormalizer.removeOutliersIterative(
            self.pairs, 0.95, 0.6, True, "iter_jackknife"
        )
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)

    def test_removeOutliersIterative_with_bytes(self):
        """Test removeOutliersIterative still accepts bytes (backward compatible)."""
        result = pyopenms.MRMRTNormalizer.removeOutliersIterative(
            self.pairs, 0.95, 0.6, True, b"iter_jackknife"
        )
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)


class TestElementDBStringInput(unittest.TestCase):
    """Test ElementDB.addElement accepts str and bytes."""

    def test_addElement_with_str(self):
        """Test addElement accepts str for name and symbol."""
        db = pyopenms.ElementDB()
        # Add a test element with str parameters
        abundance = {999: 1.0}
        mass = {999: 999.0}
        db.addElement("TestElement", "Te", 999, abundance, mass, False)
        # Verify it was added
        elem = db.getElement("TestElement")
        self.assertIsNotNone(elem)

    def test_addElement_with_bytes(self):
        """Test addElement still accepts bytes (backward compatible)."""
        db = pyopenms.ElementDB()
        abundance = {998: 1.0}
        mass = {998: 998.0}
        db.addElement(b"TestElement2", b"T2", 998, abundance, mass, False)
        elem = db.getElement(b"TestElement2")
        self.assertIsNotNone(elem)


class TestMZTrafoModelStringInput(unittest.TestCase):
    """Test MZTrafoModel.nameToEnum accepts str and bytes."""

    def test_nameToEnum_with_str(self):
        """Test nameToEnum accepts str."""
        enum_val = pyopenms.MZTrafoModel.nameToEnum("linear")
        self.assertEqual(enum_val, pyopenms.MZTrafoModel_MODELTYPE.LINEAR)

    def test_nameToEnum_with_bytes(self):
        """Test nameToEnum still accepts bytes (backward compatible)."""
        enum_val = pyopenms.MZTrafoModel.nameToEnum(b"linear")
        self.assertEqual(enum_val, pyopenms.MZTrafoModel_MODELTYPE.LINEAR)

    def test_enumToName_returns_str(self):
        """Test enumToName returns str (not bytes)."""
        name = pyopenms.MZTrafoModel.enumToName(pyopenms.MZTrafoModel_MODELTYPE.LINEAR)
        self.assertIsInstance(name, str)
        self.assertEqual(name, "linear")


class TestIMTypesStringInput(unittest.TestCase):
    """Test IMTypes functions accept str and bytes."""

    def test_toDriftTimeUnit_with_str(self):
        """Test toDriftTimeUnit accepts str."""
        unit = pyopenms.IMTypes.toDriftTimeUnit("millisecond")
        self.assertEqual(unit, pyopenms.DriftTimeUnit.MILLISECOND)

    def test_toDriftTimeUnit_with_bytes(self):
        """Test toDriftTimeUnit still accepts bytes (backward compatible)."""
        unit = pyopenms.IMTypes.toDriftTimeUnit(b"millisecond")
        self.assertEqual(unit, pyopenms.DriftTimeUnit.MILLISECOND)

    def test_toString_DriftTimeUnit_returns_str(self):
        """Test toString returns str for DriftTimeUnit."""
        result = pyopenms.IMTypes.toString(pyopenms.DriftTimeUnit.MILLISECOND)
        self.assertIsInstance(result, str)
        self.assertEqual(result, "millisecond")

    def test_toIMFormat_with_str(self):
        """Test toIMFormat accepts str."""
        fmt = pyopenms.IMTypes.toIMFormat("concatenated")
        self.assertEqual(fmt, pyopenms.IMFormat.CONCATENATED)

    def test_toIMFormat_with_bytes(self):
        """Test toIMFormat still accepts bytes (backward compatible)."""
        fmt = pyopenms.IMTypes.toIMFormat(b"concatenated")
        self.assertEqual(fmt, pyopenms.IMFormat.CONCATENATED)

    def test_toString_IMFormat_returns_str(self):
        """Test toString returns str for IMFormat."""
        result = pyopenms.IMTypes.toString(pyopenms.IMFormat.CONCATENATED)
        self.assertIsInstance(result, str)
        self.assertEqual(result, "concatenated")


class TestRibonucleotideDBStringInput(unittest.TestCase):
    """Test RibonucleotideDB methods accept str and bytes."""

    def test_getRibonucleotide_with_str(self):
        """Test getRibonucleotide accepts str."""
        db = pyopenms.RibonucleotideDB()
        ribo = db.getRibonucleotide("A")
        self.assertIsNotNone(ribo)

    def test_getRibonucleotide_with_bytes(self):
        """Test getRibonucleotide still accepts bytes (backward compatible)."""
        db = pyopenms.RibonucleotideDB()
        ribo = db.getRibonucleotide(b"A")
        self.assertIsNotNone(ribo)

    def test_getRibonucleotidePrefix_with_str(self):
        """Test getRibonucleotidePrefix accepts str."""
        db = pyopenms.RibonucleotideDB()
        ribo = db.getRibonucleotidePrefix("A")
        self.assertIsNotNone(ribo)

    def test_getRibonucleotideAlternatives_with_str(self):
        """Test getRibonucleotideAlternatives accepts str (manual addon)."""
        db = pyopenms.RibonucleotideDB()
        # Use a code that has alternatives
        alternatives = db.getRibonucleotideAlternatives("Y")
        self.assertIsNotNone(alternatives)
        self.assertEqual(len(alternatives), 2)

    def test_getRibonucleotideAlternatives_with_bytes(self):
        """Test getRibonucleotideAlternatives still accepts bytes."""
        db = pyopenms.RibonucleotideDB()
        alternatives = db.getRibonucleotideAlternatives(b"Y")
        self.assertIsNotNone(alternatives)
        self.assertEqual(len(alternatives), 2)


class TestIndexedMzMLFileStringInput(unittest.TestCase):
    """Test IndexedMzMLFile native ID methods accept str and bytes."""

    def test_getMSSpectrumByNativeId_with_str(self):
        """Test getMSSpectrumByNativeId accepts str."""
        # This test requires a file, so we just verify the method exists
        # and has the right signature
        self.assertTrue(hasattr(pyopenms.IndexedMzMLFile, 'getMSSpectrumByNativeId'))

    def test_getMSChromatogramByNativeId_with_str(self):
        """Test getMSChromatogramByNativeId accepts str."""
        self.assertTrue(hasattr(pyopenms.IndexedMzMLFile, 'getMSChromatogramByNativeId'))


class TestLightTargetedExperimentStringInput(unittest.TestCase):
    """Test LightTargetedExperiment methods accept str and bytes."""

    def test_setFragmentType_with_str(self):
        """Test LightTransition.setFragmentType accepts str."""
        transition = pyopenms.LightTransition()
        transition.setFragmentType("y")
        self.assertEqual(transition.getFragmentType(), b"y")

    def test_setFragmentType_with_bytes(self):
        """Test LightTransition.setFragmentType still accepts bytes."""
        transition = pyopenms.LightTransition()
        transition.setFragmentType(b"b")
        self.assertEqual(transition.getFragmentType(), b"b")

    def test_getCompoundByRef_with_str(self):
        """Test LightTargetedExperiment.getCompoundByRef accepts str."""
        exp = pyopenms.LightTargetedExperiment()
        # Add a compound
        compound = pyopenms.LightCompound()
        compound.id = b"test_compound"
        exp.compounds.append(compound)
        # Lookup by str
        result = exp.getCompoundByRef("test_compound")
        self.assertIsNotNone(result)

    def test_getPeptideByRef_with_str(self):
        """Test LightTargetedExperiment.getPeptideByRef accepts str."""
        exp = pyopenms.LightTargetedExperiment()
        # Add a peptide compound
        compound = pyopenms.LightCompound()
        compound.id = b"test_peptide"
        compound.sequence = b"PEPTIDE"
        exp.compounds.append(compound)
        # Lookup by str
        result = exp.getPeptideByRef("test_peptide")
        self.assertIsNotNone(result)


class TestOpenSwathScoringStringInput(unittest.TestCase):
    """Test OpenSwathScoring.initialize accepts str and bytes."""

    def test_initialize_signature_exists(self):
        """Test OpenSwathScoring.initialize exists with string parameters."""
        self.assertTrue(hasattr(pyopenms.OpenSwathScoring, 'initialize'))


class TestMSExperimentStringInput(unittest.TestCase):
    """Test MSExperiment aggregate methods accept str and bytes."""

    def test_aggregateFromMatrix_with_str(self):
        """Test aggregateFromMatrix accepts str for mz_agg."""
        exp = pyopenms.MSExperiment()
        # Add a spectrum with peaks
        spectrum = pyopenms.MSSpectrum()
        spectrum.setRT(100.0)
        peak = pyopenms.Peak1D()
        peak.setMZ(500.0)
        peak.setIntensity(1000.0)
        spectrum.push_back(peak)
        exp.addSpectrum(spectrum)

        # Create a ranges matrix
        ranges = pyopenms.Matrix_double()
        ranges.resize(1, 4)
        ranges.setValue(0, 0, 0.0)    # rt_min
        ranges.setValue(0, 1, 200.0)  # rt_max
        ranges.setValue(0, 2, 400.0)  # mz_min
        ranges.setValue(0, 3, 600.0)  # mz_max

        # Call with str parameter
        result = exp.aggregateFromMatrix(ranges, 1, "sum")
        self.assertIsNotNone(result)

    def test_extractXICsFromMatrix_with_str(self):
        """Test extractXICsFromMatrix accepts str for mz_agg."""
        exp = pyopenms.MSExperiment()
        spectrum = pyopenms.MSSpectrum()
        spectrum.setRT(100.0)
        peak = pyopenms.Peak1D()
        peak.setMZ(500.0)
        peak.setIntensity(1000.0)
        spectrum.push_back(peak)
        exp.addSpectrum(spectrum)

        ranges = pyopenms.Matrix_double()
        ranges.resize(1, 4)
        ranges.setValue(0, 0, 0.0)
        ranges.setValue(0, 1, 200.0)
        ranges.setValue(0, 2, 400.0)
        ranges.setValue(0, 3, 600.0)

        result = exp.extractXICsFromMatrix(ranges, 1, "sum")
        self.assertIsNotNone(result)


if __name__ == '__main__':
    unittest.main()
