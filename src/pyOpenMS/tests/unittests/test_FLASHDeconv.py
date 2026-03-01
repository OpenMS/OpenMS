"""
Test suite for FLASHDeconv Python bindings.

Tests cover:
- FLASHHelperClasses (LogMzPeak, PrecalculatedAveragine, MassFeature, IsobaricQuantities)
- PeakGroup
- DeconvolvedSpectrum
- SpectralDeconvolution
- FLASHDeconvAlgorithm
- FLASHDeconvSpectrumFile
- FLASHDeconvFeatureFile
"""

import unittest
import os
import tempfile

import pyopenms


class TestFLASHHelperClasses(unittest.TestCase):
    """Test FLASHHelperClasses wrapper class and static functions."""

    def test_class_instantiation(self):
        """Test that FLASHHelperClasses can be instantiated."""
        helper = pyopenms.FLASHHelperClasses()
        self.assertIsNotNone(helper)

    def test_getLogMz_positive(self):
        """Test getLogMz for positive ionization mode."""
        mz = 1300.0
        log_mz = pyopenms.FLASHHelperClasses.getLogMz(mz, True)
        self.assertAlmostEqual(log_mz, 7.169344415063863, places=4)

    def test_getLogMz_negative(self):
        """Test getLogMz for negative ionization mode."""
        mz = 1300.0
        log_mz = pyopenms.FLASHHelperClasses.getLogMz(mz, False)
        self.assertAlmostEqual(log_mz, 7.170894071437545, places=4)

    def test_getLogMz_different_modes(self):
        """Test that positive and negative modes give different results."""
        mz = 1300.0
        log_mz_pos = pyopenms.FLASHHelperClasses.getLogMz(mz, True)
        log_mz_neg = pyopenms.FLASHHelperClasses.getLogMz(mz, False)
        self.assertNotEqual(log_mz_pos, log_mz_neg)

    def test_getChargeMass_positive(self):
        """Test getChargeMass for positive ionization mode."""
        charge_mass = pyopenms.FLASHHelperClasses.getChargeMass(True)
        # Should be approximately proton mass
        self.assertAlmostEqual(charge_mass, 1.007276, places=4)

    def test_getChargeMass_negative(self):
        """Test getChargeMass for negative ionization mode."""
        charge_mass = pyopenms.FLASHHelperClasses.getChargeMass(False)
        # Should be approximately negative proton mass
        self.assertAlmostEqual(charge_mass, -1.007276, places=4)


class TestLogMzPeak(unittest.TestCase):
    """Test LogMzPeak struct."""

    def test_default_constructor(self):
        """Test default constructor."""
        peak = pyopenms.LogMzPeak()
        self.assertIsNotNone(peak)

    def test_constructor_from_peak1d_positive(self):
        """Test constructor from Peak1D in positive mode."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        log_peak = pyopenms.LogMzPeak(p, True)

        self.assertAlmostEqual(log_peak.mz, 1125.5118055019082, places=4)
        self.assertAlmostEqual(log_peak.intensity, 443505.625, places=1)
        self.assertAlmostEqual(log_peak.logMz, 7.0250977989903145, places=4)
        self.assertTrue(log_peak.is_positive)
        self.assertEqual(log_peak.abs_charge, 0)
        self.assertEqual(log_peak.isotopeIndex, 0)

    def test_constructor_from_peak1d_negative(self):
        """Test constructor from Peak1D in negative mode."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        log_peak_pos = pyopenms.LogMzPeak(p, True)
        log_peak_neg = pyopenms.LogMzPeak(p, False)

        self.assertFalse(log_peak_neg.is_positive)
        # logMz should be different for negative mode
        self.assertNotEqual(log_peak_neg.logMz, log_peak_pos.logMz)

    def test_copy_constructor(self):
        """Test copy constructor."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        original = pyopenms.LogMzPeak(p, True)
        copy = pyopenms.LogMzPeak(original)

        self.assertAlmostEqual(copy.mz, original.mz, places=4)
        self.assertAlmostEqual(copy.intensity, original.intensity, places=1)
        self.assertAlmostEqual(copy.logMz, original.logMz, places=6)

    def test_getUnchargedMass_with_charge(self):
        """Test getUnchargedMass when charge is set."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        log_peak = pyopenms.LogMzPeak(p, True)
        log_peak.abs_charge = 2

        mass = log_peak.getUnchargedMass()
        self.assertAlmostEqual(mass, 2249.0090580702745, places=2)

    def test_getUnchargedMass_zero_charge(self):
        """Test getUnchargedMass when charge is 0."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        log_peak = pyopenms.LogMzPeak(p, True)
        log_peak.abs_charge = 0

        mass = log_peak.getUnchargedMass()
        self.assertAlmostEqual(mass, 0.0, places=4)

    def test_getUnchargedMass_preset_mass(self):
        """Test getUnchargedMass when mass is already set."""
        p = pyopenms.Peak1D()
        p.setMZ(1125.5118055019082)
        p.setIntensity(443505.625)

        log_peak = pyopenms.LogMzPeak(p, True)
        log_peak.abs_charge = 2
        log_peak.mass = 5000.0

        mass = log_peak.getUnchargedMass()
        self.assertAlmostEqual(mass, 5000.0, places=4)

    def test_comparison_operators(self):
        """Test comparison operators."""
        p = pyopenms.Peak1D()
        p.setMZ(1000.0)
        p.setIntensity(100.0)

        peak1 = pyopenms.LogMzPeak(p, True)
        peak2 = pyopenms.LogMzPeak(p, True)
        peak2.logMz = 8.0  # Make peak2 have larger logMz

        self.assertLess(peak1, peak2)
        self.assertTrue(peak2 > peak1)

    def test_equality_operator(self):
        """Test equality operator."""
        p = pyopenms.Peak1D()
        p.setMZ(1000.0)
        p.setIntensity(100.0)

        peak1 = pyopenms.LogMzPeak(p, True)
        peak2 = pyopenms.LogMzPeak(peak1)

        self.assertTrue(peak1 == peak2)

    def test_member_variables(self):
        """Test that member variables can be accessed and modified."""
        peak = pyopenms.LogMzPeak()
        peak.mz = 500.0
        peak.intensity = 1000.0
        peak.logMz = 6.0
        peak.mass = 2000.0
        peak.abs_charge = 3
        peak.is_positive = True
        peak.isotopeIndex = 2

        self.assertAlmostEqual(peak.mz, 500.0, places=4)
        self.assertAlmostEqual(peak.intensity, 1000.0, places=4)
        self.assertAlmostEqual(peak.logMz, 6.0, places=4)
        self.assertAlmostEqual(peak.mass, 2000.0, places=4)
        self.assertEqual(peak.abs_charge, 3)
        self.assertTrue(peak.is_positive)
        self.assertEqual(peak.isotopeIndex, 2)


class TestPrecalculatedAveragine(unittest.TestCase):
    """Test PrecalculatedAveragine (PrecalAveragine) class."""

    def test_default_constructor(self):
        """Test default constructor."""
        avg = pyopenms.PrecalAveragine()
        self.assertIsNotNone(avg)

    def test_constructor_with_generator(self):
        """Test constructor with CoarseIsotopePatternGenerator."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        apex_idx = avg.getApexIndex(75.0)
        mass_delta = avg.getAverageMassDelta(75.0)

        self.assertEqual(apex_idx, 0)
        self.assertAlmostEqual(mass_delta, 0.04, delta=0.3)

    def test_constructor_with_decoy_distance(self):
        """Test constructor with decoy isotope distance."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False, 1.5)

        iso = avg.get(75.0)
        self.assertGreater(iso.size(), 0)

    def test_get_isotope_distribution(self):
        """Test get method for isotope distribution."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        iso = avg.get(60.0)
        self.assertIsNotNone(iso)
        self.assertGreater(iso.size(), 0)

    def test_max_isotope_index(self):
        """Test getMaxIsotopeIndex and setMaxIsotopeIndex."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        avg.setMaxIsotopeIndex(4)
        self.assertEqual(avg.getMaxIsotopeIndex(), 4)

    def test_left_count_from_apex(self):
        """Test getLeftCountFromApex."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        left_count = avg.getLeftCountFromApex(75.0)
        self.assertIsInstance(left_count, int)

    def test_right_count_from_apex(self):
        """Test getRightCountFromApex."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        right_count = avg.getRightCountFromApex(75.0)
        self.assertIsInstance(right_count, int)

    def test_apex_index(self):
        """Test getApexIndex."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        apex_idx = avg.getApexIndex(75.0)
        self.assertEqual(apex_idx, 0)

    def test_last_index(self):
        """Test getLastIndex."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        last_idx = avg.getLastIndex(50.0)
        self.assertEqual(last_idx, 2)

    def test_average_mass_delta(self):
        """Test getAverageMassDelta."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        delta = avg.getAverageMassDelta(50.0)
        self.assertAlmostEqual(delta, 0.025, delta=0.1)

    def test_most_abundant_mass_delta(self):
        """Test getMostAbundantMassDelta."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        delta = avg.getMostAbundantMassDelta(1000.0)
        self.assertAlmostEqual(delta, 0.0, delta=0.1)

    def test_snr_multiplication_factor(self):
        """Test getSNRMultiplicationFactor."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)

        factor = avg.getSNRMultiplicationFactor(75.0)
        self.assertGreater(factor, 0)

    def test_copy_constructor(self):
        """Test copy constructor."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg1 = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)
        avg1.setMaxIsotopeIndex(4)

        avg2 = pyopenms.PrecalAveragine(avg1)

        self.assertEqual(avg2.getApexIndex(75.0), avg1.getApexIndex(75.0))
        self.assertEqual(avg2.getMaxIsotopeIndex(), avg1.getMaxIsotopeIndex())

    def test_rna_averagine(self):
        """Test RNA averagine mode."""
        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg_peptide = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, False)
        avg_rna = pyopenms.PrecalAveragine(50.0, 100.0, 25.0, generator, True)

        iso_peptide = avg_peptide.get(75.0)
        iso_rna = avg_rna.get(75.0)

        self.assertGreater(iso_peptide.size(), 0)
        self.assertGreater(iso_rna.size(), 0)


class TestMassFeature(unittest.TestCase):
    """Test MassFeature (MassFeature_FDHS) struct."""

    def test_default_constructor(self):
        """Test default constructor."""
        mf = pyopenms.MassFeature_FDHS()
        self.assertIsNotNone(mf)

    def test_member_variables(self):
        """Test that all member variables can be accessed and modified."""
        mf = pyopenms.MassFeature_FDHS()

        # Set all member variables
        mf.index = 42
        mf.iso_offset = 1
        mf.scan_number = 100
        mf.min_scan_number = 90
        mf.max_scan_number = 110
        mf.rep_charge = 5
        mf.avg_mass = 15000.0
        mf.min_charge = 3
        mf.max_charge = 10
        mf.charge_count = 8
        mf.isotope_score = 0.95
        mf.qscore = 0.88
        mf.rep_mz = 1500.5
        mf.is_decoy = False
        mf.ms_level = 1

        # Verify all values
        self.assertEqual(mf.index, 42)
        self.assertEqual(mf.iso_offset, 1)
        self.assertEqual(mf.scan_number, 100)
        self.assertEqual(mf.min_scan_number, 90)
        self.assertEqual(mf.max_scan_number, 110)
        self.assertEqual(mf.rep_charge, 5)
        self.assertAlmostEqual(mf.avg_mass, 15000.0, places=2)
        self.assertEqual(mf.min_charge, 3)
        self.assertEqual(mf.max_charge, 10)
        self.assertEqual(mf.charge_count, 8)
        self.assertAlmostEqual(mf.isotope_score, 0.95, places=4)
        self.assertAlmostEqual(mf.qscore, 0.88, places=4)
        self.assertAlmostEqual(mf.rep_mz, 1500.5, places=2)
        self.assertFalse(mf.is_decoy)
        self.assertEqual(mf.ms_level, 1)

    def test_intensity_vectors(self):
        """Test per_charge_intensity and per_isotope_intensity vectors."""
        mf = pyopenms.MassFeature_FDHS()

        mf.per_charge_intensity = [100.0, 200.0, 300.0]
        mf.per_isotope_intensity = [50.0, 100.0, 75.0, 25.0]

        self.assertEqual(len(mf.per_charge_intensity), 3)
        self.assertEqual(len(mf.per_isotope_intensity), 4)
        self.assertAlmostEqual(mf.per_charge_intensity[1], 200.0, places=1)
        self.assertAlmostEqual(mf.per_isotope_intensity[2], 75.0, places=1)

    def test_comparison_operators(self):
        """Test comparison operators."""
        mf1 = pyopenms.MassFeature_FDHS()
        mf2 = pyopenms.MassFeature_FDHS()

        mf1.avg_mass = 1000.0
        mf2.avg_mass = 2000.0

        self.assertTrue(mf1 < mf2)
        self.assertFalse(mf2 < mf1)
        self.assertTrue(mf2 > mf1)
        self.assertFalse(mf1 > mf2)

    def test_equality_operator(self):
        """Test equality operator."""
        mf1 = pyopenms.MassFeature_FDHS()
        mf2 = pyopenms.MassFeature_FDHS()

        mf1.avg_mass = 1000.0
        mf2.avg_mass = 1000.0

        self.assertTrue(mf1 == mf2)

        mf2.avg_mass = 2000.0
        self.assertFalse(mf1 == mf2)


class TestIsobaricQuantities(unittest.TestCase):
    """Test IsobaricQuantities struct."""

    def test_default_constructor(self):
        """Test default constructor."""
        iq = pyopenms.IsobaricQuantities()
        self.assertIsNotNone(iq)

    def test_empty_method(self):
        """Test empty() method."""
        iq = pyopenms.IsobaricQuantities()

        # Should be empty initially
        self.assertTrue(iq.empty())

        # Add quantities
        iq.quantities = [100.0]
        self.assertFalse(iq.empty())

        # Clear quantities
        iq.quantities = []
        self.assertTrue(iq.empty())

    def test_member_variables(self):
        """Test all member variables."""
        iq = pyopenms.IsobaricQuantities()

        iq.scan = 500
        iq.rt = 120.5
        iq.precursor_mz = 750.25
        iq.precursor_mass = 1498.48
        iq.quantities = [100.0, 200.0, 150.0, 175.0]
        iq.merged_quantities = [450.0, 375.0]

        self.assertEqual(iq.scan, 500)
        self.assertAlmostEqual(iq.rt, 120.5, places=2)
        self.assertAlmostEqual(iq.precursor_mz, 750.25, places=2)
        self.assertAlmostEqual(iq.precursor_mass, 1498.48, places=2)
        self.assertEqual(len(iq.quantities), 4)
        self.assertEqual(len(iq.merged_quantities), 2)
        self.assertAlmostEqual(iq.quantities[0], 100.0, places=2)
        self.assertAlmostEqual(iq.quantities[3], 175.0, places=2)
        self.assertAlmostEqual(iq.merged_quantities[0], 450.0, places=2)


class TestPeakGroup(unittest.TestCase):
    """Test PeakGroup class."""

    def test_default_constructor(self):
        """Test default constructor."""
        pg = pyopenms.PeakGroup()
        self.assertIsNotNone(pg)

    def test_constructor_with_charge_range(self):
        """Test constructor with charge range parameters."""
        pg = pyopenms.PeakGroup(1, 10, True)
        self.assertIsNotNone(pg)
        self.assertTrue(pg.isPositive())

    def test_copy_constructor(self):
        """Test copy constructor."""
        pg1 = pyopenms.PeakGroup(1, 10, True)
        pg1.setScanNumber(5)
        pg1.setQscore(0.9)

        pg2 = pyopenms.PeakGroup(pg1)

        self.assertEqual(pg2.getScanNumber(), pg1.getScanNumber())
        self.assertAlmostEqual(pg2.getQscore(), pg1.getQscore(), places=4)

    def test_scan_number(self):
        """Test getScanNumber and setScanNumber."""
        pg = pyopenms.PeakGroup()
        pg.setScanNumber(42)
        self.assertEqual(pg.getScanNumber(), 42)

    def test_monoisotopic_mass(self):
        """Test getMonoMass and setMonoisotopicMass."""
        pg = pyopenms.PeakGroup()
        pg.setMonoisotopicMass(5000.0)
        self.assertAlmostEqual(pg.getMonoMass(), 5000.0, places=2)

    def test_intensity(self):
        """Test getIntensity."""
        pg = pyopenms.PeakGroup()
        # Intensity is computed from peaks, so test that it returns a float
        intensity = pg.getIntensity()
        self.assertIsInstance(intensity, float)

    def test_qscore(self):
        """Test getQscore and setQscore."""
        pg = pyopenms.PeakGroup()
        pg.setQscore(0.85)
        self.assertAlmostEqual(pg.getQscore(), 0.85, places=4)

    def test_qscore_2d(self):
        """Test getQscore2D and setQscore2D."""
        pg = pyopenms.PeakGroup()
        pg.setQscore2D(0.75)
        self.assertAlmostEqual(pg.getQscore2D(), 0.75, places=4)

    def test_isotope_cosine(self):
        """Test getIsotopeCosine and setIsotopeCosine."""
        pg = pyopenms.PeakGroup()
        pg.setIsotopeCosine(0.95)
        self.assertAlmostEqual(pg.getIsotopeCosine(), 0.95, places=4)

    def test_charge_score(self):
        """Test getChargeScore and setChargeScore."""
        pg = pyopenms.PeakGroup()
        pg.setChargeScore(0.8)
        self.assertAlmostEqual(pg.getChargeScore(), 0.8, places=4)

    def test_snr(self):
        """Test getSNR and setSNR."""
        pg = pyopenms.PeakGroup()
        pg.setSNR(50.0)
        self.assertAlmostEqual(pg.getSNR(), 50.0, places=2)

    def test_rep_abs_charge(self):
        """Test getRepAbsCharge and setRepAbsCharge."""
        pg = pyopenms.PeakGroup()
        pg.setRepAbsCharge(5)
        self.assertEqual(pg.getRepAbsCharge(), 5)

    def test_index(self):
        """Test getIndex and setIndex."""
        pg = pyopenms.PeakGroup()
        pg.setIndex(10)
        self.assertEqual(pg.getIndex(), 10)

    def test_feature_index(self):
        """Test getFeatureIndex and setFeatureIndex."""
        pg = pyopenms.PeakGroup()
        pg.setFeatureIndex(20)
        self.assertEqual(pg.getFeatureIndex(), 20)

    def test_qvalue(self):
        """Test getQvalue and setQvalue."""
        pg = pyopenms.PeakGroup()
        pg.setQvalue(0.05)
        self.assertAlmostEqual(pg.getQvalue(), 0.05, places=4)

    def test_target_decoy_type(self):
        """Test getTargetDecoyType and setTargetDecoyType."""
        pg = pyopenms.PeakGroup()
        pg.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.target)
        self.assertEqual(pg.getTargetDecoyType(), pyopenms.PeakGroup.TargetDecoyType.target)

        pg.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.noise_decoy)
        self.assertEqual(pg.getTargetDecoyType(), pyopenms.PeakGroup.TargetDecoyType.noise_decoy)

    def test_is_positive(self):
        """Test isPositive."""
        pg_pos = pyopenms.PeakGroup(1, 10, True)
        pg_neg = pyopenms.PeakGroup(1, 10, False)

        self.assertTrue(pg_pos.isPositive())
        self.assertFalse(pg_neg.isPositive())

    def test_is_targeted(self):
        """Test isTargeted and setTargeted."""
        pg = pyopenms.PeakGroup()
        self.assertFalse(pg.isTargeted())

        pg.setTargeted()
        self.assertTrue(pg.isTargeted())

    def test_avg_ppm_error(self):
        """Test getAvgPPMError and setAvgPPMError."""
        pg = pyopenms.PeakGroup()
        pg.setAvgPPMError(5.5)
        self.assertAlmostEqual(pg.getAvgPPMError(), 5.5, places=2)

    def test_avg_da_error(self):
        """Test getAvgDaError."""
        pg = pyopenms.PeakGroup()
        error = pg.getAvgDaError()
        self.assertIsInstance(error, float)

    def test_isotope_da_distance(self):
        """Test getIsotopeDaDistance and setIsotopeDaDistance."""
        pg = pyopenms.PeakGroup()
        pg.setIsotopeDaDistance(1.003)
        self.assertAlmostEqual(pg.getIsotopeDaDistance(), 1.003, places=3)

    def test_charge_intensity(self):
        """Test getChargeIntensity."""
        pg = pyopenms.PeakGroup(1, 10, True)
        intensity = pg.getChargeIntensity(5)
        self.assertIsInstance(intensity, float)

    def test_charge_snr(self):
        """Test getChargeSNR and setChargeSNR."""
        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setChargeSNR(5, 10.0)
        self.assertAlmostEqual(pg.getChargeSNR(5), 10.0, places=2)

    def test_charge_isotope_cosine(self):
        """Test getChargeIsotopeCosine and setChargeIsotopeCosine."""
        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setChargeIsotopeCosine(5, 0.9)
        self.assertAlmostEqual(pg.getChargeIsotopeCosine(5), 0.9, places=2)

    def test_isotope_intensities(self):
        """Test getIsotopeIntensities."""
        pg = pyopenms.PeakGroup()
        intensities = pg.getIsotopeIntensities()
        self.assertIsInstance(intensities, list)

    def test_mass_errors(self):
        """Test getMassErrors."""
        pg = pyopenms.PeakGroup()
        errors_ppm = pg.getMassErrors(True)
        errors_da = pg.getMassErrors(False)
        self.assertIsInstance(errors_ppm, list)
        self.assertIsInstance(errors_da, list)

    def test_container_operations(self):
        """Test container operations (size, empty, push_back, reserve)."""
        pg = pyopenms.PeakGroup(1, 10, True)

        self.assertTrue(pg.empty())
        self.assertEqual(pg.size(), 0)

        # Add a LogMzPeak
        p = pyopenms.Peak1D()
        p.setMZ(1000.0)
        p.setIntensity(100.0)
        log_peak = pyopenms.LogMzPeak(p, True)
        log_peak.abs_charge = 5
        log_peak.isotopeIndex = 0

        pg.push_back(log_peak)

        self.assertFalse(pg.empty())
        self.assertEqual(pg.size(), 1)

    def test_element_access(self):
        """Test element access via operator[]."""
        pg = pyopenms.PeakGroup(1, 10, True)

        p = pyopenms.Peak1D()
        p.setMZ(1000.0)
        p.setIntensity(100.0)
        log_peak = pyopenms.LogMzPeak(p, True)
        log_peak.abs_charge = 5

        pg.push_back(log_peak)

        retrieved = pg[0]
        self.assertAlmostEqual(retrieved.mz, 1000.0, places=2)

    def test_iteration(self):
        """Test iteration over peaks."""
        pg = pyopenms.PeakGroup(1, 10, True)

        for mz in [1000.0, 1001.0, 1002.0]:
            p = pyopenms.Peak1D()
            p.setMZ(mz)
            p.setIntensity(100.0)
            log_peak = pyopenms.LogMzPeak(p, True)
            log_peak.abs_charge = 5
            pg.push_back(log_peak)

        mz_values = [peak.mz for peak in pg]
        self.assertEqual(len(mz_values), 3)
        self.assertAlmostEqual(mz_values[0], 1000.0, places=2)
        self.assertAlmostEqual(mz_values[1], 1001.0, places=2)
        self.assertAlmostEqual(mz_values[2], 1002.0, places=2)

    def test_comparison_operators(self):
        """Test comparison operators."""
        pg1 = pyopenms.PeakGroup()
        pg2 = pyopenms.PeakGroup()

        pg1.setMonoisotopicMass(1000.0)
        pg2.setMonoisotopicMass(2000.0)

        self.assertTrue(pg1 < pg2)
        self.assertTrue(pg2 > pg1)

    def test_equality_operator(self):
        """Test equality operator."""
        pg1 = pyopenms.PeakGroup()
        pg2 = pyopenms.PeakGroup(pg1)

        self.assertTrue(pg1 == pg2)

    def test_sort(self):
        """Test sort method."""
        pg = pyopenms.PeakGroup(1, 10, True)

        # Add peaks in reverse order
        for mz in [1002.0, 1000.0, 1001.0]:
            p = pyopenms.Peak1D()
            p.setMZ(mz)
            p.setIntensity(100.0)
            log_peak = pyopenms.LogMzPeak(p, True)
            log_peak.abs_charge = 5
            pg.push_back(log_peak)

        pg.sort()

        # Verify sorted order
        mz_values = [peak.mz for peak in pg]
        self.assertLess(mz_values[0], mz_values[1])
        self.assertLess(mz_values[1], mz_values[2])


class TestDeconvolvedSpectrum(unittest.TestCase):
    """Test DeconvolvedSpectrum class."""

    def test_default_constructor(self):
        """Test default constructor."""
        ds = pyopenms.DeconvolvedSpectrum()
        self.assertIsNotNone(ds)

    def test_constructor_with_scan_number(self):
        """Test constructor with scan number."""
        ds = pyopenms.DeconvolvedSpectrum(42)
        self.assertEqual(ds.getScanNumber(), 42)

    def test_copy_constructor(self):
        """Test copy constructor."""
        ds1 = pyopenms.DeconvolvedSpectrum(10)
        ds2 = pyopenms.DeconvolvedSpectrum(ds1)
        self.assertEqual(ds2.getScanNumber(), ds1.getScanNumber())

    def test_scan_number(self):
        """Test getScanNumber."""
        ds = pyopenms.DeconvolvedSpectrum(100)
        self.assertEqual(ds.getScanNumber(), 100)

    def test_original_spectrum(self):
        """Test getOriginalSpectrum and setOriginalSpectrum."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        spec = pyopenms.MSSpectrum()
        p = pyopenms.Peak1D()
        p.setMZ(500.0)
        p.setIntensity(1000.0)
        spec.push_back(p)

        ds.setOriginalSpectrum(spec)

        retrieved = ds.getOriginalSpectrum()
        self.assertEqual(retrieved.size(), 1)

    def test_precursor(self):
        """Test getPrecursor and setPrecursor."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        prec = pyopenms.Precursor()
        prec.setMZ(750.0)
        prec.setCharge(5)

        ds.setPrecursor(prec)

        retrieved = ds.getPrecursor()
        self.assertAlmostEqual(retrieved.getMZ(), 750.0, places=2)
        self.assertEqual(retrieved.getCharge(), 5)

    def test_precursor_charge(self):
        """Test getPrecursorCharge."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        prec = pyopenms.Precursor()
        prec.setCharge(7)
        ds.setPrecursor(prec)

        self.assertEqual(ds.getPrecursorCharge(), 7)

    def test_precursor_scan_number(self):
        """Test getPrecursorScanNumber and setPrecursorScanNumber."""
        ds = pyopenms.DeconvolvedSpectrum(1)
        ds.setPrecursorScanNumber(50)
        self.assertEqual(ds.getPrecursorScanNumber(), 50)

    def test_precursor_peak_group(self):
        """Test getPrecursorPeakGroup and setPrecursorPeakGroup."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setMonoisotopicMass(5000.0)

        ds.setPrecursorPeakGroup(pg)

        retrieved = ds.getPrecursorPeakGroup()
        self.assertAlmostEqual(retrieved.getMonoMass(), 5000.0, places=2)

    def test_current_mass_limits(self):
        """Test getCurrentMaxMass and getCurrentMinMass."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        max_mass = ds.getCurrentMaxMass(10000.0)
        min_mass = ds.getCurrentMinMass(500.0)

        self.assertIsInstance(max_mass, float)
        self.assertIsInstance(min_mass, float)

    def test_current_max_abs_charge(self):
        """Test getCurrentMaxAbsCharge."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        max_charge = ds.getCurrentMaxAbsCharge(20)
        self.assertIsInstance(max_charge, int)

    def test_quantities(self):
        """Test getQuantities and setQuantities."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        iq = pyopenms.IsobaricQuantities()
        iq.scan = 100
        iq.quantities = [100.0, 200.0]

        ds.setQuantities(iq)

        retrieved = ds.getQuantities()
        self.assertEqual(retrieved.scan, 100)

    def test_is_decoy(self):
        """Test isDecoy."""
        ds = pyopenms.DeconvolvedSpectrum(1)
        is_decoy = ds.isDecoy()
        self.assertIsInstance(is_decoy, bool)

    def test_container_operations(self):
        """Test container operations (size, empty, clear, reserve, push_back, pop_back)."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        self.assertTrue(ds.empty())
        self.assertEqual(ds.size(), 0)

        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setMonoisotopicMass(5000.0)

        ds.push_back(pg)

        self.assertFalse(ds.empty())
        self.assertEqual(ds.size(), 1)

        ds.pop_back()

        self.assertTrue(ds.empty())
        self.assertEqual(ds.size(), 0)

    def test_set_peak_groups(self):
        """Test setPeakGroups."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        pgs = []
        for mass in [1000.0, 2000.0, 3000.0]:
            pg = pyopenms.PeakGroup(1, 10, True)
            pg.setMonoisotopicMass(mass)
            pgs.append(pg)

        ds.setPeakGroups(pgs)

        self.assertEqual(ds.size(), 3)

    def test_element_access(self):
        """Test element access via operator[]."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setMonoisotopicMass(5000.0)
        ds.push_back(pg)

        retrieved = ds[0]
        self.assertAlmostEqual(retrieved.getMonoMass(), 5000.0, places=2)

    def test_iteration(self):
        """Test iteration over peak groups."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        for mass in [1000.0, 2000.0, 3000.0]:
            pg = pyopenms.PeakGroup(1, 10, True)
            pg.setMonoisotopicMass(mass)
            ds.push_back(pg)

        masses = [pg.getMonoMass() for pg in ds]
        self.assertEqual(len(masses), 3)
        self.assertAlmostEqual(masses[0], 1000.0, places=2)
        self.assertAlmostEqual(masses[1], 2000.0, places=2)
        self.assertAlmostEqual(masses[2], 3000.0, places=2)

    def test_sort(self):
        """Test sort method."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        for mass in [3000.0, 1000.0, 2000.0]:
            pg = pyopenms.PeakGroup(1, 10, True)
            pg.setMonoisotopicMass(mass)
            ds.push_back(pg)

        ds.sort()

        masses = [pg.getMonoMass() for pg in ds]
        self.assertLess(masses[0], masses[1])
        self.assertLess(masses[1], masses[2])

    def test_sort_by_qscore(self):
        """Test sortByQscore method."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        for qscore in [0.5, 0.9, 0.7]:
            pg = pyopenms.PeakGroup(1, 10, True)
            pg.setQscore(qscore)
            ds.push_back(pg)

        ds.sortByQscore()

        qscores = [pg.getQscore() for pg in ds]
        # sortByQscore sorts in descending order
        self.assertGreaterEqual(qscores[0], qscores[1])
        self.assertGreaterEqual(qscores[1], qscores[2])

    def test_to_spectrum(self):
        """Test toSpectrum conversion."""
        ds = pyopenms.DeconvolvedSpectrum(1)

        pg = pyopenms.PeakGroup(1, 10, True)
        pg.setMonoisotopicMass(5000.0)
        ds.push_back(pg)

        # Original spectrum needed for conversion
        spec = pyopenms.MSSpectrum()
        spec.setRT(100.0)
        ds.setOriginalSpectrum(spec)

        result = ds.toSpectrum(1, 10.0, False)
        # MSSpectrum may be wrapped as _MSSpectrumDF, check for size method instead
        self.assertTrue(hasattr(result, 'size'))

    def test_comparison_operators(self):
        """Test comparison operators."""
        ds1 = pyopenms.DeconvolvedSpectrum(1)
        ds2 = pyopenms.DeconvolvedSpectrum(2)

        # Comparison is based on scan number
        self.assertTrue(ds1 < ds2)
        self.assertTrue(ds2 > ds1)

    def test_equality_operator(self):
        """Test equality operator."""
        ds1 = pyopenms.DeconvolvedSpectrum(1)
        ds2 = pyopenms.DeconvolvedSpectrum(ds1)

        self.assertTrue(ds1 == ds2)


class TestSpectralDeconvolution(unittest.TestCase):
    """Test SpectralDeconvolution class."""

    def test_default_constructor(self):
        """Test default constructor."""
        sd = pyopenms.SpectralDeconvolution()
        self.assertIsNotNone(sd)

    def test_copy_constructor(self):
        """Test copy constructor."""
        sd1 = pyopenms.SpectralDeconvolution()
        sd2 = pyopenms.SpectralDeconvolution(sd1)
        self.assertIsNotNone(sd2)

    def test_default_param_handler_interface(self):
        """Test DefaultParamHandler interface."""
        sd = pyopenms.SpectralDeconvolution()

        # Test that we can get and set parameters
        params = sd.getParameters()
        self.assertIsNotNone(params)

        # Test that default parameters exist
        defaults = sd.getDefaults()
        self.assertIsNotNone(defaults)

    def test_get_nominal_mass(self):
        """Test static getNominalMass method."""
        mass1 = 10000.0
        mass2 = 25000.0

        nominal1 = pyopenms.SpectralDeconvolution.getNominalMass(mass1)
        nominal2 = pyopenms.SpectralDeconvolution.getNominalMass(mass2)

        self.assertEqual(nominal1, 9995)
        self.assertEqual(nominal2, 24987)

    def test_calculate_averagine(self):
        """Test calculateAveragine method."""
        sd = pyopenms.SpectralDeconvolution()

        params = pyopenms.Param()
        params.setValue("max_mass", 2000.0)
        params.setValue("min_charge", 1)
        params.setValue("max_charge", 10)
        sd.setParameters(params)

        # Calculate peptide averagine
        sd.calculateAveragine(False)

        avg = sd.getAveragine()
        self.assertIsNotNone(avg)
        self.assertGreater(avg.getMaxIsotopeIndex(), 0)

    def test_calculate_rna_averagine(self):
        """Test calculateAveragine method with RNA mode."""
        sd = pyopenms.SpectralDeconvolution()

        params = pyopenms.Param()
        params.setValue("max_mass", 2000.0)
        params.setValue("min_charge", 1)
        params.setValue("max_charge", 10)
        sd.setParameters(params)

        # Calculate RNA averagine
        sd.calculateAveragine(True)

        avg = sd.getAveragine()
        self.assertIsNotNone(avg)
        self.assertGreater(avg.getMaxIsotopeIndex(), 0)

    def test_get_set_averagine(self):
        """Test getAveragine and setAveragine methods."""
        sd = pyopenms.SpectralDeconvolution()

        generator = pyopenms.CoarseIsotopePatternGenerator()
        avg = pyopenms.PrecalAveragine(50.0, 1000.0, 25.0, generator, False)
        avg.setMaxIsotopeIndex(10)

        sd.setAveragine(avg)

        retrieved = sd.getAveragine()
        self.assertEqual(retrieved.getMaxIsotopeIndex(), 10)

    def test_set_tolerance_estimation(self):
        """Test setToleranceEstimation method."""
        sd = pyopenms.SpectralDeconvolution()
        sd.setToleranceEstimation()
        # No direct way to verify, just check no exception

    def test_get_cosine(self):
        """Test static getCosine method."""
        a = [1.0, 2.0, 3.0, 4.0, 5.0]

        # Create isotope distribution for comparison
        generator = pyopenms.CoarseIsotopePatternGenerator()
        formula = pyopenms.EmpiricalFormula("C100H200N50O50")
        b = generator.run(formula)

        cosine = pyopenms.SpectralDeconvolution.getCosine(a, 0, 5, b, 0, 2)
        self.assertIsInstance(cosine, float)
        self.assertGreaterEqual(cosine, -1.0)
        self.assertLessEqual(cosine, 1.0)

    def test_set_target_decoy_type(self):
        """Test setTargetDecoyType method."""
        sd = pyopenms.SpectralDeconvolution()

        params = pyopenms.Param()
        params.setValue("max_mass", 2000.0)
        params.setValue("min_charge", 1)
        params.setValue("max_charge", 10)
        sd.setParameters(params)
        sd.calculateAveragine(False)

        ds = pyopenms.DeconvolvedSpectrum(1)

        sd.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.target, ds)


class TestFLASHDeconvAlgorithm(unittest.TestCase):
    """Test FLASHDeconvAlgorithm class."""

    def test_default_constructor(self):
        """Test default constructor."""
        algo = pyopenms.FLASHDeconvAlgorithm()
        self.assertIsNotNone(algo)

    def test_copy_constructor(self):
        """Test copy constructor."""
        algo1 = pyopenms.FLASHDeconvAlgorithm()
        algo2 = pyopenms.FLASHDeconvAlgorithm(algo1)
        self.assertEqual(algo2.getParameters(), algo1.getParameters())

    def test_default_param_handler_interface(self):
        """Test DefaultParamHandler interface."""
        algo = pyopenms.FLASHDeconvAlgorithm()

        params = algo.getParameters()
        self.assertIsNotNone(params)

        defaults = algo.getDefaults()
        self.assertIsNotNone(defaults)

    def test_progress_logger_interface(self):
        """Test ProgressLogger interface."""
        algo = pyopenms.FLASHDeconvAlgorithm()

        algo.setLogType(pyopenms.LogType.NONE)
        self.assertEqual(algo.getLogType(), pyopenms.LogType.NONE)

    def test_get_tolerances(self):
        """Test getTolerances method."""
        algo = pyopenms.FLASHDeconvAlgorithm()

        params = pyopenms.Param()
        params.setValue("SD:tol", [10.0, 5.0])
        algo.setParameters(params)

        tolerances = algo.getTolerances()
        self.assertEqual(len(tolerances), 2)
        self.assertAlmostEqual(tolerances[0], 10.0, places=2)
        self.assertAlmostEqual(tolerances[1], 5.0, places=2)

    def test_get_scan_number(self):
        """Test static getScanNumber method."""
        exp = pyopenms.MSExperiment()

        # Add a spectrum with native ID
        spec = pyopenms.MSSpectrum()
        spec.setNativeID("scan=42")
        exp.addSpectrum(spec)

        scan_num = pyopenms.FLASHDeconvAlgorithm.getScanNumber(exp, 0)
        self.assertIsInstance(scan_num, int)


class TestFLASHDeconvSpectrumFile(unittest.TestCase):
    """Test FLASHDeconvSpectrumFile class."""

    def test_default_constructor(self):
        """Test default constructor."""
        sf = pyopenms.FLASHDeconvSpectrumFile()
        self.assertIsNotNone(sf)

    def test_copy_constructor(self):
        """Test copy constructor."""
        sf1 = pyopenms.FLASHDeconvSpectrumFile()
        sf2 = pyopenms.FLASHDeconvSpectrumFile(sf1)
        self.assertIsNotNone(sf2)


class TestFLASHDeconvFeatureFile(unittest.TestCase):
    """Test FLASHDeconvFeatureFile class."""

    def test_default_constructor(self):
        """Test default constructor."""
        ff = pyopenms.FLASHDeconvFeatureFile()
        self.assertIsNotNone(ff)

    def test_copy_constructor(self):
        """Test copy constructor."""
        ff1 = pyopenms.FLASHDeconvFeatureFile()
        ff2 = pyopenms.FLASHDeconvFeatureFile(ff1)
        self.assertIsNotNone(ff2)


class TestTargetDecoyTypeEnum(unittest.TestCase):
    """Test TargetDecoyType enum."""

    def test_enum_values(self):
        """Test that all enum values exist."""
        # Enum members are distinct
        self.assertNotEqual(pyopenms.PeakGroup.TargetDecoyType.target, pyopenms.PeakGroup.TargetDecoyType.noise_decoy)
        self.assertNotEqual(pyopenms.PeakGroup.TargetDecoyType.noise_decoy, pyopenms.PeakGroup.TargetDecoyType.signal_decoy)

    def test_enum_assignment(self):
        """Test enum assignment to PeakGroup."""
        pg = pyopenms.PeakGroup()

        pg.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.target)
        self.assertEqual(pg.getTargetDecoyType(), pyopenms.PeakGroup.TargetDecoyType.target)

        pg.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.noise_decoy)
        self.assertEqual(pg.getTargetDecoyType(), pyopenms.PeakGroup.TargetDecoyType.noise_decoy)

        pg.setTargetDecoyType(pyopenms.PeakGroup.TargetDecoyType.signal_decoy)
        self.assertEqual(pg.getTargetDecoyType(), pyopenms.PeakGroup.TargetDecoyType.signal_decoy)


if __name__ == '__main__':
    unittest.main()
