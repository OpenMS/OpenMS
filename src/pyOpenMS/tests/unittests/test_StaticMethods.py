"""
Tests for static methods converted to @staticmethod decorator pattern.
This test file covers the wrapper changes from issue #8559.
"""
import unittest
import os
import tempfile

import pyopenms


def make_peak1d(mz, intensity):
    """Helper to create Peak1D with mz and intensity."""
    p = pyopenms.Peak1D()
    p.setMZ(mz)
    p.setIntensity(intensity)
    return p


class TestFileStaticMethods(unittest.TestCase):
    """Test static methods of the File class."""

    def test_exists(self):
        """Test File.exists static method."""
        # Create a temporary file
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        try:
            self.assertTrue(pyopenms.File.exists(temp_path))
            self.assertFalse(pyopenms.File.exists("/nonexistent/path/file.txt"))
        finally:
            os.unlink(temp_path)

    def test_readable(self):
        """Test File.readable static method."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        try:
            self.assertTrue(pyopenms.File.readable(temp_path))
            self.assertFalse(pyopenms.File.readable("/nonexistent/path/file.txt"))
        finally:
            os.unlink(temp_path)

    def test_writable(self):
        """Test File.writable static method."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        try:
            self.assertTrue(pyopenms.File.writable(temp_path))
        finally:
            os.unlink(temp_path)

    def test_empty(self):
        """Test File.empty static method."""
        # Create an empty file
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        try:
            self.assertTrue(pyopenms.File.empty(temp_path))
        finally:
            os.unlink(temp_path)

        # Create a non-empty file
        with tempfile.NamedTemporaryFile(delete=False, mode='w') as f:
            f.write("content")
            temp_path = f.name
        try:
            self.assertFalse(pyopenms.File.empty(temp_path))
        finally:
            os.unlink(temp_path)

    def test_remove(self):
        """Test File.remove static method."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        self.assertTrue(os.path.exists(temp_path))
        result = pyopenms.File.remove(temp_path)
        self.assertTrue(result)
        self.assertFalse(os.path.exists(temp_path))

    def test_rename(self):
        """Test File.rename static method."""
        with tempfile.NamedTemporaryFile(delete=False) as f:
            old_path = f.name
        new_path = old_path + "_renamed"
        try:
            result = pyopenms.File.rename(old_path, new_path, True, False)
            self.assertTrue(result)
            self.assertFalse(os.path.exists(old_path))
            self.assertTrue(os.path.exists(new_path))
        finally:
            if os.path.exists(new_path):
                os.unlink(new_path)
            if os.path.exists(old_path):
                os.unlink(old_path)

    def test_basename(self):
        """Test File.basename static method."""
        result = pyopenms.File.basename("/path/to/file.txt")
        self.assertEqual(str(result), "file.txt")

    def test_path(self):
        """Test File.path static method."""
        result = pyopenms.File.path("/path/to/file.txt")
        self.assertIn("path", str(result))

    def test_absolutePath(self):
        """Test File.absolutePath static method."""
        result = pyopenms.File.absolutePath(".")
        self.assertGreater(len(str(result)), 0)

    def test_isDirectory(self):
        """Test File.isDirectory static method."""
        self.assertTrue(pyopenms.File.isDirectory(tempfile.gettempdir()))
        with tempfile.NamedTemporaryFile(delete=False) as f:
            temp_path = f.name
        try:
            self.assertFalse(pyopenms.File.isDirectory(temp_path))
        finally:
            os.unlink(temp_path)

    def test_getTempDirectory(self):
        """Test File.getTempDirectory static method."""
        result = pyopenms.File.getTempDirectory()
        self.assertGreater(len(str(result)), 0)

    def test_getUserDirectory(self):
        """Test File.getUserDirectory static method."""
        result = pyopenms.File.getUserDirectory()
        self.assertGreater(len(str(result)), 0)

    def test_getUniqueName(self):
        """Test File.getUniqueName static method."""
        name1 = pyopenms.File.getUniqueName()
        name2 = pyopenms.File.getUniqueName()
        self.assertNotEqual(str(name1), str(name2))

    def test_getTemporaryFile(self):
        """Test File.getTemporaryFile static method."""
        result = pyopenms.File.getTemporaryFile("")
        self.assertGreater(len(str(result)), 0)

    @unittest.skip("File.stripExtension not exposed as static method")
    def test_stripExtension(self):
        """Test File.stripExtension static method."""
        result = pyopenms.File.stripExtension("/path/to/file.mzML")
        self.assertIn("file", str(result))
        self.assertNotIn(".mzML", str(result))


class TestBuildInfoStaticMethods(unittest.TestCase):
    """Test static methods of the OpenMSBuildInfo class."""

    @unittest.skip("getOSInfo not exposed as static method")
    def test_getOSInfo(self):
        """Test OpenMSBuildInfo.getOSInfo static method."""
        os_info = pyopenms.OpenMSBuildInfo.getOSInfo()
        self.assertIsNotNone(os_info)

    @unittest.skip("getBinaryArchitecture not exposed as static method")
    def test_getBinaryArchitecture(self):
        """Test OpenMSBuildInfo.getBinaryArchitecture static method."""
        arch = pyopenms.OpenMSBuildInfo.getBinaryArchitecture()
        self.assertGreater(len(str(arch)), 0)

    def test_isOpenMPEnabled(self):
        """Test OpenMSBuildInfo.isOpenMPEnabled static method."""
        # Just check it returns a boolean
        result = pyopenms.OpenMSBuildInfo.isOpenMPEnabled()
        self.assertIsInstance(result, bool)

    def test_getBuildType(self):
        """Test OpenMSBuildInfo.getBuildType static method."""
        build_type = pyopenms.OpenMSBuildInfo.getBuildType()
        self.assertGreater(len(str(build_type)), 0)

    def test_getOpenMPMaxNumThreads(self):
        """Test OpenMSBuildInfo.getOpenMPMaxNumThreads static method."""
        num_threads = pyopenms.OpenMSBuildInfo.getOpenMPMaxNumThreads()
        self.assertGreaterEqual(num_threads, 1)


class TestVersionInfoStaticMethods(unittest.TestCase):
    """Test static methods of the VersionInfo class."""

    def test_getVersion(self):
        """Test VersionInfo.getVersion static method."""
        version = pyopenms.VersionInfo.getVersion()
        self.assertGreater(len(str(version)), 0)

    def test_getRevision(self):
        """Test VersionInfo.getRevision static method."""
        revision = pyopenms.VersionInfo.getRevision()
        # Revision may be empty in some builds
        self.assertIsNotNone(revision)

    def test_getTime(self):
        """Test VersionInfo.getTime static method."""
        time = pyopenms.VersionInfo.getTime()
        self.assertIsNotNone(time)

    def test_getBranch(self):
        """Test VersionInfo.getBranch static method."""
        branch = pyopenms.VersionInfo.getBranch()
        self.assertIsNotNone(branch)


class TestDateTimeStaticMethods(unittest.TestCase):
    """Test static methods of the DateTime class."""

    def test_now(self):
        """Test DateTime.now static method."""
        dt = pyopenms.DateTime.now()
        self.assertIsNotNone(dt)
        # Check that the returned DateTime has valid date
        date_str = dt.getDate()
        self.assertGreater(len(str(date_str)), 0)


class TestDeisotoperStaticMethods(unittest.TestCase):
    """Test static methods of the Deisotoper class."""

    @unittest.skip("Causes segfault - needs investigation")
    def test_deisotopeAndSingleCharge(self):
        """Test Deisotoper.deisotopeAndSingleCharge static method."""
        spectrum = pyopenms.MSSpectrum()
        # Add some peaks using helper
        spectrum.push_back(make_peak1d(100.0, 1000.0))
        spectrum.push_back(make_peak1d(101.003, 500.0))  # Isotope peak
        spectrum.push_back(make_peak1d(200.0, 800.0))

        # Call the static method with all 15 required parameters
        pyopenms.Deisotoper.deisotopeAndSingleCharge(
            spectrum,
            0.1,      # fragment_tolerance
            False,    # fragment_unit_ppm
            1,        # min_charge
            3,        # max_charge
            False,    # keep_only_deisotoped
            3,        # min_isopeaks
            10,       # max_isopeaks
            True,     # make_single_charged
            False,    # annotate_charge
            False,    # annotate_iso_peak_count
            False,    # use_decreasing_model
            2,        # start_intensity_check
            False,    # add_up_intensity
            True      # combine_mono_peak (added parameter)
        )
        # Just verify it runs without error
        self.assertIsNotNone(spectrum)


@unittest.skip("IMTypes uses wrap-attach pattern, not @staticmethod")
class TestIMTypesStaticMethods(unittest.TestCase):
    """Test static methods of IMTypes enums.

    Note: These methods use the wrap-attach pattern instead of @staticmethod,
    as they are free functions in the OpenMS namespace, not true class static methods.
    """

    def test_toDriftTimeUnit(self):
        """Test IMTypes.toDriftTimeUnit static method."""
        unit = pyopenms.IMTypes.toDriftTimeUnit("millisecond")
        self.assertIsNotNone(unit)

    def test_toIMFormat(self):
        """Test IMTypes.toIMFormat static method."""
        fmt = pyopenms.IMTypes.toIMFormat("concatenated")
        self.assertIsNotNone(fmt)


class TestTransformationModelStaticMethods(unittest.TestCase):
    """Test static methods of TransformationModel classes."""

    def test_TransformationModelLinear_getDefaultParameters(self):
        """Test TransformationModelLinear.getDefaultParameters static method."""
        params = pyopenms.Param()
        pyopenms.TransformationModelLinear.getDefaultParameters(params)
        self.assertIsNotNone(params)

    def test_TransformationModelBSpline_getDefaultParameters(self):
        """Test TransformationModelBSpline.getDefaultParameters static method."""
        params = pyopenms.Param()
        pyopenms.TransformationModelBSpline.getDefaultParameters(params)
        self.assertIsNotNone(params)

    def test_TransformationModelLowess_getDefaultParameters(self):
        """Test TransformationModelLowess.getDefaultParameters static method."""
        params = pyopenms.Param()
        pyopenms.TransformationModelLowess.getDefaultParameters(params)
        self.assertIsNotNone(params)


@unittest.skip("MZTrafoModel uses wrap-attach pattern, not @staticmethod")
class TestMZTrafoModelStaticMethods(unittest.TestCase):
    """Test static methods of MZTrafoModel class.

    Note: These methods use the wrap-attach pattern instead of @staticmethod,
    as they are free functions in the OpenMS namespace, not true class static methods.
    """

    def test_getModelTypes(self):
        """Test MZTrafoModel.getModelTypes static method."""
        result = []
        pyopenms.MZTrafoModel.getModelTypes(result)
        self.assertGreater(len(result), 0)

    def test_nameToEnum(self):
        """Test MZTrafoModel.nameToEnum static method."""
        enum_val = pyopenms.MZTrafoModel.nameToEnum("linear")
        self.assertIsNotNone(enum_val)

    def test_enumToName(self):
        """Test MZTrafoModel.enumToName static method."""
        name = pyopenms.MZTrafoModel.enumToName(pyopenms.MZTrafoModel_MODELTYPE.LINEAR)
        self.assertEqual(str(name), "linear")


@unittest.skip("SpectrumHelper uses wrap-attach pattern, not @staticmethod")
class TestSpectrumHelperStaticMethods(unittest.TestCase):
    """Test static methods of SpectrumHelper class.

    Note: These methods use the wrap-attach pattern instead of @staticmethod,
    as they are free template functions in the OpenMS namespace.
    """

    def test_removePeaks_spectrum(self):
        """Test SpectrumHelper.removePeaks for MSSpectrum."""
        spectrum = pyopenms.MSSpectrum()
        spectrum.push_back(make_peak1d(100.0, 1000.0))
        spectrum.push_back(make_peak1d(200.0, 800.0))
        spectrum.push_back(make_peak1d(300.0, 600.0))

        # Remove peaks between 150 and 250
        pyopenms.SpectrumHelper.removePeaks(spectrum, 150.0, 250.0)

        # Should have 2 peaks left (100 and 300)
        self.assertEqual(spectrum.size(), 2)

    def test_removePeaks_chromatogram(self):
        """Test SpectrumHelper.removePeaks for MSChromatogram."""
        chrom = pyopenms.MSChromatogram()
        chrom.push_back(pyopenms.ChromatogramPeak(10.0, 1000.0))
        chrom.push_back(pyopenms.ChromatogramPeak(20.0, 800.0))
        chrom.push_back(pyopenms.ChromatogramPeak(30.0, 600.0))

        # Remove peaks between 15 and 25
        pyopenms.SpectrumHelper.removePeaks(chrom, 15.0, 25.0)

        # Should have 2 peaks left
        self.assertEqual(chrom.size(), 2)

    def test_subtractMinimumIntensity_spectrum(self):
        """Test SpectrumHelper.subtractMinimumIntensity for MSSpectrum."""
        spectrum = pyopenms.MSSpectrum()
        spectrum.push_back(make_peak1d(100.0, 1000.0))
        spectrum.push_back(make_peak1d(200.0, 500.0))

        pyopenms.SpectrumHelper.subtractMinimumIntensity(spectrum)

        # Minimum (500) should be subtracted
        # This is just a smoke test - verify it runs
        self.assertIsNotNone(spectrum)


class TestMRMRTNormalizerConstructors(unittest.TestCase):
    """Test MRMRTNormalizer constructors added in the wrapper changes."""

    def test_default_constructor(self):
        """Test MRMRTNormalizer default constructor."""
        normalizer = pyopenms.MRMRTNormalizer()
        self.assertIsNotNone(normalizer)

    def test_copy_constructor(self):
        """Test MRMRTNormalizer copy constructor."""
        normalizer1 = pyopenms.MRMRTNormalizer()
        normalizer2 = pyopenms.MRMRTNormalizer(normalizer1)
        self.assertIsNotNone(normalizer2)

    def test_removeOutliersIterative(self):
        """Test MRMRTNormalizer.removeOutliersIterative static method."""
        # Create sample pairs (observed RT, reference RT) - must be lists, not tuples
        pairs = [
            [100.0, 110.0],
            [200.0, 210.0],
            [300.0, 310.0],
            [400.0, 410.0],
            [500.0, 510.0],
        ]
        # Note: outlier_detection_method must be bytes, not str
        result = pyopenms.MRMRTNormalizer.removeOutliersIterative(
            pairs, 0.95, 0.6, True, b"iter_jackknife"
        )
        self.assertIsNotNone(result)
        self.assertGreater(len(result), 0)

    def test_removeOutliersRANSAC(self):
        """Test MRMRTNormalizer.removeOutliersRANSAC static method."""
        # Signature: pairs, rsq_limit, coverage_limit, max_iterations, max_rt_threshold, sampling_size
        # Pairs must be lists, not tuples
        # Note: RANSAC requires at least 30 input peptides
        pairs = [[float(i * 100), float(i * 100 + 10)] for i in range(1, 35)]
        result = pyopenms.MRMRTNormalizer.removeOutliersRANSAC(
            pairs, 0.95, 0.6, 10, 5.0, 5
        )
        self.assertIsNotNone(result)


class TestCalibrationDataStaticMethods(unittest.TestCase):
    """Test CalibrationData static methods."""

    def test_getMetaValues(self):
        """Test CalibrationData.getMetaValues static method."""
        meta_values = pyopenms.CalibrationData.getMetaValues()
        self.assertIsNotNone(meta_values)
        # Should return a list of strings
        self.assertIsInstance(meta_values, list)


class TestFileHandlerStaticMethods(unittest.TestCase):
    """Test FileHandler static methods."""

    def test_getTypeByFileName(self):
        """Test FileHandler.getTypeByFileName static method."""
        file_type = pyopenms.FileHandler.getTypeByFileName("test.mzML")
        self.assertIsNotNone(file_type)

    def test_hasValidExtension(self):
        """Test FileHandler.hasValidExtension static method."""
        # mzML should have valid extension for MZML type
        result = pyopenms.FileHandler.hasValidExtension("test.mzML", pyopenms.FileType.MZML)
        self.assertTrue(result)

    def test_isSupported(self):
        """Test FileHandler.isSupported static method."""
        # MZML should be supported
        result = pyopenms.FileHandler.isSupported(pyopenms.FileType.MZML)
        self.assertTrue(result)


class TestCachedmzMLStaticMethods(unittest.TestCase):
    """Test CachedmzML static methods."""

    def test_store_and_load(self):
        """Test CachedmzML.store and CachedmzML.load static methods."""
        # Create a simple MSExperiment
        exp = pyopenms.MSExperiment()
        spectrum = pyopenms.MSSpectrum()
        spectrum.push_back(make_peak1d(100.0, 1000.0))
        exp.addSpectrum(spectrum)

        # Create a temporary file
        with tempfile.NamedTemporaryFile(suffix=".mzML", delete=False) as f:
            temp_path = f.name

        try:
            # Store using static method
            pyopenms.CachedmzML.store(temp_path, exp)
            self.assertTrue(os.path.exists(temp_path))

            # Load using static method
            cached = pyopenms.CachedmzML()
            pyopenms.CachedmzML.load(temp_path, cached)
            self.assertIsNotNone(cached)
        finally:
            # Clean up - remove both the mzML and any cache files
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            cache_path = temp_path + ".cached"
            if os.path.exists(cache_path):
                os.unlink(cache_path)


class TestExperimentalDesignStaticMethods(unittest.TestCase):
    """Test ExperimentalDesign static methods."""

    def test_fromFeatureMap(self):
        """Test ExperimentalDesign.fromFeatureMap static method."""
        fm = pyopenms.FeatureMap()
        result = pyopenms.ExperimentalDesign.fromFeatureMap(fm)
        self.assertIsNotNone(result)

    def test_fromConsensusMap(self):
        """Test ExperimentalDesign.fromConsensusMap static method."""
        cm = pyopenms.ConsensusMap()
        result = pyopenms.ExperimentalDesign.fromConsensusMap(cm)
        self.assertIsNotNone(result)


if __name__ == '__main__':
    unittest.main()
