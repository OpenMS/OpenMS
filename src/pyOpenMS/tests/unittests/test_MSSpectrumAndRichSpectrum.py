import unittest
import os
import pyopenms
from pyopenms import Precursor, PeptideIdentification, PeptideHit, AASequence, DriftTimeUnit
import numpy as np


class TestMSSpectrumDataDict(unittest.TestCase):
    """Tests for MSSpectrum.get_data_dict() method."""

    def test_basic_spectrum(self):
        """Test get_data_dict with a basic spectrum containing peaks."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 3000.0]))
        spec.setRT(123.45)
        spec.setMSLevel(1)
        spec.setNativeID("scan=1")

        data = spec.get_data_dict()

        # Check peaks
        self.assertEqual(len(data["mz"]), 3)
        np.testing.assert_array_almost_equal(data["mz"], [100.0, 200.0, 300.0])
        np.testing.assert_array_almost_equal(data["intensity"], [1000.0, 2000.0, 3000.0], decimal=0)

        # Check metadata
        np.testing.assert_array_almost_equal(data["rt"], [123.45, 123.45, 123.45])
        np.testing.assert_array_equal(data["ms_level"], [1, 1, 1])
        self.assertEqual(data["native_id"][0], "scan=1")

    def test_ion_mobility(self):
        """Test get_data_dict with ion mobility data."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))

        # Add ion mobility array
        fda = pyopenms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.push_back(0.5)
        spec.setFloatDataArrays([fda])
        spec.setDriftTimeUnit(DriftTimeUnit.VSSC)

        data = spec.get_data_dict()
        np.testing.assert_array_almost_equal(data["ion_mobility"], [0.5])
        # DriftTimeUnit.VSSC is converted to string "1/K0" (inverse reduced mobility)
        self.assertEqual(data["ion_mobility_unit"][0], "1/K0")

    def test_precursors(self):
        """Test get_data_dict with precursor information."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0], [10.0, 20.0]))
        spec.setMSLevel(2)

        p = Precursor()
        p.setMZ(300.0)
        p.setCharge(2)
        spec.setPrecursors([p])

        data = spec.get_data_dict()
        np.testing.assert_array_almost_equal(data["precursor_mz"], [300.0, 300.0])
        np.testing.assert_array_equal(data["precursor_charge"], [2, 2])

    def test_metadata_handling(self):
        """Test get_data_dict with various meta value types."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))
        spec.setMetaValue("int_val", 5)
        spec.setMetaValue("float_val", 3.14)
        spec.setMetaValue("str_val", "test")

        data = spec.get_data_dict()
        self.assertEqual(data["int_val"][0], 5)
        self.assertAlmostEqual(data["float_val"][0], 3.14, places=5)
        self.assertEqual(data["str_val"][0], "test")

    def test_empty_spectrum(self):
        """Test get_data_dict with an empty spectrum."""
        spec = pyopenms.MSSpectrum()
        data = spec.get_data_dict()
        self.assertEqual(len(data["mz"]), 0)
        self.assertEqual(len(data["intensity"]), 0)
        # Empty arrays should still have the keys
        self.assertIn("precursor_mz", data)
        self.assertIn("precursor_charge", data)

    def test_export_meta_values_flag(self):
        """Test that export_meta_values=False excludes meta values."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))
        spec.setMetaValue(b"test", 42)

        data = spec.get_data_dict(export_meta_values=False)
        self.assertNotIn("test", data)

        data_with_meta = spec.get_data_dict(export_meta_values=True)
        self.assertIn("test", data_with_meta)

    def test_no_precursor(self):
        """Test get_data_dict when no precursor is set."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))
        spec.setMSLevel(1)

        data = spec.get_data_dict()
        self.assertTrue(np.isnan(data["precursor_mz"][0]))
        self.assertEqual(data["precursor_charge"][0], 0)

    def test_no_ion_mobility(self):
        """Test get_data_dict when no ion mobility is present."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))

        data = spec.get_data_dict()
        self.assertTrue(np.isnan(data["ion_mobility"][0]))
        self.assertEqual(data["ion_mobility_unit"][0], "")

    def test_ion_annotations(self):
        """Test get_data_dict with ion annotations from StringDataArray."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 3000.0]))

        # Add ion annotations via StringDataArray
        sda = pyopenms.StringDataArray()
        sda.setName("IonNames")
        sda.push_back("y1")
        sda.push_back("b2")
        sda.push_back("y3-H2O")
        spec.setStringDataArrays([sda])

        data = spec.get_data_dict()
        self.assertIn("ion_annotation", data)
        self.assertEqual(len(data["ion_annotation"]), 3)
        self.assertEqual(data["ion_annotation"][0], "y1")
        self.assertEqual(data["ion_annotation"][1], "b2")
        self.assertEqual(data["ion_annotation"][2], "y3-H2O")

    def test_no_ion_annotations(self):
        """Test get_data_dict when no ion annotations are present."""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0], [1000.0, 2000.0]))

        data = spec.get_data_dict()
        self.assertIn("ion_annotation", data)
        self.assertEqual(len(data["ion_annotation"]), 2)
        self.assertEqual(data["ion_annotation"][0], "")
        self.assertEqual(data["ion_annotation"][1], "")


class TestMSSpectrumGetDF(unittest.TestCase):
    """Tests for MSSpectrum.get_df() method (pandas DataFrame export)."""

    def test_get_df_basic(self):
        """Test get_df returns a valid pandas DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")

        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0], [1000.0, 2000.0]))
        spec.setRT(123.45)
        spec.setMSLevel(1)

        df = spec.get_df()

        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 2)
        self.assertIn("mz", df.columns)
        self.assertIn("intensity", df.columns)
        self.assertIn("rt", df.columns)
        self.assertIn("ms_level", df.columns)

    def test_get_df_with_precursor(self):
        """Test get_df includes precursor information."""
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")

        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [1000.0]))
        spec.setMSLevel(2)

        p = Precursor()
        p.setMZ(500.0)
        p.setCharge(3)
        spec.setPrecursors([p])

        df = spec.get_df()

        self.assertEqual(df["precursor_mz"].iloc[0], 500.0)
        self.assertEqual(df["precursor_charge"].iloc[0], 3)

    def test_get_df_with_meta_values(self):
        """Test get_df includes meta values."""
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")

        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [1000.0]))
        spec.setMetaValue("score", 0.95)

        df = spec.get_df(export_meta_values=True)
        self.assertIn("score", df.columns)
        self.assertAlmostEqual(df["score"].iloc[0], 0.95, places=5)

    def test_get_df_without_meta_values(self):
        """Test get_df excludes meta values when requested."""
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")

        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [1000.0]))
        spec.setMetaValue("score", 0.95)

        df = spec.get_df(export_meta_values=False)
        self.assertNotIn("score", df.columns)

    def test_get_df_empty_spectrum(self):
        """Test get_df on empty spectrum returns empty DataFrame."""
        try:
            import pandas as pd
        except ImportError:
            self.skipTest("pandas not installed")

        spec = pyopenms.MSSpectrum()
        df = spec.get_df()

        self.assertIsInstance(df, pd.DataFrame)
        self.assertEqual(len(df), 0)


class TestMSSpectrumBasic(unittest.TestCase):
    """Basic tests for MSSpectrum functionality."""

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))

    def testMSSpectrum(self):
        """Test basic MSSpectrum operations."""
        spec = pyopenms.MSSpectrum()
        p = pyopenms.Peak1D()
        p.setMZ(500.0)
        p.setIntensity(1e5)
        spec.push_back(p)

        p_back, = list(spec)
        self.assertIsInstance(p_back, pyopenms.Peak1D)
        self.assertEqual(p_back.getMZ(), 500.0)
        self.assertEqual(p_back.getIntensity(), 1e5)

        spec.updateRanges()
        self.assertIsInstance(spec.getMinMZ(), float)
        self.assertIsInstance(spec.getMaxMZ(), float)
        self.assertIsInstance(spec.getMinIntensity(), float)
        self.assertIsInstance(spec.getMaxIntensity(), float)

        self.assertEqual(spec.getMinIntensity(), 1e5)
        self.assertEqual(spec.getMaxIntensity(), 1e5)

    def test_set_get_peaks(self):
        """Test set_peaks and get_peaks roundtrip."""
        spec = pyopenms.MSSpectrum()
        mzs = np.array([100.0, 200.0, 300.0], dtype=np.float64)
        intensities = np.array([1000.0, 2000.0, 3000.0], dtype=np.float32)

        spec.set_peaks((mzs, intensities))
        mzs_out, intensities_out = spec.get_peaks()

        np.testing.assert_array_almost_equal(mzs_out, mzs)
        np.testing.assert_array_almost_equal(intensities_out, intensities)


if __name__ == '__main__':
    unittest.main()
