import unittest
import os

import numpy as np
import pyopenms

class TestMSSpectrumAndRichSpectrum(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))


    def testMSSpectrum(self):
        spec = pyopenms.MSSpectrum()
        p = pyopenms.Peak1D()
        p.setMZ(500.0)
        p.setIntensity(1e5)
        spec.push_back(p)

        p_back, = list(spec)
        assert isinstance(p_back, pyopenms.Peak1D)
        assert p_back.getMZ() == 500.0
        assert p_back.getIntensity() == 1e5

        spec.updateRanges()
        assert isinstance(spec.getMinMZ(), float)
        assert isinstance(spec.getMaxMZ(), float)
        assert isinstance(spec.getMinIntensity(), float)
        assert isinstance(spec.getMaxIntensity(), float)

        assert spec.getMinIntensity() == 1e5
        assert spec.getMaxIntensity() == 1e5

    def testMSSpectrumGetPeaks(self):
        """Test optimized get_peaks method"""
        spec = pyopenms.MSSpectrum()
        mz_exp = [100.0, 200.0, 300.0]
        int_exp = [1000.0, 2000.0, 500.0]
        spec.set_peaks((mz_exp, int_exp))

        mz, intensities = spec.get_peaks()

        self.assertEqual(len(mz), 3)
        self.assertEqual(len(intensities), 3)
        self.assertEqual(mz.dtype, np.float64)
        self.assertEqual(intensities.dtype, np.float32)

        for m, e in zip(mz, mz_exp):
            self.assertAlmostEqual(m, e)
        for i, e in zip(intensities, int_exp):
            self.assertAlmostEqual(i, e, places=1)

    def testMSSpectrumGetPeaksEmpty(self):
        """Test get_peaks on empty spectrum"""
        spec = pyopenms.MSSpectrum()
        mz, intensities = spec.get_peaks()

        self.assertEqual(len(mz), 0)
        self.assertEqual(len(intensities), 0)
        self.assertEqual(mz.dtype, np.float64)
        self.assertEqual(intensities.dtype, np.float32)

    def testMSSpectrumGetMzArray(self):
        """Test get_mz_array method"""
        spec = pyopenms.MSSpectrum()
        mz_exp = [100.0, 200.0, 300.0]
        int_exp = [1000.0, 2000.0, 500.0]
        spec.set_peaks((mz_exp, int_exp))

        mz = spec.get_mz_array()
        self.assertEqual(len(mz), 3)
        self.assertEqual(mz.dtype, np.float64)
        for m, e in zip(mz, mz_exp):
            self.assertAlmostEqual(m, e)

    def testMSSpectrumGetIntensityArray(self):
        """Test get_intensity_array method"""
        spec = pyopenms.MSSpectrum()
        mz_exp = [100.0, 200.0, 300.0]
        int_exp = [1000.0, 2000.0, 500.0]
        spec.set_peaks((mz_exp, int_exp))

        intensities = spec.get_intensity_array()
        self.assertEqual(len(intensities), 3)
        self.assertEqual(intensities.dtype, np.float32)
        for i, e in zip(intensities, int_exp):
            self.assertAlmostEqual(i, e, places=1)

    def testMSSpectrumDriftTimeNoIM(self):
        """Test drift time methods when no IM data present"""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0], [1000.0, 2000.0]))

        # Should return None when no IM data
        self.assertFalse(spec.containsIMData())
        self.assertIsNone(spec.get_drift_time_array())
        self.assertIsNone(spec.get_drift_time_array_mv())
        self.assertIsNone(spec.get_drift_time_unit())

    def testMSSpectrumDriftTimeWithIM(self):
        """Test drift time methods with ion mobility data"""
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0, 200.0, 300.0], [1000.0, 2000.0, 500.0]))

        # Add ion mobility data via FloatDataArray
        fda = pyopenms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.push_back(1.5)
        fda.push_back(2.5)
        fda.push_back(3.5)
        spec.setFloatDataArrays([fda])

        # Should now have IM data
        self.assertTrue(spec.containsIMData())

        # Test get_drift_time_array (copy)
        drift = spec.get_drift_time_array()
        self.assertIsNotNone(drift)
        self.assertEqual(len(drift), 3)
        self.assertEqual(drift.dtype, np.float32)
        self.assertAlmostEqual(drift[0], 1.5, places=1)
        self.assertAlmostEqual(drift[1], 2.5, places=1)
        self.assertAlmostEqual(drift[2], 3.5, places=1)

        # Test get_drift_time_unit
        unit = spec.get_drift_time_unit()
        self.assertIsNotNone(unit)

    def testFloatDataArrayGetData(self):
        """Test FloatDataArray get_data method (returns copy - safe)"""
        fda = pyopenms.FloatDataArray()
        fda.push_back(1.0)
        fda.push_back(2.0)
        fda.push_back(3.0)

        # Get a copy (safe default)
        data_copy = fda.get_data()
        self.assertEqual(len(data_copy), 3)
        self.assertEqual(data_copy.dtype, np.float32)
        self.assertAlmostEqual(data_copy[0], 1.0, places=1)
        self.assertAlmostEqual(data_copy[1], 2.0, places=1)
        self.assertAlmostEqual(data_copy[2], 3.0, places=1)

        # Modify copy - original should be unchanged
        data_copy[0] = 100.0
        self.assertAlmostEqual(fda[0], 1.0, places=1)

    def testFloatDataArrayGetDataMv(self):
        """Test FloatDataArray get_data_mv method (memory view - fast, unsafe)"""
        fda = pyopenms.FloatDataArray()
        fda.push_back(1.0)
        fda.push_back(2.0)
        fda.push_back(3.0)

        # Get a view (fast but unsafe)
        data_view = fda.get_data_mv()
        self.assertEqual(len(data_view), 3)
        self.assertAlmostEqual(data_view[0], 1.0, places=1)
        self.assertAlmostEqual(data_view[1], 2.0, places=1)
        self.assertAlmostEqual(data_view[2], 3.0, places=1)

        # Modify view - original SHOULD change (it's a view)
        data_view[0] = 100.0
        self.assertAlmostEqual(fda[0], 100.0, places=1)

    def testFloatDataArrayEmpty(self):
        """Test FloatDataArray methods on empty array"""
        fda = pyopenms.FloatDataArray()

        # get_data returns empty array
        data_copy = fda.get_data()
        self.assertEqual(len(data_copy), 0)

        # get_data_mv returns None for empty
        data_view = fda.get_data_mv()
        self.assertIsNone(data_view)

if __name__ == '__main__':
    unittest.main()
