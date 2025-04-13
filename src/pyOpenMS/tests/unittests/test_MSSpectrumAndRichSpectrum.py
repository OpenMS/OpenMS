import unittest
import os
import pyopenms
from pyopenms import Precursor, PeptideIdentification, PeptideHit, AASequence , DriftTimeUnit
import numpy as np

class TestMSSpectrumAndRichSpectrum(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
    
    def test_ion_mobility(self):
        spec = pyopenms.MSSpectrum()
        spec.set_peaks(([100.0], [10.0]))
        
        # Add ion mobility array
        fda = pyopenms.FloatDataArray()
        fda.setName("Ion Mobility")
        fda.push_back(0.5)
        spec.setFloatDataArrays([fda])
        spec.setDriftTimeUnit(DriftTimeUnit.VSSC)

        data = spec.get_data_dict()
        self.assertTrue(np.array_equal(data["ion_mobility"], [0.5]))
        self.assertEqual(data["ion_mobility_unit"][0], "VSSC")

    def test_precursors(self):
        spec = pyopenms.MSSpectrum()
        p = Precursor()
        p.setMZ(300.0)
        p.setCharge(2)
        spec.setPrecursors([p])
        
        data = spec.get_data_dict()
        self.assertEqual(data["precursor_mz"][0], 300.0)
        self.assertEqual(data["precursor_charge"][0], 2)

    def test_metadata_handling(self):
        spec = pyopenms.MSSpectrum()
        spec.setMetaValue("int_val", 5) 
        spec.setMetaValue("float_val", 3.14)
        spec.setMetaValue("str_val", "test")
        spec.setMetaValue("bool_val", True)
        
        data = spec.get_data_dict()
        self.assertEqual(data["int_val"][0], 5)
        self.assertEqual(data["float_val"][0], 3.14)
        self.assertEqual(data["str_val"][0], "test")
        self.assertEqual(data["bool_val"][0], True)

    def test_empty_spectrum(self):
        spec = pyopenms.MSSpectrum()
        data = spec.get_data_dict()
        self.assertEqual(len(data["mz"]), 0)
        self.assertEqual(data["precursor_mz"][0], np.nan)

    def test_export_meta_values_flag(self):
        spec = pyopenms.MSSpectrum()
        spec.setMetaValue(b"test", 42)
        data = spec.get_data_dict(export_meta_values=False)
        self.assertNotIn("test", data)


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

if __name__ == '__main__':
    unittest.main()
