import unittest
import pyopenms


class TestMath(unittest.TestCase):

    def test_getPPM(self):
        """Test Math.getPPM function"""
        # Test basic PPM calculation
        # If observed m/z is 1000.001 and reference is 1000.0,
        # PPM = (1000.001 - 1000.0) / 1000.0 * 1e6 = 1.0
        mz_obs = 1000.001
        mz_ref = 1000.0
        ppm = pyopenms.Math.getPPM(mz_obs, mz_ref)
        self.assertAlmostEqual(ppm, 1.0, places=5)

        # Test negative PPM (observed < reference)
        mz_obs = 999.999
        mz_ref = 1000.0
        ppm = pyopenms.Math.getPPM(mz_obs, mz_ref)
        self.assertAlmostEqual(ppm, -1.0, places=5)

    def test_getPPMAbs(self):
        """Test Math.getPPMAbs function"""
        # Test absolute PPM calculation
        mz_obs = 1000.001
        mz_ref = 1000.0
        ppm_abs = pyopenms.Math.getPPMAbs(mz_obs, mz_ref)
        self.assertAlmostEqual(ppm_abs, 1.0, places=5)

        # Test that negative PPM is returned as positive
        mz_obs = 999.999
        mz_ref = 1000.0
        ppm_abs = pyopenms.Math.getPPMAbs(mz_obs, mz_ref)
        self.assertAlmostEqual(ppm_abs, 1.0, places=5)

    def test_ppmToMass(self):
        """Test Math.ppmToMass function"""
        # Test conversion from PPM to mass difference
        # If PPM = 5.0 and mz_ref = 1000.0,
        # mass_diff = 5.0 / 1e6 * 1000.0 = 0.005
        ppm = 5.0
        mz_ref = 1000.0
        mass_diff = pyopenms.Math.ppmToMass(ppm, mz_ref)
        self.assertAlmostEqual(mass_diff, 0.005, places=6)

        # Test negative PPM
        ppm = -5.0
        mz_ref = 1000.0
        mass_diff = pyopenms.Math.ppmToMass(ppm, mz_ref)
        self.assertAlmostEqual(mass_diff, -0.005, places=6)

    def test_ppmToMassAbs(self):
        """Test Math.ppmToMassAbs function"""
        # Test absolute mass difference
        ppm = 5.0
        mz_ref = 1000.0
        mass_diff_abs = pyopenms.Math.ppmToMassAbs(ppm, mz_ref)
        self.assertAlmostEqual(mass_diff_abs, 0.005, places=6)

        # Test that negative PPM returns positive mass difference
        ppm = -5.0
        mz_ref = 1000.0
        mass_diff_abs = pyopenms.Math.ppmToMassAbs(ppm, mz_ref)
        self.assertAlmostEqual(mass_diff_abs, 0.005, places=6)

    def test_roundtrip(self):
        """Test roundtrip conversion PPM <-> mass"""
        mz_obs = 1000.005
        mz_ref = 1000.0
        
        # Calculate PPM
        ppm = pyopenms.Math.getPPM(mz_obs, mz_ref)
        
        # Convert back to mass difference
        mass_diff = pyopenms.Math.ppmToMass(ppm, mz_ref)
        
        # Should equal the original difference
        self.assertAlmostEqual(mass_diff, mz_obs - mz_ref, places=6)


if __name__ == '__main__':
    unittest.main()
