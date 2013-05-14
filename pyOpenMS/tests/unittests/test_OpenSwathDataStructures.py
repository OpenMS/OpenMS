import unittest
import os

import pyopenms

class TestOpenSwathDataStructures(unittest.TestCase):

    def setUp(self):
        pass

    def test_spectrum(self):
        exp = pyopenms.MSExperiment()
        spectrum = pyopenms.SpectrumAccessOpenMS(exp)
        spectrum = pyopenms.Spectrum()
        mz_exp = [1,2,3]
        int_exp = [4,5,6]

        spectrum.setMZArray(mz_exp)
        spectrum.setIntensityArray(int_exp)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        for m,e in zip(mz, mz_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

if __name__ == '__main__':
    unittest.main()
