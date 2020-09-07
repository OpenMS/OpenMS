import unittest
import os

import pyopenms

class TestOpenSwathDataStructures(unittest.TestCase):

    def setUp(self):
        pass

    def test_spectrum(self):
        spectrum = pyopenms.Interfaces.Spectrum()
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

    def test_chromatogram(self):
        chromatogram = pyopenms.Interfaces.Chromatogram()
        rt_exp = [1,2,3]
        int_exp = [4,5,6]

        chromatogram.setTimeArray(rt_exp)
        chromatogram.setIntensityArray(int_exp)
        time = chromatogram.getTimeArray()
        intensity = chromatogram.getIntensityArray()

        for m,e in zip(time, rt_exp):
            self.assertAlmostEqual(m,e)

        for i,e in zip(intensity, int_exp):
            self.assertAlmostEqual(i,e)

if __name__ == '__main__':
    unittest.main()
