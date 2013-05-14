import unittest
import os

import pyopenms

class TestSpectrumAccessOpenMS(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML")

    def test_readfile_content(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)
        saccess = pyopenms.SpectrumAccessOpenMS(exp)
        spectrum = saccess.getSpectrumById(0)
        mz = spectrum.getMZArray()
        intensity = spectrum.getIntensityArray()

        self.assertAlmostEqual(mz[0], 350.0000305)
        self.assertAlmostEqual(intensity[0], 0.0)
        self.assertAlmostEqual(mz[10], 358.075134277)
        self.assertAlmostEqual(intensity[10], 9210.931640625)

if __name__ == '__main__':
    unittest.main()
