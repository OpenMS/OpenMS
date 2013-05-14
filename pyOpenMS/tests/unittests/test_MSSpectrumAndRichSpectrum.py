import unittest
import os

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

    def testRichMSSpectrum(self):
        spec = pyopenms.RichMSSpectrum()
        p = pyopenms.RichPeak1D()
        p.setMZ(500.0)
        p.setIntensity(1e5)
        spec.push_back(p)

        p_back, = list(spec)
        assert isinstance(p_back, pyopenms.RichPeak1D)
        assert p_back.getMZ() == 500.0
        assert p_back.getIntensity() == 1e5


if __name__ == '__main__':
    unittest.main()
