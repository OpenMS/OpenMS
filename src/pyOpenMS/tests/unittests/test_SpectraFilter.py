import unittest
import os

import pyopenms

class TestSpectraFilter(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, self.exp)

    def test_map_NLargest(self):
        thisfilter = pyopenms.NLargest();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_NLargest(self):
        thisfilter = pyopenms.NLargest();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_Normalizer(self):
        thisfilter = pyopenms.Normalizer();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_Normalizer(self):
        thisfilter = pyopenms.Normalizer();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_Scaler(self):
        thisfilter = pyopenms.RankScaler();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_Scaler(self):
        thisfilter = pyopenms.RankScaler();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_SqrtScaler(self):
        thisfilter = pyopenms.SqrtScaler();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_SqrtScaler(self):
        thisfilter = pyopenms.SqrtScaler();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_ThresholdMower(self):
        thisfilter = pyopenms.ThresholdMower();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_ThresholdMower(self):
        thisfilter = pyopenms.ThresholdMower();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_WindowMower(self):
        thisfilter = pyopenms.WindowMower();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_WindowMower(self):
        thisfilter = pyopenms.WindowMower();

        new_firstspec = self.exp[0]
        thisfilter.filterPeakSpectrumForTopNInSlidingWindow(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

if __name__ == '__main__':
    unittest.main()
