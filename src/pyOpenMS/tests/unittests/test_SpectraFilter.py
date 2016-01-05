import unittest
import os

import pyopenms

class TestSpectraFilter(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, self.exp)

    def test_map_BernNorm(self):
        thisfilter = pyopenms.BernNorm();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_BernNorm(self):
        thisfilter = pyopenms.BernNorm();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_MarkerMower(self):
        thisfilter = pyopenms.MarkerMower();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        # this deletes the spectrum ... 
        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

    def test_spectrum_MarkerMower(self):
        thisfilter = pyopenms.MarkerMower();

        new_firstspec = self.exp[0]
        self.assertNotEqual(new_firstspec.size(), 0)
        thisfilter.filterSpectrum(new_firstspec)

        # this deletes the spectrum ... 
        self.assertEqual(new_firstspec.size(), 0)

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

    def test_map_ParentPeakMower(self):
        thisfilter = pyopenms.ParentPeakMower();
        param = thisfilter.getDefaults()
        param.setValue(b"window_size", 20.0, b"")
        param.setValue(b"set_to_zero", True, b"")
        thisfilter.setParameters(param)


        self.assertEqual(self.exp.size(), 5)
        self.assertEqual(self.exp[0].getMSLevel(), 1)

        map = []
        for spec in self.exp:
            if spec.getMSLevel() != 1:
                map.append(spec)

        self.exp.setSpectra(map)

        self.assertEqual(self.exp.size(), 3)
        self.assertEqual(self.exp[0].getMSLevel(), 2)
        self.assertEqual(self.exp[1].getMSLevel(), 2)
        self.assertEqual(self.exp[2].getMSLevel(), 2)

        old_firstspec = self.exp[0]
        old_secondspec = self.exp[1]

        thisfilter.filterPeakMap(self.exp)

        a = old_firstspec
        b = self.exp[0]
        self.assertEqual(a[174].getIntensity(), b[174].getIntensity())
        self.assertNotEqual(a[175].getIntensity(), b[175].getIntensity())
        self.assertNotEqual(a[201].getIntensity(), b[201].getIntensity())
        self.assertEqual(a[202].getIntensity(), b[202].getIntensity())

        # All peaks between 174 and 202 should be zero
        for i in range(175, 202, 1):
            self.assertEqual(b[i].getIntensity(), 0)
            self.assertEqual(a[i].getMZ(), b[i].getMZ())

    def test_spectrum_ParentPeakMower(self):
        thisfilter = pyopenms.ParentPeakMower();

        new_firstspec = self.exp[1]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[1])

        # TODO
        self.assertEqual(new_firstspec[10].getMZ(), self.exp[1][10].getMZ())
        self.assertEqual(new_firstspec[10].getIntensity(), self.exp[1][10].getIntensity())

    def test_map_Scaler(self):
        thisfilter = pyopenms.Scaler();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_Scaler(self):
        thisfilter = pyopenms.Scaler();

        new_firstspec = self.exp[0]
        thisfilter.filterSpectrum(new_firstspec)

        self.assertNotEqual(new_firstspec.size(), 0)
        self.assertNotEqual(new_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertNotEqual(new_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(new_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_map_SqrtMower(self):
        thisfilter = pyopenms.SqrtMower();

        old_firstspec = self.exp[0]
        thisfilter.filterPeakMap(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # in most cases, a different spectrum is returned
        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

    def test_spectrum_SqrtMower(self):
        thisfilter = pyopenms.SqrtMower();

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
