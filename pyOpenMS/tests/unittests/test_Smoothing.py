import unittest
import os

import pyopenms

class TestGaussFilter(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML")
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, self.exp)

    def test_init(self):
        thisfilter = pyopenms.GaussFilter();

    def test_run(self):
        thisfilter = pyopenms.GaussFilter();
        old_firstspec = self.exp[0]
        thisfilter.filterExperiment(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # MZ should not change, Intensity should
        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

class TestSavitzkyGolayFilter(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML")
        self.exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, self.exp)

    def test_init(self):
        thisfilter = pyopenms.SavitzkyGolayFilter();

    def test_run(self):
        thisfilter = pyopenms.SavitzkyGolayFilter();
        old_firstspec = self.exp[0]
        thisfilter.filterExperiment(self.exp)

        self.assertNotEqual(self.exp.size(), 0)
        self.assertNotEqual(old_firstspec, self.exp[0])

        # MZ should not change, Intensity should
        self.assertEqual(old_firstspec[10].getMZ(), self.exp[0][10].getMZ())
        self.assertNotEqual(old_firstspec[10].getIntensity(), self.exp[0][10].getIntensity())

class TestLowessSmoothing(unittest.TestCase):

    def setUp(self):
        pass

    def test_init(self):
        thisfilter = pyopenms.LowessSmoothing();

    def test_init(self):
        thisfilter = pyopenms.LowessSmoothing();
        x = [1.0,2.0,3.0,4.0]
        y = [10.0,11.0,12.0,13.0]
        y_smoothed = [0.0]
        thisfilter.smoothData(x,y,y_smoothed)

        self.assertNotEqual( len(y_smoothed), 0)

        # The smoothed data should be different from the input data
        self.assertNotEqual( y_smoothed, y)
        self.assertNotEqual( y_smoothed[0], y[0])

if __name__ == '__main__':
    unittest.main()
