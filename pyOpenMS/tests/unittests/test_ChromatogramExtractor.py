import unittest
import os

import pyopenms

class TestChromatogramExtractor(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.TraML")
        self.filename_mzml = os.path.join(dirname, "test2.mzML")

    def test_extractor(self):
        targeted = pyopenms.TargetedExperiment();
        tramlfile = pyopenms.TraMLFile();
        tramlfile.load(self.filename, targeted);

        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename_mzml, exp)

        trafo = pyopenms.TransformationDescription()

        tmp_out = pyopenms.MSExperiment();
        extractor = pyopenms.ChromatogramExtractor()
        extractor.extractChromatograms(exp, tmp_out, targeted, 10, False, trafo, -1, "tophat")

        # Basically test that the output is non-zero (e.g. the data is
        # correctly relayed to python)
        # The functionality is not tested here!
        self.assertEqual(len(tmp_out.getChromatograms()), len(targeted.getTransitions()))
        self.assertNotEqual(len(tmp_out.getChromatograms()), 0)
        self.assertEqual(tmp_out.getChromatograms()[0].size(), exp.size())
        self.assertNotEqual(tmp_out.getChromatograms()[0].size(), 0)
        self.assertNotEqual(tmp_out.getChromatograms()[0][0].getRT(), 0)
        self.assertNotEqual(tmp_out.getChromatograms()[0][0].getIntensity(), 0)

if __name__ == '__main__':
    unittest.main()


