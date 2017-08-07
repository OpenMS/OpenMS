import unittest
import os

import pyopenms

import env

eps = 2

class TestMRMFeatureFinderScoring(unittest.TestCase):

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.testdirname = os.path.join(env.OPEN_MS_SRC, "src/tests/topp")
        # set up files
        self.chromatograms = os.path.join(self.testdirname, "OpenSwathAnalyzer_1_input_chrom.mzML").encode()
        self.tramlfile = os.path.join(self.testdirname, "OpenSwathAnalyzer_1_input.TraML").encode()

    def test_run_mrmfeaturefinder(self):

        # load chromatograms
        chromatograms = pyopenms.MSExperiment()
        fh = pyopenms.FileHandler()
        fh.loadExperiment(self.chromatograms, chromatograms)

        # load TraML file
        targeted = pyopenms.TargetedExperiment();
        tramlfile = pyopenms.TraMLFile();
        tramlfile.load(self.tramlfile, targeted);

        # Create empty files as input and finally as output
        empty_swath = pyopenms.MSExperiment()
        trafo = pyopenms.TransformationDescription()
        output = pyopenms.FeatureMap();

        # set up featurefinder and run
        featurefinder = pyopenms.MRMFeatureFinderScoring()
        featurefinder.pickExperiment(chromatograms, output, targeted, trafo, empty_swath)

        self.assertAlmostEqual(output.size(), 3)
        self.assertAlmostEqual(output[0].getRT(), 3119.092041015, eps)
        self.assertAlmostEqual(output[0].getIntensity(), 3614.99755859375, eps)
        self.assertAlmostEqual(output[0].getMetaValue(b"var_xcorr_shape_weighted"), 0.997577965259552, eps)
        self.assertAlmostEqual(output[0].getMetaValue(b"sn_ratio"), 86.00413513183594, eps)

if __name__ == '__main__':
    unittest.main()
