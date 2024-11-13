import unittest
import os

import pyopenms

class TestChromatogramExtractor(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test.TraML").encode()
        self.filename_mzml = os.path.join(dirname, "test2.mzML").encode()

    def test_extractor(self):
        targeted = pyopenms.TargetedExperiment()
        tramlfile = pyopenms.TraMLFile()
        tramlfile.load(self.filename, targeted)

        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename_mzml, exp)

        trafo = pyopenms.TransformationDescription()
        # Prepare extraction coordinates
        output_chromatograms = []
        extraction_coordinates = []
        pyopenms.ChromatogramExtractor.prepare_coordinates(output_chromatograms, 
                                                         extraction_coordinates,
                                                         targeted,
                                                         -1, # rt_extraction_window
                                                         False, # ms1
                                                         0) # ms1_isotopes

        # Create input map
        input_map = pyopenms.SpectrumAccessOpenMS(exp)

        # Extract chromatograms
        extractor = pyopenms.ChromatogramExtractor()
        # extract full RT range (-1)
        extractor.extractChromatograms(input_map, output_chromatograms, extraction_coordinates, 10, False, -1, b"tophat")

        # Basically test that the output is non-zero (e.g. the data is
        # correctly relayed to python)
        # The functionality is not tested here!
        self.assertEqual(len(output_chromatograms), len(targeted.getTransitions()))
        self.assertNotEqual(len(output_chromatograms), 0)
        self.assertNotEqual(len(output_chromatograms[0].getIntensityArray()), 0)
        self.assertNotEqual(len(output_chromatograms[0].getTimeArray()), 0)
        # one chromatogram per transition; one chromatographic peak for each spectrum
        self.assertEqual(len(output_chromatograms[0].getIntensityArray()), exp.size())
        self.assertEqual(len(output_chromatograms[0].getTimeArray()), exp.size())

if __name__ == '__main__':
    unittest.main()


