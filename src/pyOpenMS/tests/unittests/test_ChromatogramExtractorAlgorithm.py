import unittest
import os

import pyopenms

class TestChromatogramExtractorAlgorithm(unittest.TestCase):

    def setUp(self):
        dirname = os.path.dirname(os.path.abspath(__file__))
        self.filename = os.path.join(dirname, "test2.mzML").encode()

    def test_readfile_content(self):
        exp = pyopenms.MSExperiment()
        pyopenms.MzMLFile().load(self.filename, exp)
        exp_size = exp.size()
        saccess = pyopenms.SpectrumAccessOpenMS(exp)

        ### double mz # mz around which should be extracted
        ### double rt_start # rt start of extraction (in seconds)
        ### double rt_end # rt end of extraction (in seconds)
        ### libcpp_string id # identifier
        targeted = []
        coord = pyopenms.ExtractionCoordinates()
        coord.mz = 618.31
        coord.rt_start = 4000 
        coord.rt_end = 5000
        coord.id = b"tr3"
        targeted.append(coord)

        coord = pyopenms.ExtractionCoordinates()
        coord.mz = 628.45
        coord.rt_start = 4000 
        coord.rt_end = 5000
        coord.id = b"tr1"
        targeted.append(coord)

        coord = pyopenms.ExtractionCoordinates()
        coord.mz = 654.38
        coord.rt_start = 4000 
        coord.rt_end = 5000
        coord.id = b"tr2"
        targeted.append(coord)

        trafo = pyopenms.TransformationDescription()

        # Start with length zero
        tmp_out = [ pyopenms.OSChromatogram() for i in range(len(targeted))]
        self.assertEqual(len(tmp_out[0].getIntensityArray()), 0)

        extractor = pyopenms.ChromatogramExtractorAlgorithm()
        mz_extraction_window = 10.0
        ppm = False
        extractor.extractChromatograms(saccess, tmp_out, targeted, mz_extraction_window, ppm, -1.0, b"tophat")

        # Basically test that the output is non-zero (e.g. the data is
        # correctly relayed to python)
        # The functionality is not tested here!
        self.assertEqual(len(tmp_out), len(targeted))
        self.assertNotEqual(len(tmp_out), 0)

        # End with different length
        self.assertEqual(len(tmp_out[0].getIntensityArray()), exp_size)
        self.assertNotEqual(len(tmp_out[0].getIntensityArray()), 0)
        self.assertNotEqual(len(tmp_out[0].getTimeArray()), 0)
        self.assertNotEqual(tmp_out[0].getIntensityArray()[0], 0)
        self.assertNotEqual(tmp_out[0].getTimeArray()[0], 0)

if __name__ == '__main__':
    unittest.main()
