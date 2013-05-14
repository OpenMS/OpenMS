import unittest,os,copy

import env
import pyopenms
from collections import defaultdict

eps = 2

class TestSILACAnalyzer(unittest.TestCase):
    """Emulates the behavior of SILACAnalyzer"""

    def setUp(self):
        self.testdirname = os.path.join(env.OPEN_MS_SRC, "source/TEST/TOPP")
        # set up files
        self.infile = os.path.join(self.testdirname, "SILACAnalyzer.mzML")

    def test_run_SILACAnalyzer(self):

        exp = pyopenms.MSExperiment()
        out_map = pyopenms.ConsensusMap()
        pyopenms.FileHandler().loadExperiment(self.infile, exp)
        exp.updateRanges()

        # 
        # 1. filter MS1 level (only keep MS1)
        # 
        tmp = copy.copy(exp)
        for spectrum in exp:
            if spectrum.getMSLevel() == 1:
                tmp.push_back(spectrum)
        exp = tmp
        exp.sortSpectra(True)

        selected_labels = "Lys8"
        charge_min = 2
        charge_max = 2
        missed_cleavages = 0
        isotopes_per_peptide_min = 3
        isotopes_per_peptide_max = 3

        rt_threshold = 80
        rt_min = 0
        intensity_cutoff = 10
        intensity_correlation = 0.95
        model_deviation = 10
        allow_missing_peaks = False
        label_identifiers = {"Lys8" : 8.0141988132}

        # 
        # 2. set parameters
        # 
        analyzer = pyopenms.SILACAnalyzer()
        analyzer.initialize(
              # section sample
              selected_labels,
              charge_min,
              charge_max,
              missed_cleavages,
              isotopes_per_peptide_min,
              isotopes_per_peptide_max,
              # section "algorithm"
              rt_threshold,
              rt_min,
              intensity_cutoff,
              intensity_correlation,
              model_deviation,
              allow_missing_peaks,
              # labels
              label_identifiers)

        # 
        # 3. run
        # 
        analyzer.run_all(exp, out_map)

        out_map.sortByPosition()

        self.assertEqual(out_map.size(), 3)
        self.assertEqual(out_map[0].getQuality(), 8.0)
        self.assertEqual(out_map[1].getQuality(), 8.0)
        self.assertEqual(out_map[2].getQuality(), 8.0)

        self.assertAlmostEqual(out_map[0].getRT(), 6632.409179688, eps)
        self.assertAlmostEqual(out_map[1].getRT(), 6635.169433594, eps)
        self.assertAlmostEqual(out_map[2].getRT(), 6657.56445312, eps)

        self.assertAlmostEqual(out_map[0].getMZ(), 668.321350097656, eps)
        self.assertAlmostEqual(out_map[1].getMZ(), 670.894470214844, eps)
        self.assertAlmostEqual(out_map[2].getMZ(), 668.8262329102, eps)

if __name__ == '__main__':
    unittest.main()
