import unittest,os

import pdb
import pyopenms

def simple_find_best_feature(output, pairs, targeted):
  f_map = {}
  for f in output:
    key = f.getMetaValue("PeptideRef").toString()
    if f_map.has_key(key):
      f_map[key].append(f)
    else: 
      f_map[key] = [f]
  
  
  for v in f_map.values():
    bestscore = -10000
    for feature in v:
      score = feature.getMetaValue("main_var_xx_lda_prelim_score").toDouble()
      if score > bestscore:
        best = feature
        bestscore = score
    
    pep = targeted.getPeptideByRef( feature.getMetaValue("PeptideRef").toString()  )
    pairs.append( [best.getRT(), pep.getRetentionTime() ] )

class TestMRMRTNormalizer(unittest.TestCase):
    """Emulates the behavior of OpenSwathMRMRTNormalizer"""

    def setUp(self):
        self.dirname = os.path.dirname(os.path.abspath(__file__))
        self.testdirname = os.path.join(self.dirname, "../../../source/TEST/TOPP/")
        # set up files
        self.chromatograms = os.path.join(self.testdirname, "OpenSwathRTNormalizer_1_input.mzML")
        self.tramlfile = os.path.join(self.testdirname, "OpenSwathRTNormalizer_1_input.TraML")

    def test_run_mrmrtnormalizer(self):

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
        # set the correct rt use values
        scoring_params = pyopenms.MRMFeatureFinderScoring().getDefaults();
        scoring_params.setValue("Scores:use_rt_score", pyopenms.DataValue('false'), '')
        featurefinder.setParameters(scoring_params);
        featurefinder.pickExperiment(chromatograms, output, targeted, trafo, empty_swath)

        # get the pairs
        pairs=[]
        simple_find_best_feature(output, pairs, targeted)
        pairs_corrected = pyopenms.MRMRTNormalizer().rm_outliers( pairs, 0.95, 0.6) 
        pairs_corrected = [ list(p) for p in pairs_corrected] 

        expected = [(1497.56884765625, 1881.0),
             (2045.9776611328125, 2409.0),
             (2151.4814453125, 2509.0),
             (1924.0750732421875, 2291.0),
             (612.9832153320312, 990.0),
             (1086.2474365234375, 1470.0),
             (1133.89404296875, 1519.0),
             (799.5291137695312, 1188.0),
             (1397.1541748046875, 1765.0)]

        for exp,res in zip(expected, pairs_corrected):
            self.assertAlmostEqual(exp[0], res[0]) 
            self.assertAlmostEqual(exp[1], res[1]) 

if __name__ == '__main__':
    unittest.main()
