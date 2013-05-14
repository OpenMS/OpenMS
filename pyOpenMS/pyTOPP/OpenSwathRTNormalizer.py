import unittest
import os,sys

import pdb
import pyopenms

"""
python pyTOPP/OpenSwathRTNormalizer.py --in ../source/TEST/TOPP/OpenSwathRTNormalizer_1_input.mzML \
        --tr ../source/TEST/TOPP/OpenSwathRTNormalizer_1_input.TraML --out OpenSwathRTNormalizer.tmp.out \
        && diff OpenSwathRTNormalizer.tmp.out ../source/TEST/TOPP/OpenSwathRTNormalizer_1_output.trafoXML 

"""

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

def main(options):

    # load chromatograms
    chromatograms = pyopenms.MSExperiment()
    fh = pyopenms.FileHandler()
    fh.loadExperiment(options.infile, chromatograms)

    # load TraML file
    targeted = pyopenms.TargetedExperiment();
    tramlfile = pyopenms.TraMLFile();
    tramlfile.load(options.traml_in, targeted);

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

    # // store transformation, using a linear model as default
    trafo_out = pyopenms.TransformationDescription()
    trafo_out.setDataPoints(pairs_corrected);
    model_params = pyopenms.Param()
    model_params.setValue("symmetric_regression", pyopenms.DataValue('false'), '');
    model_type = "linear";
    trafo_out.fitModel(model_type, model_params);
    pyopenms.TransformationXMLFile().store(options.outfile, trafo_out);

def handle_args():
    import argparse

    usage = "" 
    usage += "\nMRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", help = 'An input file containing chromatograms')
    parser.add_argument("--tr", dest="traml_in", help="TraML input file containt the transitions")
    parser.add_argument("--out", dest="outfile", help="Output trafoXML file with annotated chromatograms")

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)
