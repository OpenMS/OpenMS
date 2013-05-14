import unittest
import os,sys

import pdb
import pyopenms

"""

python pyTOPP/OpenSwathChromatogramExtractor.py --in ../source/TEST/TOPP/OpenSwathChromatogramExtractor_input.mzML\
        --tr ../source/TEST/TOPP/OpenSwathChromatogramExtractor_input.TraML --out OpenSwathChromatogramExtractor.tmp.out \
        && diff OpenSwathChromatogramExtractor.tmp.out ../source/TEST/TOPP/OpenSwathChromatogramExtractor_output.mzML

"""

def main(options):

    # load TraML file
    targeted = pyopenms.TargetedExperiment();
    pyopenms.TraMLFile().load(options.traml_in, targeted);

    # Create empty files as input and finally as output
    empty_swath = pyopenms.MSExperiment()
    trafo = pyopenms.TransformationDescription()
    output = pyopenms.MSExperiment();

    # load input
    for infile in options.infiles:
        exp = pyopenms.MSExperiment()
        pyopenms.FileHandler().loadExperiment(infile, exp)

        transition_exp_used = pyopenms.TargetedExperiment();

        do_continue = True
        if options.is_swath:
            do_continue = pyopenms.OpenSwathHelper().checkSwathMapAndSelectTransitions(exp, targeted, transition_exp_used, options.min_upper_edge_dist)
        else:
            transition_exp_used = targeted

        if do_continue:
            # set up extractor and run
            tmp_out = pyopenms.MSExperiment();
            extractor = pyopenms.ChromatogramExtractor()
            extractor.extractChromatograms(exp, tmp_out, targeted, options.extraction_window, options.ppm, trafo, options.rt_extraction_window, options.extraction_function)
            # add all chromatograms to the output
            for chrom in tmp_out.getChromatograms():
                output.addChromatogram(chrom)

    dp = pyopenms.DataProcessing()
    pa = pyopenms.ProcessingAction().SMOOTHING
    dp.setProcessingActions(set([pa]))

    chromatograms = output.getChromatograms();
    for chrom in chromatograms:
        this_dp = chrom.getDataProcessing()
        this_dp.append(dp)
        chrom.setDataProcessing(this_dp)

    output.setChromatograms(chromatograms);

    pyopenms.MzMLFile().store(options.outfile, output);

def handle_args():
    import argparse

    usage = "" 
    usage += "\nExtract chromatograms (XIC) from a MS2 map file."

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infiles", nargs = '+', help = 'An input file containing spectra')
    parser.add_argument("--tr", dest="traml_in", help="TraML input file containt the transitions")
    parser.add_argument("--out", dest="outfile", help="Output chrom.mzML file with chromatograms")
    parser.add_argument("--extraction_window", dest="extraction_window", default=0.05, help="Extraction window in Th", metavar='0.05', type=float)
    parser.add_argument("--min_upper_edge_dist", dest="min_upper_edge_dist", default=0.0, help="Minimal distance to the edge to still consider a precursor, in Thomson (for Swath)", metavar='0.0', type=float)
    parser.add_argument('--ppm', action='store_true', default=False, help="use ppm instead of Th")
    parser.add_argument('--is_swath', action='store_true', default=False, help="The input file is a SWATH file")
    parser.add_argument("--rt_extraction_window", dest="rt_extraction_window", default=-1, help="Extraction window in RT", metavar='-1', type=float)
    parser.add_argument("--extraction_function", dest="extraction_function", default="tophat", help="Extraction function (tophat or bartlett)", metavar="tophat")

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)

