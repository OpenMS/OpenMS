import copy, sys
import pyopenms

"""
python pyTOPP/SILACAnalyzer.py --in ../source/TEST/TOPP/SILACAnalyzer.mzML --out /tmp/pySILACout.consensusXML --test --selected_labels "Lys8" --rt_threshold 80 --rt_min 0 --intensity_cutoff 10 --model_deviation 10 --charge_max 2 --isotopes_per_peptide_max 3 --intensity_correlation 0.95 && ../bin/FuzzyDiff "-ini" "../source/TEST/TOPP/FuzzyDiff.ini" "-in1" /tmp/pySILACout.consensusXML "-in2" "../source/TEST/TOPP/SILACAnalyzer.consensusXML" && echo "Success"
"""

def main(options):

    # make sure that the ids are "correct" for the testcase
    date_time = pyopenms.DateTime();
    if options.test:
        date_time.set("1999-12-31 23:59:59");
        pyopenms.UniqueIdGenerator().setSeed(date_time);
    else:
        date_time = pyopenms.DateTime.now();

    exp = pyopenms.MSExperiment()
    out_map = pyopenms.ConsensusMap()
    pyopenms.FileHandler().loadExperiment(options.infile, exp)
    exp.updateRanges()

    # 
    # 1. filter MS1 level (only keep MS1)
    # 
    tmp = copy.copy(exp)
    tmp.clear(False); 
    for spectrum in exp:
        if spectrum.getMSLevel() == 1:
            tmp.push_back(spectrum)
    exp = tmp
    exp.sortSpectra(True)

    # 
    # 2. set parameters
    # 
    analyzer = pyopenms.SILACAnalyzer()
    analyzer.initialize(
          # section sample
          options.selected_labels,
          options.charge_min,
          options.charge_max,
          options.missed_cleavages,
          options.isotopes_per_peptide_min,
          options.isotopes_per_peptide_max,
          # section "algorithm"
          options.rt_threshold,
          options.rt_min,
          options.intensity_cutoff,
          options.intensity_correlation,
          options.model_deviation,
          options.allow_missing_peaks,
          # labels
          options.label_identifiers)

    # 
    # 3. run
    # 
    analyzer.run_all(exp, out_map)

    # 
    # 4. set dataprocessing and output meta information
    # 
    out_map.sortByPosition()

    dp = out_map.getDataProcessing() 
    p = pyopenms.DataProcessing()
    p.setProcessingActions(set([ pyopenms.ProcessingAction().DATA_PROCESSING, 
                                  pyopenms.ProcessingAction().PEAK_PICKING,
                                  pyopenms.ProcessingAction().FILTERING,
                                  pyopenms.ProcessingAction().QUANTITATION]))
    p.setCompletionTime(date_time)

    sw = p.getSoftware()
    sw.setName("SILACAnalyzer")
    if options.test:
        sw.setVersion("version_string")
        p.setSoftware(sw)
        p.setMetaValue("parameter: mode", "test_mode")
    else:
        sw.setVersion("pyTOPP v1.10")
        p.setSoftware(sw)
    dp.append(p)
    out_map.setDataProcessing(dp)

    # 
    # 5. write output
    # 
    analyzer.writeConsensus(pyopenms.String(options.outfile), out_map)

def handle_args():
    import argparse

    usage = "" 
    usage += "\nSILAC Analyzer"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", help = 'An input file containing raw MS data (mzML)')
    parser.add_argument("--out", dest="outfile", help="Output file (consensus XML)")

    parser.add_argument("--rt_threshold", dest="rt_threshold", default=30, help="Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)", metavar='30', type=float)
    parser.add_argument("--rt_min", dest="rt_min", default=0.0, help="Lower bound for the retention time [s].", metavar='0.0', type=float)
    parser.add_argument("--intensity_cutoff", dest="intensity_cutoff", default=1000, help="Lower bound for the intensity of isotopic peaks in a SILAC pattern.", metavar='1000', type=float)
    parser.add_argument("--intensity_correlation", dest="intensity_correlation", default=0.7, help="Lower bound for the Pearson correlation coefficient, which measures how well intensity profiles of different isotopic peaks correlate.", metavar='0.7', type=float)
    parser.add_argument("--model_deviation", dest="model_deviation", default=3.0, help="Upper bound on the factor by which the ratios of observed isotopic peaks are allowed to differ from the ratios of the theoretic averagine model, i.e. ( theoretic_ratio / model_deviation ) < observed_ratio < ( theoretic_ratio * model_deviation ).", metavar='3.0', type=float)
    parser.add_argument('--allow_missing_peaks', action='store_true', default=False, help="")

    parser.add_argument("--selected_labels", dest="selected_labels", help="Labels used for labelling the sample. [...] specifies the labels for a single sample. For example, [Lys4,Arg6][Lys8,Arg10] describes a mixtures of three samples. One of them unlabelled, one labelled with Lys4 and Arg6 and a third one with Lys8 and Arg10. For permitted labels see \'advanced parameters\', section \'labels\'. If left empty the tool identifies singlets, i.e. acts as peptide feature finder.")
    parser.add_argument("--charge_min", dest="charge_min", default=2, help="Minimal charge used in the sample", type=int)
    parser.add_argument("--charge_max", dest="charge_max", default=4, help="Maximal charge used in the sample", type=int)
    parser.add_argument("--missed_cleavages", dest="missed_cleavages", default=0, help="", type=int)
    parser.add_argument("--isotopes_per_peptide_min", dest="isotopes_per_peptide_min", default=3, help="Minimal isotopic peaks to consider", metavar='0.1', type=int)
    parser.add_argument("--isotopes_per_peptide_max", dest="isotopes_per_peptide_max", default=5, help="Maximal isotopic peaks to consider", metavar='0.1', type=int)
    parser.add_argument("--label_identifiers", dest="label_identifiers_", default="Lys8:8.0141988132", help="Name of the labels and their weights (muse correspond to option 'selected_labels'). Multiple labels are separted by comma", metavar='Lys8:8.0141988132,...')

    parser.add_argument('--test', action='store_true', default=False, help="Testing mode")

    args = parser.parse_args(sys.argv[1:])

    args.label_identifiers = {}
    for l in args.label_identifiers_.split(","):
        args.label_identifiers[l.split(":")[0]] = float(l.split(":")[1])

    print args.label_identifiers_
    print args.label_identifiers

    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)




