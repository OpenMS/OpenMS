
import sys
import pyopenms

def getTransitionGroup(exp, targeted, key, trgr_ids, chrom_map):
    r = pyopenms.LightMRMTransitionGroupCP()
    r.setTransitionGroupID(key)
    for tr in trgr_ids:
        transition = targeted.getTransitions()[tr]
        chrom_idx = chrom_map[ transition.getNativeID() ]
        chrom = exp.getChromatograms()[ chrom_idx ]
        chrom.setMetaValue("product_mz", transition.getProductMZ() )
        chrom.setMetaValue("precursor_mz", transition.getPrecursorMZ() )
        r.addTransition( transition, transition.getNativeID() )
        r.addChromatogram( chrom, chrom.getNativeID() )
    return r


def algorithm(exp, targeted, picker, scorer, trafo):

    output = pyopenms.FeatureMap()

    scorer.prepareProteinPeptideMaps_(targeted)

    chrom_map = {}
    pepmap = {}
    trmap = {}
    for i, chrom in enumerate(exp.getChromatograms()):
        chrom_map[ chrom.getNativeID() ] = i
    for i, pep in enumerate(targeted.getCompounds() ):
        pepmap[ pep.id ] = i
    for i, tr in enumerate(targeted.getTransitions() ):
        tmp = trmap.get( tr.getPeptideRef() , [])
        tmp.append( i )
        trmap[ tr.getPeptideRef() ] = tmp

    swath_maps_dummy = []
    for key, value in trmap.iteritems():
        try:
            transition_group = getTransitionGroup(exp, targeted, key, value, chrom_map)
        except Exception:
            print "Skip ", key, value
            continue
        picker.pickTransitionGroup(transition_group);
        scorer.scorePeakgroups(transition_group, trafo, swath_maps_dummy, output, False);

    return output

def main(options):
    out = options.outfile
    chromat_in = options.infile
    traml_in = options.traml_in
    trafo_in = options.trafo_in

    pp = pyopenms.MRMTransitionGroupPicker()


    metabolomics = False
    # this is an important weight for RT-deviation -- the larger the value, the less importance will be given to exact RT matches
    # for proteomics data it tends to be a good idea to set it to the length of
    # the RT space (e.g. for 100 second RT space, set it to 100)
    rt_normalization_factor = 100.0

    pp_params = pp.getDefaults();
    pp_params.setValue("PeakPickerMRM:remove_overlapping_peaks", options.remove_overlapping_peaks, '')
    pp_params.setValue("PeakPickerMRM:method", options.method, '')
    if (metabolomics):
        # Need to change those for metabolomics and very short peaks!
        pp_params.setValue("PeakPickerMRM:signal_to_noise", 0.01, '')
        pp_params.setValue("PeakPickerMRM:peak_width", 0.1, '')
        pp_params.setValue("PeakPickerMRM:gauss_width", 0.1, '')
        pp_params.setValue("resample_boundary", 0.05, '')
        pp_params.setValue("compute_peak_quality", "true", '')
    pp.setParameters(pp_params)

    scorer = pyopenms.MRMFeatureFinderScoring()
    scoring_params = scorer.getDefaults();
    # Only report the top 5 features
    scoring_params.setValue("stop_report_after_feature", 5, '')
    scoring_params.setValue("rt_normalization_factor", rt_normalization_factor, '')
    scorer.setParameters(scoring_params);

    chromatograms = pyopenms.MSExperiment()
    fh = pyopenms.FileHandler()
    fh.loadExperiment(chromat_in, chromatograms)
    targeted = pyopenms.TargetedExperiment();
    tramlfile = pyopenms.TraMLFile();
    tramlfile.load(traml_in, targeted);

    trafoxml = pyopenms.TransformationXMLFile()
    trafo = pyopenms.TransformationDescription()
    if trafo_in is not None:
        model_params = pyopenms.Param()
        model_params.setValue("symmetric_regression", "false", "", [])
        model_type = "linear"
        trafoxml.load(trafo_in, trafo, True)
        trafo.fitModel(model_type, model_params);


    light_targeted = pyopenms.LightTargetedExperiment();
    pyopenms.OpenSwathDataAccessHelper().convertTargetedExp(targeted, light_targeted)
    output = algorithm(chromatograms, light_targeted, pp, scorer, trafo)

    pyopenms.FeatureXMLFile().store(out, output);

def handle_args():
    import argparse

    usage = ""
    usage += "\nMRMTransitionGroupScorer picks transition groups in measured chromatograms (mzML) and scores them"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", help = 'An input file containing chromatograms')
    parser.add_argument("--tr", dest="traml_in", help="TraML input file containt the transitions")
    parser.add_argument("--trafo", dest="trafo_in", help="Trafo input file")
    parser.add_argument("--out", dest="outfile", help="Output file with annotated chromatograms")
    parser.add_argument("--remove_overlapping_peaks", dest="remove_overlapping_peaks", default="false", help="true/false", metavar='0.1', type=str)
    parser.add_argument("--method", dest="method", default="legacy", help="legacy/corrected", metavar='0.1', type=str)

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)

