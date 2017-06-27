
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


def algorithm(exp, targeted, picker, scorer):

    output = pyopenms.FeatureMap()

    trafo = pyopenms.TransformationDescription()

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
        transition_group = getTransitionGroup(exp, targeted, key, value, chrom_map)
        picker.pickTransitionGroup(transition_group);
        scorer.scorePeakgroups(transition_group, trafo, swath_maps_dummy, output, False);

    return output

def main(options):
    out = options.outfile
    chromat_in = options.infile
    traml_in = options.traml_in

    pp = pyopenms.MRMTransitionGroupPicker()

    # scoring_params = pyopenms.MRMFeatureFinderScoring().getDefaults();

    pp_params = pp.getDefaults();
    pp_params.setValue("PeakPickerMRM:remove_overlapping_peaks", options.remove_overlapping_peaks, '')
    pp_params.setValue("PeakPickerMRM:method", options.method, '')
    pp.setParameters(pp_params);

    scorer = pyopenms.MRMFeatureFinderScoring()

    chromatograms = pyopenms.MSExperiment()
    fh = pyopenms.FileHandler()
    fh.loadExperiment(chromat_in, chromatograms)
    targeted = pyopenms.TargetedExperiment();
    tramlfile = pyopenms.TraMLFile();
    tramlfile.load(traml_in, targeted);

    light_targeted = pyopenms.LightTargetedExperiment();
    pyopenms.OpenSwathDataAccessHelper().convertTargetedExp(targeted, light_targeted)
    output = algorithm(chromatograms, light_targeted, pp, scorer)

    pyopenms.FeatureXMLFile().store(out, output);

def handle_args():
    import argparse

    usage = ""
    usage += "\nMRMTransitionGroupScorer picks transition groups in measured chromatograms (mzML) and scores them"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", help = 'An input file containing chromatograms')
    parser.add_argument("--tr", dest="traml_in", help="TraML input file containt the transitions")
    parser.add_argument("--out", dest="outfile", help="Output file with annotated chromatograms")
    parser.add_argument("--remove_overlapping_peaks", dest="remove_overlapping_peaks", default="false", help="true/false", metavar='0.1', type=str)
    parser.add_argument("--method", dest="method", default="legacy", help="legacy/corrected", metavar='0.1', type=str)

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)

