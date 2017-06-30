
import sys
import pyopenms

def getTransitionGroup(exp, targeted, key, trgr_ids, chrom_map):
    r = pyopenms.MRMTransitionGroupCP()
    for tr in trgr_ids:
        transition = targeted.getTransitions()[tr]
        chrom_idx = chrom_map[ transition.getNativeID() ]
        chrom = exp.getChromatograms()[ chrom_idx ]
        chrom.setMetaValue("product_mz", transition.getProductMZ() )
        chrom.setMetaValue("precursor_mz", transition.getPrecursorMZ() )
        r.addTransition( transition, transition.getNativeID() )
        r.addChromatogram( chrom, chrom.getNativeID() )
    return r


def algorithm(exp, targeted, picker):

    output = pyopenms.FeatureMap()

    chrom_map = {}
    pepmap = {}
    trmap = {}
    for i, chrom in enumerate(exp.getChromatograms()):
        chrom_map[ chrom.getNativeID() ] = i
    for i, pep in enumerate(targeted.getPeptides() ):
        pepmap[ pep.id ] = i
    for i, tr in enumerate(targeted.getTransitions() ):
        tmp = trmap.get( tr.getPeptideRef() , [])
        tmp.append( i )
        trmap[ tr.getPeptideRef() ] = tmp

    for key, value in trmap.iteritems():
        print key, value
        transition_group = getTransitionGroup(exp, targeted, key, value, chrom_map)
        picker.pickTransitionGroup(transition_group);
        for mrmfeature in transition_group.getFeatures():
            features = mrmfeature.getFeatures()
            for f in features:
                # TODO
                # f.getConvexHulls().clear()
                f.ensureUniqueId()

            mrmfeature.setSubordinates(features) # add all the subfeatures as subordinates
            output.push_back(mrmfeature)

    return output

def main(options):
    out = options.outfile
    chromat_in = options.infile
    traml_in = options.traml_in

    pp = pyopenms.MRMTransitionGroupPicker()

    pp_params = pp.getDefaults();
    pp_params.setValue("PeakPickerMRM:remove_overlapping_peaks", options.remove_overlapping_peaks, '')
    pp_params.setValue("PeakPickerMRM:method", options.method, '')
    pp.setParameters(pp_params);

    chromatograms = pyopenms.MSExperiment()
    fh = pyopenms.FileHandler()
    fh.loadExperiment(chromat_in, chromatograms)
    targeted = pyopenms.TargetedExperiment();
    tramlfile = pyopenms.TraMLFile();
    tramlfile.load(traml_in, targeted);

    output = algorithm(chromatograms, targeted, pp)

    pyopenms.FeatureXMLFile().store(out, output);

def handle_args():
    import argparse

    usage = ""
    usage += "\nMRMTransitionGroupPicker picks transition groups in measured chromatograms (mzML)"

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

