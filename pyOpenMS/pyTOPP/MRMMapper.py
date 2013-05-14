import copy, sys
import pyopenms

"""
python pyTOPP/MRMMapper.py --in ../source/TEST/TOPP/MRMMapping_input.chrom.mzML  \
        --tr ../source/TEST/TOPP/MRMMapping_input.TraML --out MRMMapper.tmp.out \
        && diff MRMMapper.tmp.out ../source/TEST/TOPP/MRMMapping_output.chrom.mzML 

"""

def algorithm(chromatogram_map, targeted, precursor_tolerance, product_tolerance, allow_unmapped=True, allow_double_mappings=False):

    # copy all meta data from old chromatogram
    # TODO how to copy this! improve this!
    output = pyopenms.MSExperiment();
    output.fromExperiment(chromatogram_map)
    output.clear(False); 
    empty_chromats = []
    output.setChromatograms(empty_chromats);

    notmapped = 0
    for chrom in chromatogram_map.getChromatograms():
        mapped_already = False
        for transition in targeted.getTransitions():
            if (abs(chrom.getPrecursor().getMZ() - transition.getPrecursorMZ()) < precursor_tolerance and
                abs(chrom.getProduct().getMZ()  -  transition.getProductMZ()) < product_tolerance):
                if mapped_already:
                    this_peptide = targeted.getPeptideByRef(transition.getPeptideRef() ).sequence
                    other_peptide = chrom.getPrecursor().getMetaValue("peptide_sequence").toString()
                    print "Found mapping of", chrom.getPrecursor().getMZ(), "/", chrom.getProduct().getMZ(), "to", transition.getPrecursorMZ(), "/",transition.getProductMZ()
                    print "Of peptide", this_peptide
                    print "But the chromatogram is already mapped to", other_peptide
                    if not allow_double_mappings: raise Exception("Cannot map twice")
                mapped_already = True
                precursor = chrom.getPrecursor();
                peptide = targeted.getPeptideByRef(transition.getPeptideRef() )
                precursor.setMetaValue("peptide_sequence", pyopenms.DataValue(peptide.sequence) )
                chrom.setPrecursor(precursor)
                chrom.setNativeID(transition.getNativeID())
        if not mapped_already:
            notmapped += 1
            print "Did not find a mapping for chromatogram", chrom.getNativeID()
            if not allow_unmapped: raise Exception("No mapping")
        else:
            output.addChromatogram(chrom)

    if notmapped > 0:
        print "Could not find mapping for", notmapped, "chromatogram(s)" 


    dp = pyopenms.DataProcessing()
    # dp.setProcessingActions(ProcessingAction:::FORMAT_CONVERSION)
    pa = pyopenms.ProcessingAction().FORMAT_CONVERSION
    dp.setProcessingActions(set([pa]))

    chromatograms = output.getChromatograms();
    for chrom in chromatograms:
        this_dp = chrom.getDataProcessing()
        this_dp.append(dp)
        chrom.setDataProcessing(this_dp)

    output.setChromatograms(chromatograms);
    return output

def main(options):
    precursor_tolerance = options.precursor_tolerance
    product_tolerance = options.product_tolerance
    out = options.outfile
    chromat_in = options.infile
    traml_in = options.traml_in

    # precursor_tolerance = 0.05
    # product_tolerance = 0.05
    # out = "/tmp/out.mzML"
    # chromat_in = "../source/TEST/TOPP/MRMMapping_input.chrom.mzML"
    # traml_in = "../source/TEST/TOPP/MRMMapping_input.TraML"

    ff = pyopenms.MRMFeatureFinderScoring()
    chromatogram_map = pyopenms.MSExperiment()
    fh = pyopenms.FileHandler()
    fh.loadExperiment(chromat_in, chromatogram_map)
    targeted = pyopenms.TargetedExperiment();
    tramlfile = pyopenms.TraMLFile();
    tramlfile.load(traml_in, targeted);
     
    output = algorithm(chromatogram_map, targeted, precursor_tolerance, product_tolerance)

    pyopenms.MzMLFile().store(out, output);

def handle_args():
    import argparse

    usage = "" 
    usage += "\nMRMMapper maps measured chromatograms (mzML) and the transitions used (TraML)"

    parser = argparse.ArgumentParser(description = usage )
    parser.add_argument('--in', dest="infile", help = 'An input file containing chromatograms')
    parser.add_argument("--tr", dest="traml_in", help="TraML input file containt the transitions")
    parser.add_argument("--out", dest="outfile", help="Output file with annotated chromatograms")
    parser.add_argument("--precursor_tolerance", dest="precursor_tolerance", default=0.1, help="Precursor tolerance when mapping (in Th)", metavar='0.1', type=float)
    parser.add_argument("--product_tolerance", dest="product_tolerance", default=0.1, help="Product tolerance when mapping (in Th)", metavar='0.1', type=float)

    args = parser.parse_args(sys.argv[1:])
    return args

if __name__ == '__main__':
    options = handle_args()
    main(options)
