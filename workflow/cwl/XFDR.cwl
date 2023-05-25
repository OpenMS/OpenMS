inputs:
  in:
    doc: Crosslink Identifications in either xquest.xml, idXML, or mzIdentML format (as produced by OpenPepXL)
    type: File?
  in_type:
    doc: Type of input file provided with -in. If omitted, the file type is guessed from the file extension.
    type: string?
  out_idXML:
    doc: Output as idXML file
    type: string
  out_mzIdentML:
    doc: Output as mzIdentML file
    type: string
  out_xquest:
    doc: Output as xquest.xml file
    type: string
  decoy_string:
    doc: Prefix of decoy protein ids. The correspondig target protein id should be retrievable by deleting this prefix.
    type: string?
  minborder:
    doc: Filter for minimum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.
    type: double?
  maxborder:
    doc: Filter for maximum precursor mass error (ppm) before FDR estimation. Values outside of the tolerance window of the original search will effectively disable this filter.
    type: double?
  mindeltas:
    doc: Filter for delta score, 0 disables the filter. Minimum delta score required, hits are rejected if larger or equal. The delta score is a ratio of the score of a hit and the score of the next best hit to the same spectrum, so the value range is between 0 and 1 with 1.0 meaning the scores are equal and 0.5 meaning the next best score is half as high as the current one.
    type: double?
  minionsmatched:
    doc: Filter for minimum matched ions per peptide.
    type: long?
  uniquexl:
    doc: Calculate statistics based only on unique IDs. For a set of IDs from equal candidates (same pair of peptides, modifications and cross-linked positions), only the highest scoring hit will be considered. By default the score distribution will be estimated using all 1st ranked candidates.
    type: boolean?
  no_qvalues:
    doc: Do not transform simple FDR to q-values
    type: boolean?
  minscore:
    doc: Minimum score to be considered for FDR calculation. A number lower than the lowest score will effectively disable this filter.
    type: double?
  binsize:
    doc: Bin size for the cumulative histograms for score distributions. Should be about the same size as the smallest expected difference between scores. Smaller numbers will make XFDR more robust, but much slower. Negative numbers are not allowed. Should only be changed if the range of the main score changes or another score than the OpenPepXL score is used.
    type: double?
  log:
    doc: Name of log file (created only when specified)
    type: string?
  debug:
    doc: Sets the debug level
    type: long?
  threads:
    doc: Sets the number of threads allowed to be used by the TOPP tool
    type: long?
  no_progress:
    doc: Disables progress logging to command line
    type: boolean?
  force:
    doc: Overrides tool-specific checks
    type: boolean?
outputs:
  out_idXML:
    type: File
    outputBinding:
      glob: $(inputs.out_idXML)
  out_mzIdentML:
    type: File
    outputBinding:
      glob: $(inputs.out_mzIdentML)
  out_xquest:
    type: File
    outputBinding:
      glob: $(inputs.out_xquest)
label: XFDR
doc: Calculates false discovery rate estimates on crosslink identifications
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - XFDR
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json