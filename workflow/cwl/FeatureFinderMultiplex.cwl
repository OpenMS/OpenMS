inputs:
  in:
    doc: LC-MS dataset in either centroid or profile mode
    type: File
  out:
    doc: Output file containing the individual peptide features.
    type: string
  out_multiplets:
    doc: Optional output file containing all detected peptide groups (i.e. peptide pairs or triplets or singlets or ..). The m/z-RT positions correspond to the lightest peptide in each group.
    type: string
  out_blacklist:
    doc: Optional output file containing all peaks which have been associated with a peptide feature (and subsequently blacklisted).
    type: string
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
  algorithm__labels:
    doc: "Labels used for labelling the samples. If the sample is unlabelled (i.e. you want to detect only single peptide features) please leave this parameter empty. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL"
    type: string?
  algorithm__charge:
    doc: "Range of charge states in the sample, i.e. min charge : max charge."
    type: string?
  algorithm__isotopes_per_peptide:
    doc: "Range of isotopes per peptide in the sample. For example 3:6, if isotopic peptide patterns in the sample consist of either three, four, five or six isotopic peaks. "
    type: string?
  algorithm__rt_typical:
    doc: Typical retention time [s] over which a characteristic peptide elutes. (This is not an upper bound. Peptides that elute for longer will be reported.)
    type: double?
  algorithm__rt_band:
    doc: The algorithm searches for characteristic isotopic peak patterns, spectrum by spectrum. For some low-intensity peptides, an important peak might be missing in one spectrum but be present in one of the neighbouring ones. The algorithm takes a bundle of neighbouring spectra with width rt_band into account. For example with rt_band = 0, all characteristic isotopic peaks have to be present in one and the same spectrum. As rt_band increases, the sensitivity of the algorithm but also the likelihood of false detections increases.
    type: double?
  algorithm__rt_min:
    doc: Lower bound for the retention time [s]. (Any peptides seen for a shorter time period are not reported.)
    type: double?
  algorithm__mz_tolerance:
    doc: m/z tolerance for search of peak patterns.
    type: double?
  algorithm__mz_unit:
    doc: Unit of the 'mz_tolerance' parameter.
    type: string?
  algorithm__intensity_cutoff:
    doc: Lower bound for the intensity of isotopic peaks.
    type: double?
  algorithm__peptide_similarity:
    doc: Two peptides in a multiplet are expected to have the same isotopic pattern. This parameter is a lower bound on their similarity.
    type: double?
  algorithm__averagine_similarity:
    doc: The isotopic pattern of a peptide should resemble the averagine model at this m/z position. This parameter is a lower bound on similarity between measured isotopic pattern and the averagine model.
    type: double?
  algorithm__averagine_similarity_scaling:
    doc: Let x denote this scaling factor, and p the averagine similarity parameter. For the detection of single peptides, the averagine parameter p is replaced by p' = p + x(1-p), i.e. x = 0 -> p' = p and x = 1 -> p' = 1. (For knock_out = true, peptide doublets and singlets are detected simultaneously. For singlets, the peptide similarity filter is irreleavant. In order to compensate for this 'missing filter', the averagine parameter p is replaced by the more restrictive p' when searching for singlets.)
    type: double?
  algorithm__missed_cleavages:
    doc: Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)
    type: long?
  algorithm__spectrum_type:
    doc: Type of MS1 spectra in input mzML file. 'automatic' determines the spectrum type directly from the input mzML file.
    type: string?
  algorithm__averagine_type:
    doc: The type of averagine to use, currently RNA, DNA or peptide
    type: string?
  algorithm__knock_out:
    doc: Is it likely that knock-outs are present? (Supported for doublex, triplex and quadruplex experiments only.)
    type: boolean?
  labels__Arg6:
    doc: "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188"
    type: double?
  labels__Arg10:
    doc: "Label:13C(6)15N(4)  |  C(-6) 13C(6) N(-4) 15N(4)  |  unimod #267"
    type: double?
  labels__Lys4:
    doc: "Label:2H(4)  |  H(-4) 2H(4)  |  unimod #481"
    type: double?
  labels__Lys6:
    doc: "Label:13C(6)  |  C(-6) 13C(6)  |  unimod #188"
    type: double?
  labels__Lys8:
    doc: "Label:13C(6)15N(2)  |  C(-6) 13C(6) N(-2) 15N(2)  |  unimod #259"
    type: double?
  labels__Leu3:
    doc: "Label:2H(3)  |  H(-3) 2H(3)  |  unimod #262"
    type: double?
  labels__Dimethyl0:
    doc: "Dimethyl  |  H(4) C(2)  |  unimod #36"
    type: double?
  labels__Dimethyl4:
    doc: "Dimethyl:2H(4)  |  2H(4) C(2)  |  unimod #199"
    type: double?
  labels__Dimethyl6:
    doc: "Dimethyl:2H(4)13C(2)  |  2H(4) 13C(2)  |  unimod #510"
    type: double?
  labels__Dimethyl8:
    doc: "Dimethyl:2H(6)13C(2)  |  H(-2) 2H(6) 13C(2)  |  unimod #330"
    type: double?
  labels__ICPL0:
    doc: "ICPL  |  H(3) C(6) N O  |  unimod #365"
    type: double?
  labels__ICPL4:
    doc: "ICPL:2H(4)  |  H(-1) 2H(4) C(6) N O  |  unimod #687"
    type: double?
  labels__ICPL6:
    doc: "ICPL:13C(6)  |  H(3) 13C(6) N O  |  unimod #364"
    type: double?
  labels__ICPL10:
    doc: "ICPL:13C(6)2H(4)  |  H(-1) 2H(4) 13C(6) N O  |  unimod #866"
    type: double?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_multiplets:
    type: File
    outputBinding:
      glob: $(inputs.out_multiplets)
  out_blacklist:
    type: File
    outputBinding:
      glob: $(inputs.out_blacklist)
label: FeatureFinderMultiplex
doc: Determination of peak ratios in LC-MS data
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderMultiplex
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json