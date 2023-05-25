inputs:
  sirius_executable:
    doc: The Sirius executable. Provide a full or relative path, or make sure it can be found in your PATH environment.
    type: File?
  in:
    doc: MzML input file(s) used for assay library generation
    type: File[]
  in_id:
    doc: FeatureXML input file(s) containing identification information (e.g. AccurateMassSearch)
    type: File[]
  out:
    doc: Assay library output file
    type: string
  fragment_annotation:
    doc: Fragment annotation method
    type: string?
  ambiguity_resolution_mz_tolerance:
    doc: Mz tolerance for the resolution of identification ambiguity over multiple files
    type: double?
  ambiguity_resolution_mz_tolerance_unit:
    doc: Unit of the ambiguity_resolution_mz_tolerance
    type: string?
  ambiguity_resolution_rt_tolerance:
    doc: RT tolerance in seconds for the resolution of identification ambiguity over multiple files
    type: double?
  total_occurrence_filter:
    doc: Filter compound based on total occurrence in analysed samples
    type: double?
  fragment_annotation_score_threshold:
    doc: Filters annotations based on the explained intensity of the peaks in a spectrum
    type: double?
  decoy_generation:
    doc: Decoys will be generated using the fragmentation tree re-rooting approach. This option does only work in combination with the fragment annotation via Sirius.
    type: boolean?
  decoy_generation_method:
    doc: Uses different methods for decoy generation. Basis for the method is the fragmentation-tree re-rooting approach ('original'). This approach can be extended by using 'resolve_overlap', which will resolve overlapping fragments of the highest intensity fragments chosen, by adding -CH2 mass to the overlapping fragments. 'Add_shift' will add a -CH2 mass shift to the target fragments and use them as additional decoys if fragmentation-tree re-rooting failed. 'Both' combines the extended methods (resolve_overlap, add_shift).
    type: string?
  method:
    doc: Spectrum with the highest precursor intensity or a consensus spectrum is used for assay library construction (if no fragment annotation is used).
    type: string?
  use_exact_mass:
    doc: Use exact mass for precursor and fragment annotations
    type: boolean?
  exclude_ms2_precursor:
    doc: Excludes precursor in ms2 from transition list
    type: boolean?
  precursor_mz_distance:
    doc: Max m/z distance of the precursor entries of two spectra to be merged in [Da].
    type: double?
  precursor_recalibration_window:
    doc: Tolerance window for precursor selection (Annotation of precursor mz and intensity)
    type: double?
  precursor_recalibration_window_unit:
    doc: Unit of the precursor_mz_tolerance_annotation
    type: string?
  consensus_spectrum_precursor_rt_tolerance:
    doc: Tolerance window (left and right) for precursor selection [seconds], for consensus spectrum generation (only available without fragment annotation)
    type: double?
  use_known_unknowns:
    doc: Use features without identification information
    type: boolean?
  min_transitions:
    doc: Minimal number of transitions
    type: long?
  max_transitions:
    doc: Maximal number of transitions
    type: long?
  cosine_similarity_threshold:
    doc: Threshold for cosine similarity of MS2 spectra from the same precursor used in consensus spectrum creation
    type: double?
  transition_threshold:
    doc: Further transitions need at least x% of the maximum intensity (default 5%)
    type: double?
  min_fragment_mz:
    doc: Minimal m/z of a fragment ion choosen as a transition
    type: double?
  max_fragment_mz:
    doc: Maximal m/z of a fragment ion choosen as a transition
    type: double?
  read_sirius_stdout:
    doc: Read and print the standard output and error of the Sirius executable, even if it succeeds.
    type: boolean?
  out_workspace_directory:
    doc: Output directory for SIRIUS workspace
    type: string?
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
  deisotoping__use_deisotoper:
    doc: Use Deisotoper (if no fragment annotation is used)
    type: boolean?
  deisotoping__fragment_tolerance:
    doc: Tolerance used to match isotopic peaks
    type: double?
  deisotoping__fragment_unit:
    doc: Unit of the fragment tolerance
    type: string?
  deisotoping__min_charge:
    doc: The minimum charge considered
    type: long?
  deisotoping__max_charge:
    doc: The maximum charge considered
    type: long?
  deisotoping__min_isopeaks:
    doc: The minimum number of isotopic peaks (at least 2) required for an isotopic cluster
    type: long?
  deisotoping__max_isopeaks:
    doc: The maximum number of isotopic peaks (at least 2) considered for an isotopic cluster
    type: long?
  deisotoping__keep_only_deisotoped:
    doc: Only monoisotopic peaks of fragments with isotopic pattern are retained
    type: boolean?
  deisotoping__annotate_charge:
    doc: Annotate the charge to the peaks
    type: boolean?
  preprocessing__filter_by_num_masstraces:
    doc: Number of mass traces each feature has to have to be included. To use this parameter, setting the feature_only flag is necessary
    type: long?
  preprocessing__precursor_mz_tolerance:
    doc: Tolerance window for precursor selection (Feature selection in regard to the precursor)
    type: double?
  preprocessing__precursor_mz_tolerance_unit:
    doc: Unit of the precursor_mz_tolerance
    type: string?
  preprocessing__precursor_rt_tolerance:
    doc: Tolerance window (left and right) for precursor selection [seconds]
    type: double?
  preprocessing__isotope_pattern_iterations:
    doc: Number of iterations that should be performed to extract the C13 isotope pattern. If no peak is found (C13 distance) the function will abort. Be careful with noisy data - since this can lead to wrong isotope patterns
    type: long?
  preprocessing__feature_only:
    doc: Uses the feature information from in_featureinfo to reduce the search space to MS2 associated with a feature
    type: boolean?
  preprocessing__no_masstrace_info_isotope_pattern:
    doc: Use this flag if the masstrace information from a feature should be discarded and the isotope_pattern_iterations should be used instead
    type: boolean?
  project__maxmz:
    doc: "Just consider compounds with a precursor mz lower or equal\nthis maximum mz. All other compounds in the input file\nare ignored."
    type: long?
  project__processors:
    doc: Number of cpu cores to use. If not specified SIRIUS uses all available cores.
    type: long?
  project__loglevel:
    doc: "Set logging level of the Jobs SIRIUS will execute.\nValid values: SEVERE, WARNING, INFO, FINER, ALL\nDefault: WARNING"
    type: string?
  project__ignore_formula:
    doc: Ignore given molecular formula in internal .ms format, while processing.
    type: boolean?
  project__q:
    doc: Suppress shell output
    type: boolean?
  sirius__ppm_max:
    doc: Maximum allowed mass deviation in ppm for decomposing masses [ppm].
    type: double?
  sirius__ppm_max_ms2:
    doc: "Maximum allowed mass deviation in ppm for decomposing masses in MS2 [ppm].If not specified, the same value as for the MS1 is used. "
    type: double?
  sirius__tree_timeout:
    doc: Time out in seconds per fragmentation tree computations. 0 for an infinite amount of time
    type: long?
  sirius__compound_timeout:
    doc: Maximal computation time in seconds for a single compound. 0 for an infinite amount of time.
    type: long?
  sirius__no_recalibration:
    doc: Disable recalibration of input spectra
    type: boolean?
  sirius__profile:
    doc: Name of the configuration profile
    type: string?
  sirius__formulas:
    doc: Specify the neutral molecular formula of the measured compound to compute its tree or a list of candidate formulas the method should discriminate. Omit this option if you want to consider all possible molecular formulas
    type: string?
  sirius__ions_enforced:
    doc: "The iontype/adduct of the MS/MS data. Example: [M+H]+, \n[M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a \ncomma separated list of adducts."
    type: string?
  sirius__candidates:
    doc: The number of formula candidates in the SIRIUS output
    type: long?
  sirius__candidates_per_ion:
    doc: Minimum number of candidates in the output for each ionization. Set to force output of results for each possible ionization, even if not part of highest ranked results.
    type: long?
  sirius__elements_considered:
    doc: Set the allowed elements for rare element detection. Write SBrClBSe to allow the elements S,Br,Cl,B and Se.
    type: string?
  sirius__elements_enforced:
    doc: "Enforce elements for molecular formula determination. Write CHNOPSCl to allow the elements C, H, N, O, P, S and Cl. Add numbers in brackets to restrict the minimal and maximal allowed occurrence of these elements: CHNOP[5]S[8]Cl[1-2]. When one number is given then it is interpreted as upper bound."
    type: string?
  sirius__no_isotope_score:
    doc: Disable isotope pattern score.
    type: boolean?
  sirius__no_isotope_filter:
    doc: Disable molecular formula filter. When filtering is enabled, molecular formulas are excluded if their theoretical isotope pattern does not match the theoretical one, even if their MS/MS pattern has high score.
    type: boolean?
  sirius__ions_considered:
    doc: "the iontype/adduct of the MS/MS data. Example: [M+H]+, [M-H]-, [M+Cl]-, [M+Na]+, [M]+. You can also provide a comma separated list of adducts."
    type: string?
  sirius__db:
    doc: "Search formulas in the Union of the given databases db-name1,db-name2,db-name3. If no database is given all possible molecular formulas will be respected (no database is used). Example: possible DBs: ALL,BIO,PUBCHEM,MESH,HMDB,KNAPSACK,CHEBI,PUBMED,KEGG,HSDB,MACONDA,METACYC,GNPS,ZINCBIO,UNDP,YMDB,PLANTCYC,NORMAN,ADDITIONAL,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBCHEMANNOTATIONFOOD,KEGGMINE,ECOCYCMINE,YMDBMINE"
    type: string?
  sirius__solver:
    doc: "For GUROBI and CPLEX environment variables need to be configured. \n(see SIRIUS manual: https://boecker-lab.github.io/docs.sirius.github.io/install/)."
    type: string?
  fingerid__db:
    doc: "Search structures in the Union of the given databases db-name1,db-name2,db-name3. If no database is given all possible molecular formulas will be respected (no database is used). Example: possible DBs: ALL,BIO,PUBCHEM,MESH,HMDB,KNAPSACK,CHEBI,PUBMED,KEGG,HSDB,MACONDA,METACYC,GNPS,ZINCBIO,UNDP,YMDB,PLANTCYC,NORMAN,ADDITIONAL,PUBCHEMANNOTATIONBIO,PUBCHEMANNOTATIONDRUG,PUBCHEMANNOTATIONSAFETYANDTOXIC,PUBCHEMANNOTATIONFOOD,KEGGMINE,ECOCYCMINE,YMDBMINE"
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: AssayGeneratorMetabo
doc: Assay library generation from DDA data (Metabolomics)
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - AssayGeneratorMetabo
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json