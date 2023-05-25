inputs:
  in:
    doc: featureXML or consensusXML file
    type: File
  out:
    doc: mzTab file
    type: string
  out_annotation:
    doc: A copy of the input file, annotated with matching hits from the database.
    type: string
  positive_adducts:
    doc: This file contains the list of potential positive adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts. By default CHEMISTRY/PositiveAdducts.tsv in OpenMS/share is used.
    type: File
  negative_adducts:
    doc: This file contains the list of potential negative adducts that will be looked for in the database. Edit the list if you wish to exclude/include adducts. By default CHEMISTRY/NegativeAdducts.tsv in OpenMS/share is used.
    type: File
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
  db__mapping:
    doc: Database input file(s), containing three tab-separated columns of mass, formula, identifier. If 'mass' is 0, it is re-computed from the molecular sum formula. By default CHEMISTRY/HMDBMappingFile.tsv in OpenMS/share is used! If empty, the default will be used.
    type: File[]
  db__struct:
    doc: Database input file(s), containing four tab-separated columns of identifier, name, SMILES, INCHI.The identifier should match with mapping file. SMILES and INCHI are reported in the output, but not used otherwise. By default CHEMISTRY/HMDB2StructMapping.tsv in OpenMS/share is used! If empty, the default will be used.
    type: File[]
  algorithm__mass_error_value:
    doc: Tolerance allowed for accurate mass search.
    type: double?
  algorithm__mass_error_unit:
    doc: Unit of mass error (ppm or Da)
    type: string?
  algorithm__ionization_mode:
    doc: Positive or negative ionization mode? If 'auto' is used, the first feature of the input map must contain the meta-value 'scan_polarity'. If its missing, the tool will exit with error.
    type: string?
  algorithm__isotopic_similarity:
    doc: Computes a similarity score for each hit (only if the feature exhibits at least two isotopic mass traces).
    type: string?
  algorithm__use_feature_adducts:
    doc: Whether to filter AMS candidates mismatching available feature adduct annotation.
    type: string?
  algorithm__keep_unidentified_masses:
    doc: Keep features that did not yield any DB hit.
    type: string?
  algorithm__id_format:
    doc: Use legacy (ProteinID/PeptideID based storage of metabolomics data) with mzTab-v1.0.0 as output format or novel Identification Data (ID) with mzTab-v2.0.0-M as output format (ID and its MzTab-M output is currently only support for featureXML files).
    type: string?
  algorithm__mzTab__exportIsotopeIntensities:
    doc: "[featureXML input only] Export column with available isotope trace intensities (opt_global_MTint)"
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  out_annotation:
    type: File
    outputBinding:
      glob: $(inputs.out_annotation)
label: AccurateMassSearch
doc: Match MS signals to molecules from a database by mass.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - AccurateMassSearch
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json