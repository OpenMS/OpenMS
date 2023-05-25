inputs:
  in:
    doc: "input file in mzData format.\nNote: In mode 'mascot_out' a Mascot results file (.mascotXML) is read"
    type: File
  out:
    doc: "output file in idXML format.\nNote: In mode 'mascot_in' Mascot generic format is written."
    type: string
  out_type:
    doc: output file type (for TOPPAS)
    type: string?
  instrument:
    doc: the instrument that was used to measure the spectra
    type: string?
  precursor_mass_tolerance:
    doc: the precursor mass tolerance
    type: double?
  peak_mass_tolerance:
    doc: the peak mass tolerance
    type: double?
  taxonomy:
    doc: the taxonomy
    type: string?
  modifications:
    doc: the modifications i.e. Carboxymethyl (C)
    type: string[]?
  variable_modifications:
    doc: the variable modifications i.e. Carboxymethyl (C)
    type: string[]?
  charges:
    doc: the different charge states
    type: string[]?
  db:
    doc: the database to search in
    type: string?
  hits:
    doc: the number of hits to report
    type: string?
  cleavage:
    doc: The enzyme descriptor to the enzyme used for digestion. (Trypsin is default, None would be best for peptide input or unspecific digestion, for more please refer to your mascot server).
    type: string?
  missed_cleavages:
    doc: number of allowed missed cleavages
    type: long?
  sig_threshold:
    doc: significance threshold
    type: double?
  pep_homol:
    doc: peptide homology threshold
    type: double?
  pep_ident:
    doc: peptide ident threshold
    type: double?
  pep_rank:
    doc: peptide rank
    type: long?
  prot_score:
    doc: protein score
    type: double?
  pep_score:
    doc: peptide score
    type: double?
  pep_exp_z:
    doc: peptide expected charge
    type: long?
  show_unassigned:
    doc: show_unassigned
    type: long?
  first_dim_rt:
    doc: additional information which is added to every peptide identification as metavalue if set > 0
    type: double?
  boundary:
    doc: MIME boundary for mascot output format
    type: string?
  mass_type:
    doc: mass type
    type: string?
  mascot_directory:
    doc: the directory in which mascot is located
    type: string?
  temp_data_directory:
    doc: a directory in which some temporary files can be stored
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
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: MascotAdapter
doc: Annotates MS/MS spectra using Mascot.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MascotAdapter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json