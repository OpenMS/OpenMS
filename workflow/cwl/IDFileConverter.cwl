inputs:
  in:
    doc: "Input file or directory containing the data to convert. This may be:\n- a single file in OpenMS database format (.oms),\n- a single file in a multi-purpose XML format (.idXML, .mzid, .pepXML, .protXML),\n- a single file in a search engine-specific format (Mascot: .mascotXML, OMSSA: .omssaXML, X! Tandem: .xml, Percolator: .psms, xQuest: .xquest.xml),\n- a single file in fasta format (can only be used to generate a theoretical mzML),\n- a single text file (tab separated) with one line for all peptide sequences matching a spectrum (top N hits),\n- for Sequest results, a directory containing .out files.\n"
    type: File
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type (default: determined from file extension)"
    type: string?
  mz_file:
    doc: "[pepXML, Sequest, Mascot, X! Tandem, mzid, Percolator only] Retention times and native spectrum ids (spectrum_references) will be looked up in this file"
    type: File?
  mz_name:
    doc: "[pepXML only] Experiment filename/path (extension will be removed) to match in the pepXML file ('base_name' attribute). Only necessary if different from 'mz_file'."
    type: string?
  peptideprophet_analyzed:
    doc: "[pepXML output only] Write output in the format of a PeptideProphet analysis result. By default a 'raw' pepXML is produced that contains only search engine results."
    type: boolean?
  score_type:
    doc: "[Percolator only] Which of the Percolator scores to report as 'the' score for a peptide hit"
    type: string?
  ignore_proteins_per_peptide:
    doc: "[Sequest only] Workaround to deal with .out files that contain e.g. \"+1\" in references column,\nbut do not list extra references in subsequent lines (try -debug 3 or 4)"
    type: boolean?
  scan_regex:
    doc: "[Mascot, pepXML, Percolator only] Regular expression used to extract the scan number or retention time. See documentation for details."
    type: string?
  no_spectra_data_override:
    doc: "[+mz_file only] Avoid overriding 'spectra_data' in protein identifications if 'mz_file' is given and 'spectrum_reference's are added/updated. Use only if you are sure it is absolutely the same 'mz_file' as used for identification."
    type: boolean?
  no_spectra_references_override:
    doc: "[+mz_file only] Avoid overriding 'spectrum_reference' in peptide identifications if 'mz_file' is given and a 'spectrum_reference' is already present."
    type: boolean?
  add_ionmatch_annotation:
    doc: "[+mz_file only] Annotate the identifications with ion matches from spectra in 'mz_file' using the given tolerance (in Da). This will take quite some time."
    type: double?
  concatenate_peptides:
    doc: "[FASTA output only] Will concatenate the top peptide hits to one peptide sequence, rather than write a new peptide for each hit."
    type: boolean?
  number_of_hits:
    doc: "[FASTA output only] Controls how many peptide hits will be exported. A value of 0 or less exports all hits."
    type: long?
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
  fasta_to_mzml__isotope_model:
    doc: Model to use for isotopic peaks ('none' means no isotopic peaks are added, 'coarse' adds isotopic peaks in unit mass distance, 'fine' uses the hyperfine isotopic generator to add accurate isotopic peaks. Note that adding isotopic peaks is very slow.
    type: string?
  fasta_to_mzml__max_isotope:
    doc: Defines the maximal isotopic peak which is added if 'isotope_model' is 'coarse'
    type: long?
  fasta_to_mzml__max_isotope_probability:
    doc: Defines the maximal isotopic probability to cover if 'isotope_model' is 'fine'
    type: double?
  fasta_to_mzml__add_metainfo:
    doc: Adds the type of peaks as metainfo to the peaks, like y8+, [M-H2O+2H]++
    type: boolean?
  fasta_to_mzml__add_losses:
    doc: Adds common losses to those ion expect to have them, only water and ammonia loss is considered
    type: boolean?
  fasta_to_mzml__sort_by_position:
    doc: Sort output by position
    type: string?
  fasta_to_mzml__add_precursor_peaks:
    doc: Adds peaks of the unfragmented precursor ion to the spectrum
    type: boolean?
  fasta_to_mzml__add_all_precursor_charges:
    doc: Adds precursor peaks with all charges in the given range
    type: boolean?
  fasta_to_mzml__add_abundant_immonium_ions:
    doc: Add most abundant immonium ions (for Proline, Cystein, Iso/Leucine, Histidin, Phenylalanin, Tyrosine, Tryptophan)
    type: boolean?
  fasta_to_mzml__add_first_prefix_ion:
    doc: If set to true e.g. b1 ions are added
    type: boolean?
  fasta_to_mzml__add_y_ions:
    doc: Add peaks of y-ions to the spectrum
    type: string?
  fasta_to_mzml__add_b_ions:
    doc: Add peaks of b-ions to the spectrum
    type: string?
  fasta_to_mzml__add_a_ions:
    doc: Add peaks of a-ions to the spectrum
    type: boolean?
  fasta_to_mzml__add_c_ions:
    doc: Add peaks of c-ions to the spectrum
    type: boolean?
  fasta_to_mzml__add_x_ions:
    doc: Add peaks of  x-ions to the spectrum
    type: boolean?
  fasta_to_mzml__add_z_ions:
    doc: Add peaks of z-ions to the spectrum
    type: boolean?
  fasta_to_mzml__y_intensity:
    doc: Intensity of the y-ions
    type: double?
  fasta_to_mzml__b_intensity:
    doc: Intensity of the b-ions
    type: double?
  fasta_to_mzml__a_intensity:
    doc: Intensity of the a-ions
    type: double?
  fasta_to_mzml__c_intensity:
    doc: Intensity of the c-ions
    type: double?
  fasta_to_mzml__x_intensity:
    doc: Intensity of the x-ions
    type: double?
  fasta_to_mzml__z_intensity:
    doc: Intensity of the z-ions
    type: double?
  fasta_to_mzml__relative_loss_intensity:
    doc: Intensity of loss ions, in relation to the intact ion intensity
    type: double?
  fasta_to_mzml__precursor_intensity:
    doc: Intensity of the precursor peak
    type: double?
  fasta_to_mzml__precursor_H2O_intensity:
    doc: Intensity of the H2O loss peak of the precursor
    type: double?
  fasta_to_mzml__precursor_NH3_intensity:
    doc: Intensity of the NH3 loss peak of the precursor
    type: double?
  fasta_to_mzml__enzyme:
    doc: Enzym used to digest the fasta proteins
    type: string?
  fasta_to_mzml__missed_cleavages:
    doc: Number of allowed missed cleavages while digesting the fasta proteins
    type: long?
  fasta_to_mzml__min_charge:
    doc: Minimum charge
    type: long?
  fasta_to_mzml__max_charge:
    doc: Maximum charge
    type: long?
  fasta_to_mzml__precursor_charge:
    doc: "Manually set precursor charge. (default: 0, meaning max_charge + 1 will be used as precursor charge)"
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: IDFileConverter
doc: Converts identification engine file formats.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - IDFileConverter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json