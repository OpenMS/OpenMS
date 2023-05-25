inputs:
  in:
    doc: Identifications from searching a target-decoy database.
    type: File
  out:
    doc: Identifications with annotated FDR
    type: string
  PSM:
    doc: Perform FDR calculation on PSM level
    type: string?
  peptide:
    doc: "Perform FDR calculation on peptide level and annotates it as meta value\n(Note: if set, also calculates FDR/q-value on PSM level.)"
    type: boolean?
  protein:
    doc: Perform FDR calculation on protein level
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
  FDR__PSM:
    doc: Filter PSMs based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)
    type: double?
  FDR__protein:
    doc: Filter proteins based on q-value (e.g., 0.05 = 5% FDR, disabled for 1)
    type: double?
  FDR__cleanup__remove_proteins_without_psms:
    doc: Remove proteins without PSMs (due to being decoy or below PSM FDR threshold).
    type: string?
  FDR__cleanup__remove_psms_without_proteins:
    doc: Remove PSMs without proteins (due to being decoy or below protein FDR threshold).
    type: string?
  FDR__cleanup__remove_spectra_without_psms:
    doc: "Remove spectra without PSMs (due to being decoy or below protein FDR threshold). Caution: if remove_psms_without_proteins is false, protein level filtering does not propagate."
    type: string?
  algorithm__no_qvalues:
    doc: If 'true' strict FDRs will be calculated instead of q-values (the default)
    type: boolean?
  algorithm__use_all_hits:
    doc: If 'true' not only the first hit, but all are used (peptides only)
    type: boolean?
  algorithm__split_charge_variants:
    doc: If 'true' charge variants are treated separately (for peptides of combined target/decoy searches only).
    type: boolean?
  algorithm__treat_runs_separately:
    doc: If 'true' different search runs are treated separately (for peptides of combined target/decoy searches only).
    type: boolean?
  algorithm__add_decoy_peptides:
    doc: If 'true' decoy peptides will be written to output file, too. The q-value is set to the closest target score.
    type: boolean?
  algorithm__add_decoy_proteins:
    doc: If 'true' decoy proteins will be written to output file, too. The q-value is set to the closest target score.
    type: boolean?
  algorithm__conservative:
    doc: If 'true' (D+1)/T instead of (D+1)/(T+D) is used as a formula.
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FalseDiscoveryRate
doc: Estimates the false discovery rate on peptide and protein level using decoy searches.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FalseDiscoveryRate
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json