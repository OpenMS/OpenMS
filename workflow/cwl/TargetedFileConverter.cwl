inputs:
  in:
    doc: "Input file to convert.\n See http://www.openms.de/current_doxygen/html/UTILS_TargetedFileConverter.html for format of OpenSWATH transition TSV file or SpectraST MRM file."
    type: File
  in_type:
    doc: "input file type -- default: determined from file extension or content\n"
    type: string?
  out:
    doc: Output file
    type: string
  out_type:
    doc: "Output file type -- default: determined from file extension or content\nNote: not all conversion paths work or make sense."
    type: string?
  legacy_traml_id:
    doc: "PQP to TraML: Should legacy TraML IDs be used?"
    type: boolean?
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
  algorithm__retentionTimeInterpretation:
    doc: How to interpret the provided retention time (the retention time column can either be interpreted to be in iRT, minutes or seconds)
    type: string?
  algorithm__override_group_label_check:
    doc: Override an internal check that assures that all members of the same PeptideGroupLabel have the same PeptideSequence (this ensures that only different isotopic forms of the same peptide can be grouped together in the same label group). Only turn this off if you know what you are doing.
    type: boolean?
  algorithm__force_invalid_mods:
    doc: Force reading even if invalid modifications are encountered (OpenMS may not recognize the modification)
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: TargetedFileConverter
doc: Converts different transition files for targeted proteomics / metabolomics analysis.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - TargetedFileConverter
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json