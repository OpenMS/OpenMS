inputs:
  in:
    doc: Peptide multiplets with assigned sequence information
    type: File
  in_blacklist:
    doc: Optional input containing spectral peaks blacklisted during feature detection. Needed for generation of dummy features.
    type: File?
  out:
    doc: Complete peptide multiplets.
    type: string
  out_conflicts:
    doc: Optional output containing peptide multiplets without ID annotation or with conflicting quant/ID information.
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
    doc: "Labels used for labelling the samples. [...] specifies the labels for a single sample. For example\n\n[][Lys8,Arg10]        ... SILAC\n[][Lys4,Arg6][Lys8,Arg10]        ... triple-SILAC\n[Dimethyl0][Dimethyl6]        ... Dimethyl\n[Dimethyl0][Dimethyl4][Dimethyl8]        ... triple Dimethyl\n[ICPL0][ICPL4][ICPL6][ICPL10]        ... ICPL"
    type: string?
  algorithm__missed_cleavages:
    doc: Maximum number of missed cleavages due to incomplete digestion. (Only relevant if enzymatic cutting site coincides with labelling site. For example, Arg/Lys in the case of trypsin digestion and SILAC labelling.)
    type: long?
  algorithm__mass_tolerance:
    doc: Mass tolerance in Da for matching the mass shifts in the detected peptide multiplet to the theoretical mass shift pattern.
    type: double?
  algorithm__mz_tolerance:
    doc: m/z tolerance in ppm for checking if dummy feature vicinity was blacklisted.
    type: long?
  algorithm__rt_tolerance:
    doc: Retention time tolerance in seconds for checking if dummy feature vicinity was blacklisted.
    type: long?
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
  out_conflicts:
    type: File
    outputBinding:
      glob: $(inputs.out_conflicts)
label: MultiplexResolver
doc: Completes peptide multiplets and resolves conflicts within them.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MultiplexResolver
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json