inputs:
  in:
    doc: Input file
    type: File
  out:
    doc: Output file
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
  algorithm__rt_estimate:
    doc: If 'true' the optimal RT pair distance and deviation are estimated by fitting a gaussian distribution to the histogram of pair distance. Note that this works only datasets with a significant amount of pairs! If 'false' the parameters 'rt_pair_dist', 'rt_dev_low' and 'rt_dev_high' define the optimal distance.
    type: string?
  algorithm__rt_pair_dist:
    doc: optimal pair distance in RT [sec] from light to heavy feature
    type: double?
  algorithm__rt_dev_low:
    doc: maximum allowed deviation below optimal retention time distance
    type: double?
  algorithm__rt_dev_high:
    doc: maximum allowed deviation above optimal retention time distance
    type: double?
  algorithm__mz_pair_dists:
    doc: optimal pair distances in m/z [Th] for features with charge +1 (adapted to +2, +3, .. by division through charge)
    type: double[]?
  algorithm__mz_dev:
    doc: "maximum allowed deviation from optimal m/z distance\n"
    type: double?
  algorithm__mrm:
    doc: this option should be used if the features correspond mrm chromatograms (additionally the precursor is taken into account)
    type: boolean?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FeatureLinkerLabeled
doc: Groups corresponding isotope-labeled features in a feature map.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureLinkerLabeled
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json