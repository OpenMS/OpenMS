inputs:
  in:
    doc: Input peak file
    type: File
  out:
    doc: "Output file "
    type: string
  rscript_executable:
    doc: "Path to the Rscript executable (default: 'Rscript')."
    type: File?
  ppm_match_tolerance:
    doc: Finding calibrants in raw data uses this tolerance (for lock masses and ID's).
    type: double?
  ms_level:
    doc: Target MS levels to apply the transformation onto. Does not affect calibrant collection.
    type: long[]?
  RT_chunking:
    doc: RT window (one-sided, i.e. left->center, or center->right) around an MS scan in which calibrants are collected to build a model. Set to -1 to use ALL calibrants for all scans, i.e. a global model.
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
  cal__id_in:
    doc: Identifications or features whose peptide ID's serve as calibration masses.
    type: File?
  cal__lock_in:
    doc: "Input file containing reference m/z values (text file with each line as: m/z ms-level charge) which occur in all scans."
    type: File?
  cal__lock_out:
    doc: Optional output file containing peaks from 'in' which were matched to reference m/z values. Useful to see which peaks were used for calibration.
    type: string
  cal__lock_fail_out:
    doc: Optional output file containing lock masses which were NOT found or accepted(!) in data from 'in'. Useful to see which peaks were used for calibration.
    type: string
  cal__lock_require_mono:
    doc: Require all lock masses to be monoisotopic, i.e. not the iso1, iso2 etc ('charge' column is used to determine the spacing). Peaks which are not mono-isotopic are not used.
    type: boolean?
  cal__lock_require_iso:
    doc: Require all lock masses to have at least the +1 isotope. Peaks without isotope pattern are not used.
    type: boolean?
  cal__model_type:
    doc: Type of function to be fitted to the calibration points.
    type: string?
  RANSAC__enabled:
    doc: Apply RANSAC to calibration points to remove outliers before fitting a model.
    type: boolean?
  RANSAC__threshold:
    doc: Threshold for accepting inliers (instrument precision (not accuracy!) as ppm^2 distance)
    type: double?
  RANSAC__pc_inliers:
    doc: Minimum percentage (of available data) of inliers (<threshold away from model) to accept the model.
    type: long?
  RANSAC__iter:
    doc: "Maximal # iterations."
    type: long?
  goodness__median:
    doc: The median ppm error of calibrated masses must be smaller than this threshold.
    type: double?
  goodness__MAD:
    doc: The median absolute deviation of the ppm error of calibrated masses must be smaller than this threshold.
    type: double?
  quality_control__models:
    doc: Table of model parameters for each spectrum.
    type: string
  quality_control__models_plot:
    doc: Plot image of model parameters for each spectrum.
    type: string
  quality_control__residuals:
    doc: Table of pre- and post calibration errors.
    type: string
  quality_control__residuals_plot:
    doc: Plot image of pre- and post calibration errors.
    type: string
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  cal__lock_out:
    type: File
    outputBinding:
      glob: $(inputs.cal__lock_out)
  cal__lock_fail_out:
    type: File
    outputBinding:
      glob: $(inputs.cal__lock_fail_out)
  quality_control__models:
    type: File
    outputBinding:
      glob: $(inputs.quality_control__models)
  quality_control__models_plot:
    type: File
    outputBinding:
      glob: $(inputs.quality_control__models_plot)
  quality_control__residuals:
    type: File
    outputBinding:
      glob: $(inputs.quality_control__residuals)
  quality_control__residuals_plot:
    type: File
    outputBinding:
      glob: $(inputs.quality_control__residuals_plot)
label: InternalCalibration
doc: Applies an internal mass recalibration.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - InternalCalibration
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json