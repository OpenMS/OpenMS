inputs:
  in:
    doc: Input files to align (all must have the same file type)
    type: File[]
  copy_data:
    doc: Copy data (faster, more memory required) or reload data (slower, less memory required) when aligning many files.
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
  algorithm__model_type:
    doc: Options to control the modeling of retention time transformations from data
    type: string?
  algorithm__model__type:
    doc: Type of model
    type: string?
  algorithm__model__linear__symmetric_regression:
    doc: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.
    type: boolean?
  algorithm__model__linear__x_weight:
    doc: Weight x values
    type: string?
  algorithm__model__linear__y_weight:
    doc: Weight y values
    type: string?
  algorithm__model__linear__x_datum_min:
    doc: Minimum x value
    type: double?
  algorithm__model__linear__x_datum_max:
    doc: Maximum x value
    type: double?
  algorithm__model__linear__y_datum_min:
    doc: Minimum y value
    type: double?
  algorithm__model__linear__y_datum_max:
    doc: Maximum y value
    type: double?
  algorithm__model__b_spline__wavelength:
    doc: Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.
    type: double?
  algorithm__model__b_spline__num_nodes:
    doc: Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.
    type: long?
  algorithm__model__b_spline__extrapolate:
    doc: "Method to use for extrapolation beyond the original data range. 'linear': Linear extrapolation using the slope of the B-spline at the corresponding endpoint. 'b_spline': Use the B-spline (as for interpolation). 'constant': Use the constant value of the B-spline at the corresponding endpoint. 'global_linear': Use a linear fit through the data (which will most probably introduce discontinuities at the ends of the data range)."
    type: string?
  algorithm__model__b_spline__boundary_condition:
    doc: "Boundary condition at B-spline endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)"
    type: long?
  algorithm__model__lowess__span:
    doc: Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.
    type: double?
  algorithm__model__lowess__num_iterations:
    doc: Number of robustifying iterations for lowess fitting.
    type: long?
  algorithm__model__lowess__delta:
    doc: Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.
    type: double?
  algorithm__model__lowess__interpolation_type:
    doc: "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation"
    type: string?
  algorithm__model__lowess__extrapolation_type:
    doc: "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation."
    type: string?
  algorithm__model__interpolated__interpolation_type:
    doc: Type of interpolation to apply.
    type: string?
  algorithm__model__interpolated__extrapolation_type:
    doc: "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border."
    type: string?
  algorithm__align_algorithm__score_type:
    doc: Name of the score type to use for ranking and filtering (.oms input only). If left empty, a score type is picked automatically.
    type: string?
  algorithm__align_algorithm__score_cutoff:
    doc: Use only IDs above a score cut-off (parameter 'min_score') for alignment?
    type: boolean?
  algorithm__align_algorithm__min_score:
    doc: "If 'score_cutoff' is 'true': Minimum score for an ID to be considered.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides."
    type: double?
  algorithm__align_algorithm__min_run_occur:
    doc: "Minimum number of runs (incl. reference, if any) in which a peptide must occur to be used for the alignment.\nUnless you have very few runs or identifications, increase this value to focus on more informative peptides."
    type: long?
  algorithm__align_algorithm__max_rt_shift:
    doc: "Maximum realistic RT difference for a peptide (median per run vs. reference). Peptides with higher shifts (outliers) are not used to compute the alignment.\nIf 0, no limit (disable filter); if > 1, the final value in seconds; if <= 1, taken as a fraction of the range of the reference RT scale."
    type: double?
  algorithm__align_algorithm__use_unassigned_peptides:
    doc: Should unassigned peptide identifications be used when computing an alignment of feature or consensus maps? If 'false', only peptide IDs assigned to features will be used.
    type: string?
  algorithm__align_algorithm__use_feature_rt:
    doc: "When aligning feature or consensus maps, don't use the retention time of a peptide identification directly; instead, use the retention time of the centroid of the feature (apex of the elution profile) that the peptide was matched to. If different identifications are matched to one feature, only the peptide closest to the centroid in RT is used.\nPrecludes 'use_unassigned_peptides'."
    type: string?
  algorithm__align_algorithm__use_adducts:
    doc: If IDs contain adducts, treat differently adducted variants of the same molecule as different.
    type: string?
outputs:
  {}
label: MapAlignerTreeGuided
doc: Tree guided correction of retention time distortions between maps.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MapAlignerTreeGuided
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json