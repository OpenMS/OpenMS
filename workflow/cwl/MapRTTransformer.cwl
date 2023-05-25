inputs:
  in:
    doc: Input file to transform (separated by blanks)
    type: File?
  out:
    doc: Output file (same file type as 'in'). This option or 'trafo_out' has to be provided; they can be used together.
    type: string
  trafo_in:
    doc: Transformation to apply
    type: File
  trafo_out:
    doc: Transformation output file. This option or 'out' has to be provided; they can be used together.
    type: string
  invert:
    doc: Invert transformation (approximatively) before applying it
    type: boolean?
  store_original_rt:
    doc: Store the original retention times (before transformation) as meta data in the output file
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
  model__type:
    doc: Type of model
    type: string?
  model__linear__symmetric_regression:
    doc: Perform linear regression on 'y - x' vs. 'y + x', instead of on 'y' vs. 'x'.
    type: boolean?
  model__linear__x_weight:
    doc: Weight x values
    type: string?
  model__linear__y_weight:
    doc: Weight y values
    type: string?
  model__linear__x_datum_min:
    doc: Minimum x value
    type: double?
  model__linear__x_datum_max:
    doc: Maximum x value
    type: double?
  model__linear__y_datum_min:
    doc: Minimum y value
    type: double?
  model__linear__y_datum_max:
    doc: Maximum y value
    type: double?
  model__b_spline__wavelength:
    doc: Determines the amount of smoothing by setting the number of nodes for the B-spline. The number is chosen so that the spline approximates a low-pass filter with this cutoff wavelength. The wavelength is given in the same units as the data; a higher value means more smoothing. '0' sets the number of nodes to twice the number of input points.
    type: double?
  model__b_spline__num_nodes:
    doc: Number of nodes for B-spline fitting. Overrides 'wavelength' if set (to two or greater). A lower value means more smoothing.
    type: long?
  model__b_spline__extrapolate:
    doc: "Method to use for extrapolation beyond the original data range. 'linear': Linear extrapolation using the slope of the B-spline at the corresponding endpoint. 'b_spline': Use the B-spline (as for interpolation). 'constant': Use the constant value of the B-spline at the corresponding endpoint. 'global_linear': Use a linear fit through the data (which will most probably introduce discontinuities at the ends of the data range)."
    type: string?
  model__b_spline__boundary_condition:
    doc: "Boundary condition at B-spline endpoints: 0 (value zero), 1 (first derivative zero) or 2 (second derivative zero)"
    type: long?
  model__lowess__span:
    doc: Fraction of datapoints (f) to use for each local regression (determines the amount of smoothing). Choosing this parameter in the range .2 to .8 usually results in a good fit.
    type: double?
  model__lowess__num_iterations:
    doc: Number of robustifying iterations for lowess fitting.
    type: long?
  model__lowess__delta:
    doc: Nonnegative parameter which may be used to save computations (recommended value is 0.01 of the range of the input, e.g. for data ranging from 1000 seconds to 2000 seconds, it could be set to 10). Setting a negative value will automatically do this.
    type: double?
  model__lowess__interpolation_type:
    doc: "Method to use for interpolation between datapoints computed by lowess. 'linear': Linear interpolation. 'cspline': Use the cubic spline for interpolation. 'akima': Use an akima spline for interpolation"
    type: string?
  model__lowess__extrapolation_type:
    doc: "Method to use for extrapolation outside the data range. 'two-point-linear': Uses a line through the first and last point to extrapolate. 'four-point-linear': Uses a line through the first and second point to extrapolate in front and and a line through the last and second-to-last point in the end. 'global-linear': Uses a linear regression to fit a line through all data points and use it for interpolation."
    type: string?
  model__interpolated__interpolation_type:
    doc: Type of interpolation to apply.
    type: string?
  model__interpolated__extrapolation_type:
    doc: "Type of extrapolation to apply: two-point-linear: use the first and last data point to build a single linear model, four-point-linear: build two linear models on both ends using the first two / last two points, global-linear: use all points to build a single linear model. Note that global-linear may not be continuous at the border."
    type: string?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
  trafo_out:
    type: File
    outputBinding:
      glob: $(inputs.trafo_out)
label: MapRTTransformer
doc: Applies retention time transformations to maps.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - MapRTTransformer
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json