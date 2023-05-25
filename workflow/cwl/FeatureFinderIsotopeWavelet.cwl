inputs:
  in:
    doc: input file
    type: File
  out:
    doc: output file
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
  algorithm__max_charge:
    doc: The maximal charge state to be considered.
    type: long?
  algorithm__intensity_threshold:
    doc: "The final threshold t' is build upon the formula: t' = av+t*sd, where t is the intensity_threshold, av the average intensity within the wavelet transformed signal and sd the standard deviation of the transform. If you set intensity_threshold=-1, t' will be zero.\nAs the 'optimal' value for this parameter is highly data dependent, we would recommend to start with -1, which will also extract features with very low signal-to-noise ratio. Subsequently, one might increase the threshold to find an optimized trade-off between false positives and true positives. Depending on the dynamic range of your spectra, suitable value ranges include: -1, [0:10], and if your data features even very high intensity values, t can also adopt values up to around 30. Please note that this parameter is not of an integer type, s.t. you can also use t:=0.1, e.g."
    type: double?
  algorithm__intensity_type:
    doc: Determines the intensity type returned for the identified features. 'ref' (default) returns the sum of the intensities of each isotopic peak within an isotope pattern. 'trans' refers to the intensity of the monoisotopic peak within the wavelet transform. 'corrected' refers also to the transformed intensity with an attempt to remove the effects of the convolution. While the latter ones might be preferable for qualitative analyses, 'ref' might be the best option to obtain quantitative results. Please note that intensity values might be spoiled (in particular for the option 'ref'), as soon as patterns overlap (see also the explanations given in the class documentation of FeatureFinderAlgorihtmIsotopeWavelet).
    type: string?
  algorithm__check_ppm:
    doc: Enables/disables a ppm test vs. the averagine model, i.e. potential peptide masses are checked for plausibility. In addition, a heuristic correcting potential mass shifts induced by the wavelet is applied.
    type: boolean?
  algorithm__hr_data:
    doc: Must be true in case of high-resolution data, i.e. for spectra featuring large m/z-gaps (present in FTICR and Orbitrap data, e.g.). Please check a single MS scan out of your recording, if you are unsure.
    type: boolean?
  algorithm__sweep_line__rt_votes_cutoff:
    doc: Defines the minimum number of subsequent scans where a pattern must occur to be considered as a feature.
    type: long?
  algorithm__sweep_line__rt_interleave:
    doc: Defines the maximum number of scans (w.r.t. rt_votes_cutoff) where an expected pattern is missing. There is usually no reason to change the default value.
    type: long?
outputs:
  out:
    type: File
    outputBinding:
      glob: $(inputs.out)
label: FeatureFinderIsotopeWavelet
doc: Detects two-dimensional features in LC-MS data.
cwlVersion: v1.2
class: CommandLineTool
baseCommand:
  - FeatureFinderIsotopeWavelet
requirements:
  InlineJavascriptRequirement: {}
  InitialWorkDirRequirement:
    listing:
      - entryname: cwl_inputs.json
        entry: $(JSON.stringify(inputs))
arguments:
  - -ini
  - cwl_inputs.json