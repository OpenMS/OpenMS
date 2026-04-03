// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// Consumers
#include <OpenMS/FORMAT/DATAACCESS/MSDataWritingConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataSqlConsumer.h>
#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>

// Files
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/SwathFile.h>
#include <OpenMS/FORMAT/DATAACCESS/MSDataTransformingConsumer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathWindowLoader.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathQC.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionTSVFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/config.h>
#include <filesystem>
#include <OpenMS/SYSTEM/File.h>

// Kernel and implementations
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMS.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessTransforming.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SpectrumAccessOpenMSInMemory.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/SwathMap.h>

// Helpers
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

// Algorithms
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>
#include <OpenMS/ANALYSIS/OPENSWATH/ChromatogramExtractor.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/FORMAT/TransformationXMLFile.h>

#include <algorithm>
#include <cassert>
#include <limits>
#include <memory>

// #define OPENSWATH_WORKFLOW_DEBUG

using namespace OpenMS;

// OpenMS base classes
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/OpenSwathBase.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>


#include <unordered_map>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page TOPP_OpenSwathWorkflow OpenSwathWorkflow

@brief Complete workflow to run OpenSWATH

This implements the OpenSWATH workflow as described in Rost and Rosenberger
et al. (Nature Biotechnology, 2014) and provides a complete, integrated
analysis tool without the need to run multiple tools consecutively. See also
http://openswath.org/ for additional documentation.

It executes the following steps in order, which is implemented in @ref OpenMS::OpenSwathWorkflow "OpenSwathWorkflow":

<ul>
  <li>Reading of input SRM/MRM/PRM/DIA(PASEF) mzML files or Bruker .d (TDF) directories (requires WITH_OPENTIMS)</li>
  <li>Computing the retention time transformation, mass-to-charge and ion mobility correction using calibrant peptides</li>
  <li>Reading of the transition list</li>
  <li>Extracting the specified transitions</li>
  <li>Scoring the peak groups in the extracted ion chromatograms (XIC)</li>
  <li>Reporting the peak groups and the chromatograms</li>
</ul>


See below or have a look at the INI file (via "OpenSwathWorkflow -write_ini myini.ini") for available parameters and more functionality.

<h3>Input</h3>

<h4>DIA (including diaPASEF)</h4>
DIA data (commonly referred to as SWATH) can be provided as full-scan mzML
files produced by the instrument. In single-file mode the tool expects the
acquisition to contain consistent cycles of MS1 and MS2 spectra (for example:
MS1, MS2 [400-425], MS2 [425-450], MS1, MS2 [400-425], MS2 [425-450]). Each
MS2 window is identified by its isolation lower/upper bounds. OpenSwathWorkflow
will attempt to read the SWATH windows from the file; if this is not possible
please supply a tab-separated SWATH windows file via the
`-swath_windows_file` parameter. Note: extraction windows must not overlap —
overlapping windows lead to peptides being extracted from more than one
window and produce ambiguous assignments.

Split-file input is supported via the `-split_file_input` flag: in this mode
each SWATH window (and optionally the MS1 map) may be provided as a separate
mzML file (n+1 files). This is useful for very large datasets or when the
instrument exports individual window files.

Since files can be large, it is recommended to avoid loading whole datasets
into memory; use the `-readOptions` (for example `cacheWorkingInMemory`) to
cache or reduce data in advance.

<h4>Bruker .d (diaPASEF, requires WITH_OPENTIMS)</h4>
Bruker TimsTOF .d directories containing DIA-PASEF data can be passed directly
without prior mzML conversion when OpenMS is built with the WITH_OPENTIMS option.
The tool automatically discovers SWATH windows from the TDF metadata and
partitions spectra accordingly. Ion mobility lower/upper limits are attached as
spectrum meta values so that PASEF windows sharing the same m/z range but
differing in ion mobility are correctly distinguished.

<h4>PRM</h4>
PRM (parallel reaction monitoring) is a targeted MS2 acquisition mode. PRM
data typically contains MS2 traces targeted to specific precursors and may
have windows that overlap or differ in ordering from classic SWATH runs.
When analyzing PRM-like data consider using the `-matching_window_only` option
to restrict extraction to the best-matching window for each assay. PRM
datasets can be provided as single-file mzMLs or as split files as described
above.

<h4>SRM / MRM (targeted chromatogram input)</h4>
SRM/MRM are supported as targeted chromatogram input. These inputs often consist of chromatogram-only mzMLfiles (exported XICs) rather than full spectra. The targeted data loader in OpenMS will read chromatogram-only mzMLs and attempt to map chromatograms to assays using the supplied transition list. When providing pre-extracted chromatograms, note that MS1 extraction (`-enable_ms1`) does not apply — the tool will operate on the supplied chromatograms and use the mapping logic to associate traces with transitions.


The assay library (transition list) is provided through the @p -tr parameter and can be in one of the following formats:

  <ul>
    <li> @ref OpenMS::TransitionPQPFile "OpenSWATH PQP SQLite files" (Recommended) </li>
    <li> @ref OpenMS::TransitionTSVFile "OpenSWATH TSV transition lists" </li>
    <li> @ref OpenMS::TransitionParquetFile "OpenSWATH Parquet library (.oswpq)" </li>
    <li> @ref OpenMS::TraMLFile "TraML" </li>
    <li> SpectraST MRM transition lists </li>
    <li> Skyline transition lists </li>
    <li> Spectronaut transition lists </li>
  </ul>

<h3>Parameters</h3>
The current parameters are optimized for 2 hour gradients on SCIEX 5600 /
6600 TripleTOF instruments with a peak width of around 30 seconds using iRT
peptides.  If your chromatography differs, please consider adjusting
@p -Scoring:TransitionGroupPicker:min_peak_width  to allow for smaller or larger
peaks and adjust the @p -rt_extraction_window to use a different extraction
window for the retention time. In m/z domain, consider adjusting
@p -mz_extraction_window to your instrument resolution, which can be in Th or
ppm.

Furthermore, if you wish to use MS1 information, use the @p -enable_ms1 flag
and provide an MS1 map in addition to the SWATH data.

If you encounter issues with peak picking, try to disable peak filtering by
setting @p -Scoring:TransitionGroupPicker:compute_peak_quality false which will
disable the filtering of peaks by chromatographic quality. Furthermore, you
can adjust the smoothing parameters for the peak picking, by adjusting
@p -Scoring:TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length or using a
Gaussian smoothing based on your estimated peak width. Adjusting the signal
to noise threshold will make the peaks wider or smaller.

<h3>Output: Feature list, chromatograms, and ion mobilograms </h3>
The output of the OpenSwathWorkflow is a feature list, either as FeatureXML,
a @ref OpenMS::OSWFile "OpenSWATH SQLite file", or an OpenSWATH Parquet output
(use @p -out_features) while the SQLite output is more memory
friendly and can be directly used as input to other tools such as pyProphet (a Python
re-implementation of mProphet) software tool, see Reiter et al (2011, Nature
Methods).
If you analyze large datasets, it is recommended to only use the @ref OpenMS::OSWFile "OSWFile format".
For downstream analysis (e.g. using pyProphet) the @ref OpenMS::OSWFile "OSWFile format" is recommended.

In addition, the extracted chromatograms can be written out using the
@p -out_chrom parameter.

When processing ion mobility (diaPASEF) data, the extracted ion mobilograms
(XIMs) can optionally be saved to a Parquet file using the @p -out_mobilogram
parameter. The output file must have the @p .xim extension. The resulting file can be
read back using the @ref OpenMS::XIMParquetFile "XIMParquetFile" class.

<h4> Feature list output format </h4>

For more information on the feature tables in the @ref OpenMS::OSWFile "OpenSWATH SQLite file output", see @ref OpenMS::OpenSwathOSWWriter "the OpenSwathOSWWriter class".

<h3>Execution flow:</h3>

The overall execution flow for this tool is implemented in @ref OpenMS::OpenSwathWorkflow "OpenSwathWorkflow".

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathWorkflow.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathWorkflow.html

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES
class TOPPOpenSwathWorkflow
  : public TOPPOpenSwathBase
{
public:

  TOPPOpenSwathWorkflow()
    : TOPPOpenSwathBase("OpenSwathWorkflow", "Complete workflow to run OpenSWATH", true,
                        {
                          {"Roest, H.L. et al.",
                           "OpenSWATH enables automated, targeted analysis of data-independent acquisition MS data",
                           "Nature Biotechnology volume 32, pages 219–223 (2014)",
                           "https://doi.org/10.1038/nbt.2841"},
                          {"Rosenberger, G. et al.",
                           "Inference and quantification of peptidoforms in large sample cohorts by SWATH-MS",
                           "Nature Biotechnology volume 35, pages 781–788 (2017)",
                           "https://doi.org/10.1038/nbt.3908"},
                          {"Meier, F. et al.",
                           "diaPASEF: parallel accumulation–serial fragmentation combined with data-independent acquisition",
                           "Nature Methods volume 17, pages 1229–1236 (2020)",
                           "https://doi.org/10.1038/s41592-020-00998-0"}
                        })
  {
  }

protected:

  void registerOptionsAndFlags_() override
  {
    registerInputFileList_("in", "<files>", StringList(), "Input files separated by blank");
    StringList in_formats = {"mzML", "mzXML", "sqMass"};
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
    setValidFormats_("in", in_formats);

    registerInputFile_("tr", "<file>", "", "transition file ('TraML','tsv','pqp','oswpq')");
    StringList tr_formats = {"traML", "tsv", "pqp"};
    tr_formats.push_back("oswpq");
    setValidFormats_("tr", tr_formats);
    registerStringOption_("tr_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    setValidStrings_("tr_type", tr_formats);

    registerInputFile_("swath_windows_file", "<file>", "", "Optional, tab-separated file containing the SWATH windows for extraction: lower_offset upper_offset. Note that the first line is a header and will be skipped.", false);
    registerFlag_("sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from swath_windows_file", true);

    registerStringOption_("enable_ms1", "<true|false>", "true", "Extract the precursor ion trace(s) and use for scoring if present", false, true);
    setValidStrings_("enable_ms1", ListUtils::create<String>("true,false"));

    registerStringOption_("enable_ipf", "<true|false>", "true", "Enable additional scoring of identification assays using IPF (see online documentation)", false, true);
    setValidStrings_("enable_ipf", ListUtils::create<String>("true,false"));

    registerOutputFile_("out_features", "<file>", "", "feature output file, either .osw (PyProphet-compatible SQLite file), .oswpq, or .featureXML", false);
    std::vector<String> out_feature_formats = {"osw", "featureXML"};
    out_feature_formats.push_back("oswpq");
    setValidFormats_("out_features", out_feature_formats);

    registerStringOption_("out_features_type", "<type>", "", "input file type -- default: determined from file extension or content\n", false);
    setValidStrings_("out_features_type", out_feature_formats);

    registerOutputFile_("out_chrom", "<file>", "", "Also output all computed chromatograms output in mzML (chrom.mzML), sqMass (SQLite format) or xic (Parquet)", false, true);
    setValidFormats_("out_chrom", ListUtils::create<String>("mzML,sqMass,xic"));

    registerOutputFile_("out_mobilogram", "<file>", "", "Also output extracted ion mobilograms in Parquet format", false, true);
    setValidFormats_("out_mobilogram", ListUtils::create<String>("xim"));

    // additional QC data
    registerOutputFile_("out_qc", "<file>", "", "Optional QC meta data (charge distribution in MS1). Only works with mzML input files.", false, true);
    setValidFormats_("out_qc", ListUtils::create<String>("json"));


    // misc options
    registerDoubleOption_("min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the upper edge of a Swath window to still consider a precursor, in Thomson", false, true);
    registerFlag_("pasef", "data is PASEF data");

    registerDoubleOption_("rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract over the whole range, a value of 600 means to extract around +/- 300 s of the expected elution).", false);
    registerDoubleOption_("extra_rt_extraction_window", "<double>", 0.0, "Output an XIC with a RT-window by this much larger (e.g. to visually inspect a larger area of the chromatogram)", false, true);
    setMinFloat_("extra_rt_extraction_window", 0.0);
    registerDoubleOption_("ion_mobility_window", "<double>", -1, "Extraction window in ion mobility dimension (in 1/K0, milliseconds, or CCS depending on library). This is the full window size, e.g. a value of 10 milliseconds would extract 5 milliseconds on either side. -1 means extract over the whole range or ion mobility is not present. (Default for diaPASEF data: 0.06 1/K0, or ~10 CCS)", false);
    registerDoubleOption_("mz_extraction_window", "<double>", 50, "Extraction window in Thomson or ppm (see mz_extraction_window_unit)", false);
    setMinFloat_("mz_extraction_window", 0.0);
    registerStringOption_("mz_extraction_window_unit", "<name>", "ppm", "Unit for mz extraction", false, true);
    setValidStrings_("mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));

    // MS1 mz windows and ion mobility
    registerDoubleOption_("mz_extraction_window_ms1", "<double>", 50, "Extraction window used in MS1 in Thomson or ppm (see mz_extraction_window_ms1_unit)", false);
    setMinFloat_("mz_extraction_window_ms1", 0.0);
    registerStringOption_("mz_extraction_window_ms1_unit", "<name>", "ppm", "Unit of the MS1 m/z extraction window", false, true);
    setValidStrings_("mz_extraction_window_ms1_unit", ListUtils::create<String>("ppm,Th"));
    registerDoubleOption_("im_extraction_window_ms1", "<double>", -1, "Extraction window in ion mobility dimension for MS1 (in 1/K0, milliseconds, or CCS depending on library). -1 means this is not ion mobility data.", false);

    registerStringOption_("use_ms1_ion_mobility", "<name>", "true", "Also perform precursor extraction using the same ion mobility window as for fragment ion extraction", false, true);
    setValidStrings_("use_ms1_ion_mobility", ListUtils::create<String>("true,false"));

    registerStringOption_("matching_window_only", "<name>", "false", "Assume the input data is targeted / PRM-like data with potentially overlapping DIA windows. Will only attempt to extract each assay from the *best* matching DIA window (instead of all matching windows).", false, true);
    setValidStrings_("matching_window_only", ListUtils::create<String>("true,false"));

    // iRT mz and IM windows
    registerDoubleOption_("irt_mz_extraction_window", "<double>", 50, "Extraction window used for iRT and m/z correction in Thomson or ppm (see irt_mz_extraction_window_unit)", false, true);
    setMinFloat_("irt_mz_extraction_window", 0.0);
    registerStringOption_("irt_mz_extraction_window_unit", "<name>", "ppm", "Unit for mz extraction", false, true);
    setValidStrings_("irt_mz_extraction_window_unit", ListUtils::create<String>("Th,ppm"));
    registerDoubleOption_("irt_im_extraction_window", "<double>", -1, "Ion mobility extraction window used for iRT (in 1/K0 or milliseconds depending on library). -1 means do not perform ion mobility calibration", false, true);

    registerFlag_("split_file_input", "The input files each contain one single SWATH (alternatively: all SWATH are in separate files)", true);
    registerFlag_("use_elution_model_score", "Turn on elution model score (EMG fit to peak)", true);

    registerFlag_("append_oswpq", "If out_features is an oswpq archive, optionally append to the existing .oswpq archive instead of overwriting. This may be useful if you run separate instances of OpenSwathWorkflow for separate input files. (default: overwrite)", true);

    registerStringOption_("readOptions", "<name>", "normal", "Whether to run OpenSWATH directly on the input data, cache data to disk first or to perform a datareduction step first. If you choose cache, make sure to also set tempDirectory", false, true);
    setValidStrings_("readOptions", ListUtils::create<String>("normal,cache,cacheWorkingInMemory,workingInMemory"));

    registerStringOption_("tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory to store cached files for example", false, true);
    registerFlag_("keep_cached_files", "If set, do not remove cached files created in tempDirectory (disable automated cleanup)", false);

    registerStringOption_("extraction_function", "<name>", "tophat", "Function used to extract the signal", false, true);
    setValidStrings_("extraction_function", ListUtils::create<String>("tophat,bartlett"));

    registerIntOption_("batchSize", "<number>", 1000, "The batch size of chromatograms to process (0 means to only have one batch, sensible values are around 250-1000)", false, true);
    setMinInt_("batchSize", 0);
    registerIntOption_("outer_loop_threads", "<number>", -1, "How many threads should be used for the outer loop (-1 use all threads, use 4 to analyze 4 SWATH windows in memory at once).", false, true);

    registerIntOption_("ms1_isotopes", "<number>", 3, "The number of MS1 isotopes used for extraction", false, true);
    setMinInt_("ms1_isotopes", 0);

    registerSubsection_("Scoring", "Scoring parameters section");
    registerSubsection_("Library", "Library parameters section");

    registerSubsection_("Calibration", "Parameters for calibrant iRT peptides for RT normalization and mass / ion mobility correction.");
    registerSubsection_("Calibration:RTNormalization", "Parameters for the RTNormalization for iRT peptides. This specifies how the RT alignment is performed and how outlier detection is applied. Outlier detection can be done iteratively (by default) which removes one outlier per iteration or using the RANSAC algorithm.");
    registerSubsection_("Calibration:MassIMCorrection", "Parameters for the m/z and ion mobility calibration.");
  registerSubsection_("MRMMapping", "Parameters for mapping chromatograms to transitions (SRM/MRM data).");

    registerTOPPSubsection_("Debugging", "Debugging");
    registerOutputFile_("Debugging:irt_mzml", "<file>", "", "Chromatogram mzML containing the iRT peptides", false);
    setValidFormats_("Debugging:irt_mzml", ListUtils::create<String>("mzML"));
    registerOutputFile_("Debugging:irt_trafo", "<file>", "", "Transformation file for RT transform", false);
    setValidFormats_("Debugging:irt_trafo", ListUtils::create<String>("trafoXML"));
    registerStringList_("Debugging:disable_features", "<list>", StringList(),
      "Selectively disable features for debugging/benchmarking. "
      "Valid values: "
      "'no_IM_calibration' (skip ion mobility calibration by iRT peptides), "
      "'no_IM_windowing' (extract chromatograms over the full IM range and disable PASEF window matching; "
      "does not affect IM calibration -- use no_IM_calibration for that). "
      "Note: IM scoring is controlled separately via -Scoring:Scores:use_ion_mobility_scores (auto/true/false).", false, true);
    setValidStrings_("Debugging:disable_features",
      ListUtils::create<String>("no_IM_calibration,no_IM_windowing"));
  }

  Param getSubsectionDefaults_(const String& name) const override
  {
    if (name == "Scoring")
    {
      // set sensible default parameters
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("stop_report_after_feature", 5);
      feature_finder_param.setValue("rt_normalization_factor", 100.0); // for iRT peptides between 0 and 100 (more or less)
      feature_finder_param.setValue("Scores:use_ms1_mi", "true");
      feature_finder_param.setValue("Scores:use_mi_score", "true");

      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_shape_metrics", "false");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");

      // Peak Picker
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false"); // no log messages
      // TODO it seems that the legacy method produces slightly larger peaks, e.g. it will not cut off peaks too early
      // however the same can be achieved by using a relatively low SN cutoff in the -Scoring:TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise 0.5
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn", -1);
      feature_finder_param.setValue("uis_threshold_peak_area", 0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:stop_after_feature");

      // IM scoring: use "auto" so OpenSwathWorkflow can auto-enable for PASEF
      // while still respecting explicit "true"/"false" from the user.
      // Resolved to "true"/"false" before passing to MRMFeatureFinderScoring.
      feature_finder_param.setValue("Scores:use_ion_mobility_scores", "auto",
        "Use ion mobility scores. 'auto' enables for PASEF data and disables otherwise; "
        "'true'/'false' override auto-detection.");
      feature_finder_param.setValidStrings("Scores:use_ion_mobility_scores", {"true", "false", "auto"});

      // EMG Scoring - turn off by default since it is very CPU-intensive
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");

      // remove these parameters
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else if (name == "Calibration")
    {
      // Use CalibrationWorkflow's defaults and add OpenSwathWorkflow-specific parameters
      CalibrationWorkflow cal_wf;
      Param p = cal_wf.getDefaults();
      
      // Add OpenSwathWorkflow-specific parameters that aren't in CalibrationWorkflow
      p.setValue("tr_irt_priority_sampling", "", "Optional custom transition file (TSV format only) containing additional priority peptides for iRT sampling. These peptides will be prioritized alongside the built-in irtkit and cirtkit peptides when `auto_irt` is enabled. Useful for including project-specific or custom iRT peptides.");
      p.setValue("rt_norm", "", "RT normalization file (how to map the RTs of this run to the ones stored in the library). If set, tr_irt may be omitted.");
      
      return p;
    }
    else if (name == "Calibration:RTNormalization")
    {
      Param p;

      p.setValue("alignmentMethod", "linear", "How to perform the alignment to the normalized RT space using anchor points. 'linear': perform linear regression (for few anchor points). 'interpolated': Interpolate between anchor points (for few, noise-free anchor points). 'lowess' Use local regression (for many, noisy anchor points). 'b_spline' use b splines for smoothing.");
      p.setValidStrings("alignmentMethod", {"linear","interpolated","lowess","b_spline"});
      p.setValue("lowess:auto_span", "true", "If true, or if 'span' is 0, automatically select LOWESS span by cross-validation.");
      p.setValidStrings("lowess:auto_span", {"true","false"});
      p.setValue("lowess:span", 0.05, "Span parameter for lowess");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("lowess:auto_span_min", 0.15,"Lower bound for auto-selected span.");
      p.setMinFloat("lowess:auto_span_min", 0.001);
      p.setValue("lowess:auto_span_max", 0.80,"Upper bound for auto-selected span.");
      p.setMaxFloat("lowess:auto_span_max", 0.99);
      p.setValue("lowess:auto_span_grid", "0.005,0.01,0.05,0.15,0.25,0.30,0.50,0.70,0.90", "Optional explicit grid of span candidates in (0,1]. Comma-separated list, e.g. '0.2,0.3,0.5'.  If empty, a default grid is used.");
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline");
      p.setMinInt("b_spline:num_nodes", 0);

      p.setValue("useIterativeChauvenet", "false", "Whether to use Chauvenet's criterion when using iterative methods. This should be used if the algorithm removes too many datapoints but it may lead to true outliers being retained.");
      p.setValidStrings("useIterativeChauvenet", {"true","false"});

      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for the RANSAC outlier detection algorithm.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum threshold in RT dimension for the RANSAC outlier detection algorithm (in percent of the total gradient). Default is set to 3% which is around +/- 4 minutes on a 120 gradient.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size of data points per iteration for the RANSAC outlier detection algorithm.");

      p.setValue("estimateBestPeptides", "false", "Whether the algorithms should try to choose the best peptides based on their peak shape for normalization. Use this option you do not expect all your peptides to be detected in a sample and too many 'bad' peptides enter the outlier removal step (e.g. due to them being endogenous peptides or using a less curated list of peptides).");
      p.setValidStrings("estimateBestPeptides", {"true","false"});

      p.setValue("InitialQualityCutoff", 0.5, "The initial overall quality cutoff for a peak to be scored (range ca. -2 to 2)");
      p.setValue("OverallQualityCutoff", 5.5, "The overall quality cutoff for a peak to go into the retention time estimation (range ca. 0 to 10)");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage. This option should be used to ensure that there is a complete coverage of the RT space (this should detect cases where only a part of the RT gradient is actually covered by normalization peptides)");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides that are required for a bin to counted as 'covered'");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered");
      return p;
    }
    else if (name == "MRMMapping")
    {
      Param p;

      p.setValue("precursor_tolerance", 0.9, "Precursor tolerance when mapping (in Th)");
      p.setValue("product_tolerance", 1.2, "Product tolerance when mapping (in Th)");

      p.setValue("irt_precursor_tolerance", 1.5, "Precursor tolerance when mapping iRT transitions (in Th)");
      p.setValue("irt_product_tolerance", 1.5, "Product tolerance when mapping iRT transitions (in Th)");
      
      p.setValue("map_multiple_assays", "false", "Allow to map multiple assays to chromatograms and duplicate these chromatograms in the output.");
      p.setValidStrings("map_multiple_assays", {"true","false"});
      
      p.setValue("error_on_unmapped", "false", "Treat remaining, unmapped chromatograms as an error");
      p.setValidStrings("error_on_unmapped", {"true","false"});

      return p;
    }
    else if (name == "Calibration:MassIMCorrection")
    {
      Param p = SwathMapMassCorrection().getDefaults();
      return p;
    }
    else
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
    }
  }

  /**
    @brief Load priority peptide sequences from TSV files (irtkit and cirtkit)
    
    Loads peptide sequences from the specified TSV files and returns them as a
    set for quick lookup. Used to prioritize common iRT peptides during sampling.
    
    @param[in] tsv_files Vector of file paths to TSV files to load
    @param[in] tsv_reader_param Parameters for the TSV reader
    
    @return Set of unique peptide sequences from the loaded files
  */
  std::unordered_set<std::string> loadPriorityPeptideSequences(
    const std::vector<String>& tsv_files,
    const Param& tsv_reader_param)
  {
    std::unordered_set<std::string> priority_sequences;
    
    for (const auto& tsv_file : tsv_files)
    {
      if (tsv_file.empty() || !File::exists(tsv_file))
      {
        OPENMS_LOG_WARN << "Priority peptide file not found: " << tsv_file << std::endl;
        continue;
      }
      
      try
      {
        FileTypes::Type file_type = FileHandler::getType(tsv_file);
        OpenSwath::LightTargetedExperiment priority_exp = loadTransitionList(file_type, tsv_file, tsv_reader_param);
        
        for (const auto& compound : priority_exp.getCompounds())
        {
          if (!compound.sequence.empty())
          {
            priority_sequences.insert(compound.sequence);
          }
        }
        
        OPENMS_LOG_INFO << "Loaded " << priority_exp.getCompounds().size() 
                        << " compounds from priority file: " << tsv_file << std::endl;
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Failed to load priority peptide file " << tsv_file 
                        << ": " << e.what() << std::endl;
      }
    }
    
    OPENMS_LOG_INFO << "Total unique priority peptide sequences: " 
                    << priority_sequences.size() << std::endl;
    
    return priority_sequences;
  }

  ExitCodes main_(int, const char **) override
  {
    ///////////////////////////////////
    // Prepare Parameters
    ///////////////////////////////////
    StringList file_list = getStringList_("in");
    String tr_file = getStringOption_("tr");
    String out_features = getStringOption_("out_features");

    //tr_file input file type
    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = FileHandler::getType(tr_file);
      writeDebug_(String("Input file type (-tr): ") + FileTypes::typeToName(tr_type), 2);
    }

    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type for '-tr' !");
      return PARSE_ERROR;
    }

    //tr_file input file type
    FileTypes::Type out_features_type = FileTypes::nameToType(getStringOption_("out_features_type"));
    if (out_features_type == FileTypes::UNKNOWN)
    {
      out_features_type = FileHandler::getType(out_features);
      writeDebug_(String("Input file type (-out): ") + FileTypes::typeToName(out_features_type), 2);
    }

    if (out_features_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type for '-out_features' !");
      return PARSE_ERROR;
    }
    String out_qc = getStringOption_("out_qc");

    Param irt_calibration_params = getParam_().copy("Calibration:", true);
    bool auto_irt = (irt_calibration_params.getValue("auto_irt:enabled").toString() == "true");

    // Extract only the parameters needed for OpenSwathWorkflow-specific validation and logic
    String irt_tr_file = irt_calibration_params.getValue("files:linear_irt_file").toString();
    String priority_sampling_irt_tr_file = irt_calibration_params.getValue("tr_irt_priority_sampling").toString();
    String trafo_in = irt_calibration_params.getValue("rt_norm").toString();
    
    // Extract parameters needed for OpenSwathWorkflow validation logic
    UInt irt_bins_lin = irt_calibration_params.getValue("auto_irt:irt_bins");
    UInt irt_pep_lin  = irt_calibration_params.getValue("auto_irt:irt_peptides_per_bin");

    // If a linear iRT file is explicitly provided, auto_irt must be disabled
    // and any priority sampling iRT file should be ignored.
    if (!irt_tr_file.empty())
    {
      if (auto_irt || !priority_sampling_irt_tr_file.empty())
      {
        OPENMS_LOG_WARN << "Calibration:files:linear_irt_file provided -> disabling auto_irt and ignoring tr_irt_priority_sampling." << std::endl;
      }
      auto_irt = false;
      irt_calibration_params.setValue("auto_irt:enabled", "false");
      // clear the priority sampling file so downstream logic won't attempt to use/validate it
      priority_sampling_irt_tr_file.clear();
      irt_calibration_params.setValue("tr_irt_priority_sampling", "");
    }
    
    String swath_windows_file = getStringOption_("swath_windows_file");

    String out_chrom = getStringOption_("out_chrom");
    String out_mobilogram = getStringOption_("out_mobilogram");
    bool split_file = getFlag_("split_file_input");
    bool use_emg_score = getFlag_("use_elution_model_score");
    bool force = getFlag_("force");
    bool pasef = getFlag_("pasef");
    bool sort_swath_maps = getFlag_("sort_swath_maps");
    bool use_ms1_traces = getStringOption_("enable_ms1") == "true";
    bool enable_uis_scoring = getStringOption_("enable_ipf") == "true";
    int batchSize = (int)getIntOption_("batchSize");
    int outer_loop_threads = (int)getIntOption_("outer_loop_threads");
    int ms1_isotopes = (int)getIntOption_("ms1_isotopes");
    Size debug_level = (Size)getIntOption_("debug");

    Param debug_params = getParam_().copy("Debugging:", true);

    StringList disable_features = getStringList_("Debugging:disable_features");
    bool disable_im_calibration = std::find(disable_features.begin(), disable_features.end(), "no_IM_calibration") != disable_features.end();
    bool disable_im_windowing   = std::find(disable_features.begin(), disable_features.end(), "no_IM_windowing")   != disable_features.end();

    String readoptions = getStringOption_("readOptions");
    bool keep_cached_files = getFlag_("keep_cached_files");

    // make sure tmp is a directory with proper separator at the end (downstream methods simply do path + filename)
    // File::absolutePath() always uses '/' separators
    String tmp_dir = File::absolutePath(getStringOption_("tempDirectory")).ensureLastChar('/');

    ///////////////////////////////////
    // Parameter validation
    ///////////////////////////////////

    bool load_into_memory = false;
    if (readoptions == "cacheWorkingInMemory")
    {
      readoptions = "cache";
      load_into_memory = true;
    }
    else if (readoptions == "workingInMemory")
    {
      readoptions = "normal";
      load_into_memory = true;
    }

    bool is_sqmass_input  = (FileHandler::getTypeByFileName(file_list[0]) == FileTypes::SQMASS);
    if (is_sqmass_input && !load_into_memory)
    {
      std::cout << "When using sqMass input files, it is highly recommended to use the workingInMemory option as otherwise data access will be very slow." << std::endl;
    }

    if (trafo_in.empty() && irt_tr_file.empty() && !auto_irt)
    {
      std::cout << "Since neither rt_norm nor tr_irt nor auto_irt is set, OpenSWATH will " <<
        "not use RT-transformation (rather a null transformation will be applied)" << std::endl;
    }

    // -----------------------------------------------------------------
    // Validate auto_irt parameters
    // -----------------------------------------------------------------
    if (auto_irt)
    {
      // linear sampling must have at least one bin and one peptide per bin
      if (irt_bins_lin == 0)
      {
        writeLogError_("Parameter error: --irt_bins must be > 0 when auto_irt is enabled.");
        return PARSE_ERROR;
      }
      if (irt_pep_lin == 0)
      {
        writeLogError_("Parameter error: --irt_peptides_per_bin must be > 0 when auto_irt is enabled.");
        return PARSE_ERROR;
      }
    }
    
    // Validate priority iRT sampling file format if provided and if auto_irt is enabled
    if (!priority_sampling_irt_tr_file.empty())
    {
      if (!auto_irt)
      {
        // auto_irt disabled (possibly due to explicit linear iRT file); ignore provided priority sampling file
        OPENMS_LOG_WARN << "Priority iRT sampling file provided but auto_irt is disabled; ignoring: " << priority_sampling_irt_tr_file << std::endl;
      }
      else
      {
        if (!File::exists(priority_sampling_irt_tr_file))
        {
          writeLogError_("Parameter error: Priority iRT file does not exist: " + priority_sampling_irt_tr_file);
          return PARSE_ERROR;
        }
        
        FileTypes::Type priority_file_type = FileHandler::getType(priority_sampling_irt_tr_file);
        if (priority_file_type != FileTypes::TSV)
        {
          writeLogError_("Parameter error: Priority iRT file must be in TSV format. Provided: " + 
                         FileTypes::typeToName(priority_file_type));
          return PARSE_ERROR;
        }
      }
    }

    // Check swath window input
    if (!swath_windows_file.empty())
    {
      OPENMS_LOG_INFO << "Validate provided Swath windows file:" << std::endl;
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);

      for (Size i = 0; i < swath_prec_lower.size(); i++)
      {
        OPENMS_LOG_DEBUG << "Read lower swath window " << swath_prec_lower[i] << " and upper window " << swath_prec_upper[i] << std::endl;
      }
    }

    double min_upper_edge_dist = getDoubleOption_("min_upper_edge_dist");
    bool use_ms1_im = getStringOption_("use_ms1_ion_mobility") == "true";
    bool prm = getStringOption_("matching_window_only") == "true";

    ChromExtractParams cp;
    cp.min_upper_edge_dist   = min_upper_edge_dist;
    cp.mz_extraction_window  = getDoubleOption_("mz_extraction_window");
    cp.ppm                   = getStringOption_("mz_extraction_window_unit") == "ppm";
    cp.rt_extraction_window  = getDoubleOption_("rt_extraction_window");
    cp.im_extraction_window  = getDoubleOption_("ion_mobility_window");
    cp.extraction_function   = getStringOption_("extraction_function");
    cp.extra_rt_extract      = getDoubleOption_("extra_rt_extraction_window");

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1; // extract the whole RT range for iRT measurements
    cp_irt.mz_extraction_window = getDoubleOption_("irt_mz_extraction_window");
    cp_irt.im_extraction_window = getDoubleOption_("irt_im_extraction_window");

    if (disable_im_calibration)
    {
      cp_irt.im_extraction_window = -1;
      OPENMS_LOG_INFO << "Debugging: no_IM_calibration active -- forcing irt_im_extraction_window = -1 (skipping IM calibration)." << std::endl;
    }
    else if ( (cp_irt.im_extraction_window == -1) && (cp.im_extraction_window != -1) )
    {
      OPENMS_LOG_WARN << "Warning: -irt_im_extraction_window is not set, this will lead to no ion mobility calibration" << std::endl;
    }

    cp_irt.ppm                  = getStringOption_("irt_mz_extraction_window_unit") == "ppm";

    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window  = getDoubleOption_("mz_extraction_window_ms1");
    cp_ms1.ppm                   = getStringOption_("mz_extraction_window_ms1_unit") == "ppm";
    cp_ms1.im_extraction_window  = (use_ms1_im) ? getDoubleOption_("im_extraction_window_ms1") : -1;

    if (disable_im_windowing)
    {
      cp.im_extraction_window = -1;
      cp_ms1.im_extraction_window = -1;
      OPENMS_LOG_INFO << "Debugging: no_IM_windowing active -- forcing im_extraction_window = -1 for MS2 and MS1 (using full IM range)." << std::endl;
    }

    Param feature_finder_param = getParam_().copy("Scoring:", true);
    feature_finder_param.setValue("use_ms1_ion_mobility", getStringOption_("use_ms1_ion_mobility"));


    Param tsv_reader_param = getParam_().copy("Library:", true);
    if (use_emg_score)
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "true");
    }
    else
    {
      feature_finder_param.setValue("Scores:use_elution_model_score", "false");
    }
    if (use_ms1_traces)
    {
      feature_finder_param.setValue("Scores:use_ms1_correlation", "true");
      feature_finder_param.setValue("Scores:use_ms1_fullscan", "true");
    }
    if (enable_uis_scoring)
    {
      feature_finder_param.setValue("Scores:use_uis_scores", "true");
    }

    bool compute_peak_shape_metrics = feature_finder_param.getValue("TransitionGroupPicker:compute_peak_shape_metrics").toBool();
    if (compute_peak_shape_metrics)
    {
      feature_finder_param.setValue("Scores:use_peak_shape_metrics", "true");
    }

    // Extract MRMMapping-related options from the CLI subsection and pass them to the workflow
    Param tmp_mrm_map_param = getParam_().copy("MRMMapping:", true);
    Param irt_mrm_map_param = OpenMS::MRMMapping().getDefaults();
    irt_mrm_map_param.setValue("precursor_tolerance", tmp_mrm_map_param.getValue("irt_precursor_tolerance"));
    irt_mrm_map_param.setValue("product_tolerance", tmp_mrm_map_param.getValue("irt_product_tolerance"));
    irt_mrm_map_param.setValue("map_multiple_assays", tmp_mrm_map_param.getValue("map_multiple_assays"));
    irt_mrm_map_param.setValue("error_on_unmapped", tmp_mrm_map_param.getValue("error_on_unmapped"));

    Param mrm_map_param = OpenMS::MRMMapping().getDefaults();
    mrm_map_param.setValue("precursor_tolerance", tmp_mrm_map_param.getValue("precursor_tolerance"));
    mrm_map_param.setValue("product_tolerance", tmp_mrm_map_param.getValue("product_tolerance"));
    mrm_map_param.setValue("map_multiple_assays", tmp_mrm_map_param.getValue("map_multiple_assays"));
    mrm_map_param.setValue("error_on_unmapped", tmp_mrm_map_param.getValue("error_on_unmapped"));

    ///////////////////////////////////
    // Load the transitions
    ///////////////////////////////////
    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(tr_type, tr_file, tsv_reader_param);
    OPENMS_LOG_INFO << "Loaded " << transition_exp.getProteins().size() << " proteins, " <<
      transition_exp.getCompounds().size() << " compounds with " << transition_exp.getTransitions().size() << " transitions." << std::endl;

    if (out_features_type == FileTypes::OSW)
    {
      if (tr_type == FileTypes::PQP)
      {
         // copy the PQP file and name it OSW file
          std::ifstream  src(tr_file.c_str(), std::ios::binary);
          std::ofstream  dst(out_features.c_str(), std::ios::binary | std::ios::trunc);
          dst << src.rdbuf();
      }
      else if (tr_type == FileTypes::TSV)
      {
        // Convert TSV to .PQP 
        TransitionTSVFile tsv_reader;
        TargetedExperiment transition_exp_heavy;
        tsv_reader.setParameters(tsv_reader_param);
        tsv_reader.convertTSVToTargetedExperiment(tr_file.c_str(), tr_type, transition_exp_heavy);
        TransitionPQPFile().convertTargetedExperimentToPQP(out_features.c_str(), transition_exp_heavy);

        // instead of reloading - edit the already loaded transition_exp to be compatible with .pqp format
        // read the PQP to traMLID mapping
        auto precursor_traml_to_pqp = TransitionPQPFile().getPQPIDToTraMLIDMap(out_features.c_str(), "PRECURSOR");
        auto transition_traml_to_pqp = TransitionPQPFile().getPQPIDToTraMLIDMap(out_features.c_str(), "TRANSITION");

        // convert tramlID in transitionExp to PQP ID
        for (auto & prec : transition_exp.getCompounds())
        {
          if (auto id = precursor_traml_to_pqp.find(prec.id); id != precursor_traml_to_pqp.end())
          {
            prec.id = id->second;
          }
        }

        for (auto & tr : transition_exp.getTransitions())
        {
          // convert transition tramlID peptide reference in transitionExp to PQP ID 
          auto pep = precursor_traml_to_pqp.find(tr.getPeptideRef());
          if (pep != precursor_traml_to_pqp.end())
          {
            tr.peptide_ref = pep->second;
          }

          // Update transition id
          auto id = transition_traml_to_pqp.find(tr.transition_name);
          if (id != transition_traml_to_pqp.end())
          {
            tr.transition_name = id->second;
          }
        }
      }
      else if (tr_type == FileTypes::OSWPQ)
      {
        // Convert parquet library to .PQP for OSW output
        TransitionPQPFile().convertLightTargetedExperimentToPQP(out_features.c_str(), transition_exp);

        auto precursor_traml_to_pqp = TransitionPQPFile().getPQPIDToTraMLIDMap(out_features.c_str(), "PRECURSOR");
        auto transition_traml_to_pqp = TransitionPQPFile().getPQPIDToTraMLIDMap(out_features.c_str(), "TRANSITION");

        for (auto & prec : transition_exp.getCompounds())
        {
          if (auto id = precursor_traml_to_pqp.find(prec.id); id != precursor_traml_to_pqp.end())
          {
            prec.id = id->second;
          }
        }

        for (auto & tr : transition_exp.getTransitions())
        {
          auto pep = precursor_traml_to_pqp.find(tr.getPeptideRef());
          if (pep != precursor_traml_to_pqp.end())
          {
            tr.peptide_ref = pep->second;
          }

          auto id = transition_traml_to_pqp.find(tr.transition_name);
          if (id != transition_traml_to_pqp.end())
          {
            tr.transition_name = id->second;
          }
        }
      }
      else if (tr_type == FileTypes::TRAML)
      {
        if (out_features_type == FileTypes::OSW)
        {
          throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, String("Conversion from TraML to OSW is not supported."));
        }
      }
    }

    // Load priority peptide sequences from irtkit and cirtkit if auto_irt is enabled
    std::unordered_set<std::string> priority_peptides;
    if (auto_irt)
    {
      String data_path = File::getOpenMSDataPath();
      std::vector<String> priority_files;
      
      String irtkit_path = data_path + "/CHEMISTRY/irtkit.tsv";
      String cirtkit_path = data_path + "/CHEMISTRY/cirtkit.tsv";
      
      if (File::exists(irtkit_path))
      {
        priority_files.push_back(irtkit_path);
      }
      else
      {
        OPENMS_LOG_WARN << "irtkit.tsv not found at: " << irtkit_path << std::endl;
      }
      
      if (File::exists(cirtkit_path))
      {
        priority_files.push_back(cirtkit_path);
      }

      // Check if tr_irt_priority_sampling is enabled
      if (!priority_sampling_irt_tr_file.empty())
      {
        priority_files.push_back(priority_sampling_irt_tr_file);
      }
      
      if (!priority_files.empty())
      {
        Param priority_tsv_param = TransitionTSVFile().getDefaults();
        priority_peptides = loadPriorityPeptideSequences(priority_files, priority_tsv_param);
      }
      else
      {
        OPENMS_LOG_WARN << "No priority peptide files found. Continuing without priority sampling." << std::endl;
      }
    }

    OPENMS_LOG_INFO << std::endl;

    ///////////////////////////////////
    // Create run groups based on split_file_input flag
    ///////////////////////////////////

    std::vector<StringList> run_groups;
    if (split_file)
    {
      // split_file_input mode: one run with all files grouped together (split-window)
      run_groups.push_back(file_list);
    }
    else
    {
      // multi-run mode: N runs, each with single file
      for (const String& file : file_list)
      {
        run_groups.push_back({file});
      }
    }

    ///////////////////////////////////
    // Iterate over each run group and process individually
    ///////////////////////////////////

    // Set up shared output objects that persist across files
    FeatureMap out_featureFile;  // accumulates features across all files
    const bool write_osw = (out_features_type == FileTypes::OSW);
    const bool write_parquet = (out_features_type == FileTypes::OSWPQ);
    String osw_out_filename = write_osw ? out_features : "";
    OpenSwathOSWWriter oswwriter(osw_out_filename, enable_uis_scoring);

    String parquet_dir = out_features;
    bool parquet_zip_output = false;
    std::unique_ptr<File::TempDir> parquet_temp_dir;
    OpenSwathOSWParquetWriter parquet_writer;
    // Configure writer append behavior from CLI flag
    parquet_writer.setPreserveExisting(getFlag_("append_oswpq"));
    if (write_parquet)
    {
      parquet_zip_output = out_features.hasSuffix(".oswpq");
      if (parquet_zip_output)
      {
        if (getFlag_("append_oswpq") && File::exists(out_features))
        {
          // Extract existing archive so prior run data is preserved when appending.
          parquet_dir = ZipArchiveFile::unzipDirectory(out_features, parquet_temp_dir);
        }
        else
        {
          parquet_temp_dir = std::make_unique<File::TempDir>();
          parquet_dir = parquet_temp_dir->getPath() + "/oswpq_output";
          File::makeDir(parquet_dir);
        }
      }
    }

    // Write DB schema once (only for first file)
    if (write_osw)
    {
      oswwriter.writeHeader();
    }

    const bool user_pasef = pasef;

    Size run_index = 0;
    for (const StringList& current_run_files : run_groups)
    {
      OPENMS_LOG_INFO << "Processing Run " << (run_index + 1) << "/" << run_groups.size() << std::endl;
      
      // Create fresh copies of parameters for each run to avoid carry-over
      pasef = user_pasef;
      ChromExtractParams cp_current = cp;
      ChromExtractParams cp_ms1_current = cp_ms1;
      ChromExtractParams cp_irt_current = cp_irt;
      Param feature_finder_param_run = feature_finder_param;
      
      ///////////////////////////////////
      // Per-run temporary cache directory (created only when using cache readOptions)
      // Use File::TempDir for RAII-based cleanup: destructor removes dir (unless keep_cached_files is true)
      String per_run_tmp = tmp_dir;
      std::unique_ptr<File::TempDir> per_run_temp_dir;
      if (readoptions == "cache")
      {
        per_run_temp_dir = std::make_unique<File::TempDir>(tmp_dir, keep_cached_files);
        per_run_tmp = per_run_temp_dir->getPath();
      }

      // Load the SWATH files (if split data, otherwise load single experiment mzML)
      ///////////////////////////////////
      std::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
      std::vector< OpenSwath::SwathMap > swath_maps;
      std::vector<String> swath_map_sources;

      StringList single_file_list = current_run_files;

      // collect some QC data (only for the first run if multiple runs)
      if (!out_qc.empty() && run_index == 0)
      {
        OpenSwath::SwathQC qc(30, 0.04);
        MSDataTransformingConsumer qc_consumer; // apply some transformation
        qc_consumer.setSpectraProcessingFunc(qc.getSpectraProcessingFunc());
        qc_consumer.setExperimentalSettingsFunc(qc.getExpSettingsFunc());
        if (!loadSwathFiles(single_file_list, exp_meta, swath_maps, swath_map_sources, split_file, per_run_tmp, readoptions,
                swath_windows_file, min_upper_edge_dist, force,
                sort_swath_maps, prm, &qc_consumer))
        {
          OPENMS_LOG_ERROR << "Failed to load SWATH files for Run " << (run_index + 1)
                           << ": " << ListUtils::concatenate(single_file_list, ", ") << std::endl
                           << "  (split_file=" << (split_file ? "true" : "false") << ", tmp_dir='" << tmp_dir << "', readoptions='" << readoptions << "')"
                           << std::endl
                           << "  Please check that the input files exist, are valid mzML/sqMass files, and that read options are correct." << std::endl;
          return PARSE_ERROR;
        }
        qc.storeJSON(out_qc);
      }
      else
      {
        if (!loadSwathFiles(single_file_list, exp_meta, swath_maps, swath_map_sources, split_file, per_run_tmp, readoptions,
                swath_windows_file, min_upper_edge_dist, force,
                sort_swath_maps, prm))
        {
          OPENMS_LOG_ERROR << "Failed to load SWATH files for Run " << (run_index + 1)
                           << ": " << ListUtils::concatenate(single_file_list, ", ") << std::endl
                           << "  (split_file=" << (split_file ? "true" : "false") << ", tmp_dir='" << tmp_dir << "', readoptions='" << readoptions << "')"
                           << std::endl
                           << "  Please check that the input files exist, are valid mzML/sqMass files, and that read options are correct." << std::endl;
          return PARSE_ERROR;
        }
      }

      // Auto-detect PASEF data from swath maps if user didn't explicitly set -pasef
      if (!pasef)
      {
        bool has_im_windows = std::any_of(swath_maps.begin(), swath_maps.end(),
          [](const OpenSwath::SwathMap& m) { return !m.ms1 && m.imLower >= 0 && m.imUpper >= 0; });
        if (has_im_windows)
        {
          OPENMS_LOG_INFO << "Auto-detected ion mobility (PASEF) data from SWATH windows. "
                          << "Enabling PASEF mode automatically." << std::endl;
          pasef = true;
        }
      }

      // Override: disable PASEF window matching if no_IM_windowing is active
      if (disable_im_windowing && pasef)
      {
        OPENMS_LOG_INFO << "Debugging: no_IM_windowing active -- forcing pasef = false (disabling PASEF window matching)." << std::endl;
        pasef = false;
      }

      // Resolve "auto" for ion mobility scoring: enable for PASEF data, disable otherwise.
      // Explicit "true"/"false" from the user is always respected.
      {
        String im_score_setting = feature_finder_param_run.getValue("Scores:use_ion_mobility_scores").toString();
        if (im_score_setting == "auto")
        {
          if (pasef)
          {
            OPENMS_LOG_INFO << "PASEF mode active: automatically enabling ion mobility scores "
                            << "(Scoring:Scores:use_ion_mobility_scores = auto -> true)." << std::endl;
            feature_finder_param_run.setValue("Scores:use_ion_mobility_scores", "true");
          }
          else
          {
            feature_finder_param_run.setValue("Scores:use_ion_mobility_scores", "false");
            if (disable_im_windowing)
            {
              OPENMS_LOG_INFO << "IM scoring auto-disabled because no_IM_windowing forced pasef = false. "
                              << "Pass -Scoring:Scores:use_ion_mobility_scores true to keep IM scoring." << std::endl;
            }
          }
        }
      }

      // Validate that transitions have IM values when PASEF mode is active
      if (pasef)
      {
        for (const auto& tr : transition_exp.getTransitions())
        {
          if (tr.precursor_im < 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Error: Transition " + tr.getNativeID() +
              " does not have a valid IM value, this must be set to use PASEF mode (auto-detected or via -pasef flag)");
          }
        }
      }

      ///////////////////////////////////
      // Get the transformation information (using iRT peptides) - per file
      ///////////////////////////////////

      // Create a basename for this run's outputs (used for multi-run scenarios)
      String file_basename;
      if (run_groups.size() > 1)
      {
        // Extract basename from input file path (remove directory and extension)
        // For multi-run mode, each run has exactly one file
        String filename = File::basename(current_run_files[0]);
        Size dot_pos = filename.find_last_of('.');
        if (dot_pos != String::npos)
        {
          file_basename = filename.substr(0, dot_pos);
        }
        else
        {
          file_basename = filename;
        }
      }

      String irt_trafo_out = debug_params.getValue("irt_trafo").toString();
      if (!irt_trafo_out.empty() && run_groups.size() > 1)
      {
        // For multi-run, use basename prefix to make unique filenames
        String base_name = irt_trafo_out.substr(0, irt_trafo_out.find_last_of('.'));
        String extension = irt_trafo_out.substr(irt_trafo_out.find_last_of('.'));
        irt_trafo_out = file_basename + "_" + base_name + extension;
      }

      String irt_mzml_out = debug_params.getValue("irt_mzml").toString();
      if (!irt_mzml_out.empty() && run_groups.size() > 1)
      {
        // For multi-run, use basename prefix to make unique filenames
        String base_name = irt_mzml_out.substr(0, irt_mzml_out.find_last_of('.'));
        String extension = irt_mzml_out.substr(irt_mzml_out.find_last_of('.'));
        irt_mzml_out = file_basename + "_" + base_name + extension;
      }

      Param irt_detection_param = getParam_().copy("Calibration:RTNormalization:", true);
      Param calibration_param = getParam_().copy("Calibration:MassIMCorrection:", true);
      calibration_param.setValue("mz_extraction_window", cp_irt_current.mz_extraction_window);
      calibration_param.setValue("mz_extraction_window_ppm", cp_irt_current.ppm ? "true" : "false");
      calibration_param.setValue("im_extraction_window", cp_irt_current.im_extraction_window);

      // Detect SRM/MRM mode: check if all swath_maps are chromatogram-only (no spectra, not MS1)
      
      bool mrm_mode = true;
      for (const auto& sm : swath_maps)
      {
        if (sm.ms1 || sm.sptr->getNrSpectra() > 0)
        {
          mrm_mode = false;
          break;
        }
      }

      ///////////////////////////////////
      // RT transformation: Load existing trafo or perform calibration
      ///////////////////////////////////
      
      TransformationDescription trafo_rtnorm;
      
      if (!trafo_in.empty())
      {
        // Load existing RT transformation file using FileHandler so any metadata is preserved
        // and fit the selected model using the RTNormalization parameters so behavior
        // matches the original implementation.
        FileHandler().loadTransformations(trafo_in, trafo_rtnorm, false, {FileTypes::TRANSFORMATIONXML});

        // Prepare model parameters from RTNormalization section
        Param model_params;
        model_params.setValue("symmetric_regression", "false");
        // copy span and num_nodes from the RTNormalization detection params
        model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
        model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
        String model_type = irt_detection_param.getValue("alignmentMethod").toString();

        // Fit the model to the loaded transformation
        trafo_rtnorm.fitModel(model_type, model_params);

        // Note: When using existing RT transformation, no m/z or IM calibration is performed
        OPENMS_LOG_WARN << "Using existing RT transformation - which has no m/z and ion mobility calibration" << std::endl;
      }
      else
      {      
        // Setup CalibrationWorkflow configuration from TOPP parameters
        CalibrationWorkflow calibration_wf;
        
        // Filter parameters to exclude MassIMCorrection and RTNormalization parameters, these are passed as their own separate parameter objects
        // Pass everything else to CalibrationWorkflow
        Param cal_params = irt_calibration_params.copy("", false);
        
        // Remove the sections that don't belong to CalibrationWorkflow
        cal_params.remove("MassIMCorrection:");
        cal_params.remove("RTNormalization:");
        // Top-level parameters handled by OpenSwathWorkflow
        cal_params.remove("rt_norm");  
        cal_params.remove("tr_irt_priority_sampling");  
        
        calibration_wf.setParameters(cal_params);
        calibration_wf.setLogType(log_type_);
        
        // Determine iRT strategy based on configured parameters and multi-run context
        IrtStrategy strategy = calibration_wf.determineIrtStrategy(
          transition_exp, run_groups.size());
        
        // Prepare iRT experiments for this run
        std::vector<String> priority_pep_strings;
        priority_pep_strings.reserve(priority_peptides.size());
        for (const auto& pep : priority_peptides)
        {
          priority_pep_strings.push_back(String(pep));
        }
        
        CalibrationWorkflow::IrtExperiments irt_experiments = calibration_wf.prepareIrtExperiments(
          strategy, transition_exp, priority_pep_strings, run_index);
        
        // Single modular calibration call - handles all scenarios  
        auto calibration_result = calibration_wf.performCalibration(
          swath_maps,         // SWATH data maps  
          transition_exp,     // Target transition experiment (may be modified for IM)
          cp_current,         // Extraction parameters (may be updated with estimates)
          cp_ms1_current,     // MS1 extraction parameters (may be updated)
          irt_experiments,    // Pre-prepared iRT experiments
          feature_finder_param_run,  // Feature finder parameters (per-run copy)
          cp_irt_current,     // iRT extraction parameters
          irt_detection_param, // iRT detection parameters
          calibration_param,  // Calibration parameters (m/z, IM correction)
          irt_mrm_map_param,  // MRM mapping parameters
          pasef,              // PASEF data flag
          load_into_memory,   // Load data into memory flag
          irt_trafo_out,      // Transformation output file
          irt_mzml_out,       // iRT chromatograms output file
          debug_level         // Debug level for iRT chromatogram output
        );
        
        // Extract results
        trafo_rtnorm = calibration_result.rt_trafo;
      }

    ///////////////////////////////////
    // perform extraction on current file
    ///////////////////////////////////

    // Create one run ID per input run and reuse it across all run-level outputs.
    UInt64 cur_run = OpenMS::UniqueIdGenerator::getUniqueId();

    // Set up chromatogram output for this file
    // Either use chrom.mzML or sqliteDB (sqMass)
    Interfaces::IMSDataConsumer* chromatogramConsumer;
    String out_chrom_current = out_chrom;
    if (!out_chrom.empty() && run_groups.size() > 1)
    {
      // Preserve parent directory when creating per-run filenames.
      // Split path and filename first, then prepend the run prefix to the filename only.
      String parent = File::path(out_chrom);
      String filename = File::basename(out_chrom);
      String stem = filename.substr(0, filename.find_last_of('.'));
      String extension = filename.substr(filename.find_last_of('.'));
      String fname_with_prefix = file_basename + "_" + stem + extension;
      out_chrom_current = (parent == "." ? fname_with_prefix : parent + "/" + fname_with_prefix);
    }
    prepareChromOutput(&chromatogramConsumer, exp_meta, transition_exp, out_chrom_current, cur_run, current_run_files[0]);

    // Prepare mobilogram output (per-run)
    std::unique_ptr<MobilogramParquetConsumer> mobilogramConsumer;
    String out_mobilogram_current = out_mobilogram;
    if (!out_mobilogram.empty() && run_groups.size() > 1)
    {
      // Preserve parent directory when creating per-run filenames.
      // Split path and filename first, then prepend the run prefix to the filename only.
      String parent = File::path(out_mobilogram);
      String filename = File::basename(out_mobilogram);
      String stem = filename.substr(0, filename.find_last_of('.'));
      String extension = filename.substr(filename.find_last_of('.'));
      String fname_with_prefix = file_basename + "_" + stem + extension;
      out_mobilogram_current = (parent == "." ? fname_with_prefix : parent + "/" + fname_with_prefix);
    }
    prepareMobilogramOutput(mobilogramConsumer, exp_meta, transition_exp, out_mobilogram_current, cur_run, current_run_files[0]);

    // Register the same run ID in OSW.
    // For OSW, use the first file in the run group as the representative filename
    if (write_osw)
    {
      oswwriter.addRun(cur_run, current_run_files[0]);
    }
    // Also register run in chromatogram consumer if it is a SQL consumer
    MSDataSqlConsumer* sql_cons = dynamic_cast<MSDataSqlConsumer*>(chromatogramConsumer);
    if (sql_cons != nullptr)
    {
      sql_cons->addRun(current_run_files[0], cur_run);
    }

    // set current run id in writer
    oswwriter.setRunId(cur_run);
    // set current run id for sqMass writer as well (reuse previous cast)
    if (sql_cons != nullptr)
    {
      sql_cons->setRunId(cur_run);
    }

    FeatureMap run_featureFile;
    FeatureMap& active_feature_map = write_parquet ? run_featureFile : out_featureFile;

    OpenSwathWorkflow wf(use_ms1_traces, use_ms1_im, prm, pasef, mrm_mode, outer_loop_threads);
    wf.setLogType(log_type_);

    // perform extraction for this file's swath maps
    wf.performExtraction(swath_maps, trafo_rtnorm, cp_current, cp_ms1_current, feature_finder_param_run, transition_exp,
             active_feature_map, true, oswwriter, chromatogramConsumer, batchSize, ms1_isotopes, load_into_memory, mrm_map_param, mobilogramConsumer.get());

    if (mobilogramConsumer)
    {
      mobilogramConsumer->finalize();
      mobilogramConsumer.reset();
    }

    //// Write out data

    delete chromatogramConsumer;

    if (write_parquet)
    {
      parquet_writer.write(parquet_dir, transition_exp, active_feature_map,
                           cur_run, current_run_files[0], enable_uis_scoring);
    }

    OPENMS_LOG_INFO << std::endl;

    ++run_index;
    } // end for each run

    if (write_parquet && parquet_zip_output)
    {
      // Stream files into the zip archive instead of unzipping/rezipping the
      // whole directory. This uses ZipArchiveFile::addOrReplaceFromFile which
      // streams from disk and avoids loading large parquet blobs into memory.
      const std::filesystem::path dirpath = std::filesystem::path(std::string(parquet_dir));
      const String output_zip_abs = File::absolutePath(out_features);
      if (File::exists(output_zip_abs))
      {
        File::remove(output_zip_abs);
      }

      for (auto it = std::filesystem::recursive_directory_iterator(dirpath); it != std::filesystem::recursive_directory_iterator(); ++it)
      {
        if (it->is_directory()) continue;
        const auto full = it->path();
        std::string rel = std::filesystem::relative(full, dirpath).generic_string();
        ZipArchiveFile::addOrReplaceFromFile(out_features, String(rel), String(full.string()));
      }
      // Write the embedded sidecar index that enables random-access reads
      // directly from the archive without extracting (RAF pattern).
      ZipArchiveFile::writeSidecarIndex(output_zip_abs);
    }

    if ( out_features_type == FileTypes::FEATUREXML )
    {
      std::cout << "Writing features ..." << std::endl;
      addDataProcessing_(out_featureFile, getProcessingInfo_(DataProcessing::QUANTITATION));
      out_featureFile.ensureUniqueId();
      FileHandler().storeFeatures(out_features, out_featureFile, {FileTypes::FEATUREXML});
    }

    return EXECUTION_OK;
  }

};

int main(int argc, const char ** argv)
{
  TOPPOpenSwathWorkflow tool;
  return tool.main(argc, argv);
}

/// @endcond
