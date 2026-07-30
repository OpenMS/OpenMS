// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/OpenSwathBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/CalibrationWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathExportConfig.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathGeneInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathLibraryPreparation.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptideInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPercolatorScoring.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathProteinInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathWorkflow.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SwathMapMassCorrection.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionListEvidenceFilter.h>
#include <OpenMS/ANALYSIS/TARGETED/MRMMapping.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <unordered_set>
#include <vector>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_OpenDIA OpenDIA

@brief High-level DIA workflow front-end that prepares peptide query parameters, performs targeted extraction, rescoring, inference, and export in one invocation.

Three library modes:
- @c auto: probe the input library for decoy transitions and automatically choose between the prepared and empirical paths
- @c prepared: accept an already prepared peptide query parameter library (`tsv`, `pqp`, `oswpq`, or `TraML`) and normalize it to an internal `prepared_library.pqp`
- @c empirical: accept an empirical transition library, run assay preparation
  and decoy generation, normalize the result to `prepared_library.pqp`, then continue

The downstream run then executes:
1. extraction/scoring to a working `workflow.osw`
2. in-process Percolator rescoring in the same working OSW
3. peptide/protein/gene/peptidoform inference in the same working OSW
4. final result and matrix export to `out_dir`

By default intermediate files are removed after success. Use
@p workflow:keep_intermediate_files to retain them, or
@p workflow:intermediate_dir to control their location.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenDIA.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenDIA.html
*/

/// @cond TOPPCLASSES
class TOPPOpenDIA final :
  public TOPPOpenSwathBase
{
public:
  TOPPOpenDIA() :
    TOPPOpenSwathBase("OpenDIA", "High-level DIA workflow front-end built on the OpenSWATH/OpenMS stack.", true)
  {
  }

protected:
  using Rescorer = OpenMS::OpenSwathPercolatorScoring;
  using RescoreLevel = Rescorer::Level;

  enum class LibraryMode
  {
    AUTO,
    PREPARED,
    EMPIRICAL
  };

  struct InferenceTask
  {
    InferenceLevel level = InferenceLevel::Peptidoform;
    std::optional<InferenceContext> context;
  };

  enum class ExportTaskType
  {
    FeatureParquet,
    Results,
    Matrix
  };

  struct ExportTask
  {
    ExportTaskType type = ExportTaskType::Results;
    std::optional<OpenSwathMatrixLevel> matrix_level;
  };

  struct WorkingDirectory
  {
    std::string path;
    std::unique_ptr<File::TempDir> temp_dir;
    bool remove_on_success = false;
  };

  static bool toBool_(const std::string& value)
  {
    return value == "true";
  }

  static StringList validRescoreLevels_()
  {
    StringList levels;
    levels.reserve(Rescorer::names_of_level.size());
    for (const auto& name : Rescorer::names_of_level)
    {
      levels.push_back(name);
    }
    return levels;
  }

  static RescoreLevel toRescoreLevel_(const std::string& value)
  {
    const auto it = std::find(Rescorer::names_of_level.begin(),
                              Rescorer::names_of_level.end(),
                              static_cast<std::string>(value));
    if (it == Rescorer::names_of_level.end())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Unsupported OpenSWATH Percolator level.",
                                    value);
    }
    return static_cast<RescoreLevel>(std::distance(Rescorer::names_of_level.begin(), it));
  }

  static std::vector<RescoreLevel> parseRescoreLevels_(const StringList& values)
  {
    std::vector<RescoreLevel> levels;
    for (const auto& value : values)
    {
      const RescoreLevel level = toRescoreLevel_(value);
      if (std::find(levels.begin(), levels.end(), level) == levels.end())
      {
        levels.push_back(level);
      }
    }
    return levels;
  }

  static OpenSwathExportFileFormat toExportFormat_(const std::string& value)
  {
    return value == "parquet" ? OpenSwathExportFileFormat::Parquet : OpenSwathExportFileFormat::TSV;
  }

  static OpenSwathIPFExportMode toIPFMode_(const std::string& value)
  {
    if (value == "augmented") return OpenSwathIPFExportMode::Augmented;
    if (value == "disable") return OpenSwathIPFExportMode::Disable;
    return OpenSwathIPFExportMode::Peptidoform;
  }

  static OpenSwathMatrixNormalization toNormalization_(const std::string& value)
  {
    if (value == "median") return OpenSwathMatrixNormalization::Median;
    if (value == "medianmedian") return OpenSwathMatrixNormalization::MedianMedian;
    return OpenSwathMatrixNormalization::None;
  }

  static void appendInferenceTask_(std::vector<InferenceTask>& tasks,
                                   const InferenceLevel level,
                                   const std::optional<InferenceContext>& context)
  {
    const auto duplicate = std::find_if(tasks.begin(), tasks.end(),
      [&](const auto& task)
      {
        return task.level == level && task.context == context;
      });
    if (duplicate == tasks.end())
    {
      tasks.push_back({level, context});
    }
  }

  static std::string inferenceTaskLabel_(const InferenceTask& task)
  {
    if (task.level == InferenceLevel::Peptidoform)
    {
      return "peptidoform inference";
    }
    return toString(task.level) + " inference in '" + toString(*task.context) + "' context";
  }

  static std::string targetLabelPlural_(const InferenceLevel level)
  {
    switch (level)
    {
      case InferenceLevel::Peptidoform: return "peptidoforms";
      case InferenceLevel::Peptide: return "peptides";
      case InferenceLevel::Protein: return "proteins";
      case InferenceLevel::Gene: return "genes";
    }
    return "targets";
  }

  static Int64 runKey_(const std::optional<Int64>& run_id)
  {
    return run_id.value_or(std::numeric_limits<Int64>::min());
  }

  void registerWorkflowOptions_()
  {
    registerInputFileList_("in", "<files>", StringList(), "Input DIA files separated by blank.");
    StringList in_formats = {"mzML", "mzXML", "sqMass"};
#ifdef WITH_OPENTIMS
    in_formats.push_back("d");
#endif
#ifdef WITH_THERMO_RAW
    in_formats.push_back("raw");
#endif
    setValidFormats_("in", in_formats);

    registerInputFile_("tr", "<file>", "", "Library input file.");
    StringList tr_formats = {"traML", "tsv", "pqp", "oswpq"};
    setValidFormats_("tr", tr_formats);
    registerStringOption_("tr_type", "<type>", "", "Library input file type -- default: determined from file extension or content.", false);
    setValidStrings_("tr_type", tr_formats);

    registerOutputDir_("out_dir", "<dir>", ".", "Directory for final exported OpenDIA outputs.", false, false);

    registerTOPPSubsection_("workflow", "Workflow options.");
    registerStringOption_("workflow:library_mode", "<choice>", "auto", "How to enter the workflow: auto-detect based on decoys already present in the input library, force prepared-library normalization, or force empirical assay/decoy preparation.", false);
    setValidStrings_("workflow:library_mode", {"auto", "prepared", "empirical"});
    registerStringOption_("workflow:keep_intermediate_files", "<true|false>", "false", "Whether to retain prepared_library.pqp and the single working workflow.osw after success.", false);
    setValidStrings_("workflow:keep_intermediate_files", {"true", "false"});
    registerOutputDir_("workflow:intermediate_dir", "<dir>", "", "Optional working directory for intermediate workflow files. Preserved automatically on failure.", false, false);
  }

  void registerTargetedDataExtractionOptions_()
  {
    registerTOPPSubsection_("TargetedDataExtraction", "Targeted data extraction options.");

    registerInputFile_("TargetedDataExtraction:swath_windows_file", "<file>", "", "Optional tab-separated SWATH window file used during extraction.", false);
    registerFlag_("TargetedDataExtraction:sort_swath_maps", "Sort input SWATH files when matching to SWATH windows from TargetedDataExtraction:swath_windows_file.", true);

    registerStringOption_("TargetedDataExtraction:enable_ms1", "<true|false>", "true", "Extract precursor traces and use them for scoring if present.", false, true);
    setValidStrings_("TargetedDataExtraction:enable_ms1", {"true", "false"});
    registerStringOption_("TargetedDataExtraction:enable_ipf", "<true|false>", "true", "Enable extraction-time UIS/IPF scoring when identifying transitions are present.", false, true);
    setValidStrings_("TargetedDataExtraction:enable_ipf", {"true", "false"});

    registerOutputFile_("TargetedDataExtraction:out_chrom", "<file>", "", "Optional chromatogram output in mzML (chrom.mzML), sqMass, or xic (Parquet).", false, true);
    setValidFormats_("TargetedDataExtraction:out_chrom", {"mzML", "sqMass", "xic"});
    registerOutputFile_("TargetedDataExtraction:out_mobilogram", "<file>", "", "Optional extracted ion mobilogram output in Parquet (.xim).", false, true);
    setValidFormats_("TargetedDataExtraction:out_mobilogram", {"xim"});

    registerDoubleOption_("TargetedDataExtraction:min_upper_edge_dist", "<double>", 0.0, "Minimal distance to the upper edge of a SWATH window to still consider a precursor, in Thomson.", false, true);
    registerFlag_("TargetedDataExtraction:pasef", "Treat the input as PASEF / diaPASEF data.");

    registerDoubleOption_("TargetedDataExtraction:rt_extraction_window", "<double>", 600.0, "Only extract RT around this value (-1 means extract the whole range; 600 means +/- 300 s around the expected elution).", false);
    registerDoubleOption_("TargetedDataExtraction:extra_rt_extraction_window", "<double>", 0.0, "Optional extra RT padding for the written chromatogram output.", false, true);
    setMinFloat_("TargetedDataExtraction:extra_rt_extraction_window", 0.0);
    registerDoubleOption_("TargetedDataExtraction:ion_mobility_window", "<double>", -1, "Extraction window in ion mobility dimension (full width; -1 extracts the whole range).", false);
    registerDoubleOption_("TargetedDataExtraction:mz_extraction_window", "<double>", 50, "MS2 extraction window in Thomson or ppm (see TargetedDataExtraction:mz_extraction_window_unit).", false);
    setMinFloat_("TargetedDataExtraction:mz_extraction_window", 0.0);
    registerStringOption_("TargetedDataExtraction:mz_extraction_window_unit", "<name>", "ppm", "Unit for the MS2 m/z extraction window.", false, true);
    setValidStrings_("TargetedDataExtraction:mz_extraction_window_unit", {"Th", "ppm"});

    registerDoubleOption_("TargetedDataExtraction:mz_extraction_window_ms1", "<double>", 50, "MS1 extraction window in Thomson or ppm (see TargetedDataExtraction:mz_extraction_window_ms1_unit).", false);
    setMinFloat_("TargetedDataExtraction:mz_extraction_window_ms1", 0.0);
    registerStringOption_("TargetedDataExtraction:mz_extraction_window_ms1_unit", "<name>", "ppm", "Unit for the MS1 m/z extraction window.", false, true);
    setValidStrings_("TargetedDataExtraction:mz_extraction_window_ms1_unit", {"ppm", "Th"});
    registerDoubleOption_("TargetedDataExtraction:im_extraction_window_ms1", "<double>", -1, "MS1 ion mobility extraction window. -1 disables MS1 ion mobility extraction.", false);
    registerStringOption_("TargetedDataExtraction:use_ms1_ion_mobility", "<name>", "true", "Also perform precursor extraction using the fragment-ion ion mobility window.", false, true);
    setValidStrings_("TargetedDataExtraction:use_ms1_ion_mobility", {"true", "false"});

    registerStringOption_("TargetedDataExtraction:matching_window_only", "<name>", "false", "Assume PRM-like overlapping DIA windows and extract each assay only from the best matching window.", false, true);
    setValidStrings_("TargetedDataExtraction:matching_window_only", {"true", "false"});

    registerDoubleOption_("TargetedDataExtraction:irt_mz_extraction_window", "<double>", 50, "Extraction window used for iRT and m/z correction in Thomson or ppm (see TargetedDataExtraction:irt_mz_extraction_window_unit).", false, true);
    setMinFloat_("TargetedDataExtraction:irt_mz_extraction_window", 0.0);
    registerStringOption_("TargetedDataExtraction:irt_mz_extraction_window_unit", "<name>", "ppm", "Unit for the iRT m/z extraction window.", false, true);
    setValidStrings_("TargetedDataExtraction:irt_mz_extraction_window_unit", {"Th", "ppm"});
    registerDoubleOption_("TargetedDataExtraction:irt_im_extraction_window", "<double>", -1, "Ion mobility extraction window used for iRT calibration (-1 disables IM calibration).", false, true);

    registerFlag_("TargetedDataExtraction:split_file_input", "Interpret the input files as one SWATH window per file.", true);
    registerFlag_("TargetedDataExtraction:use_elution_model_score", "Turn on the CPU-intensive EMG elution model score.", true);

    registerStringOption_("TargetedDataExtraction:readOptions", "<name>", "normal", "Whether to read directly, cache to disk first, or keep cached data in memory.", false, true);
    setValidStrings_("TargetedDataExtraction:readOptions", {"normal", "cache", "cacheWorkingInMemory", "workingInMemory"});
    registerStringOption_("TargetedDataExtraction:tempDirectory", "<tmp>", File::getTempDirectory(), "Temporary directory used for cached SWATH input files.", false, true);
    registerFlag_("TargetedDataExtraction:keep_cached_files", "If set, do not remove cached input files created in TargetedDataExtraction:tempDirectory.", false);

    registerStringOption_("TargetedDataExtraction:extraction_function", "<name>", "tophat", "Extraction kernel used for signal extraction.", false, true);
    setValidStrings_("TargetedDataExtraction:extraction_function", {"tophat", "bartlett"});

    registerIntOption_("TargetedDataExtraction:batchSize", "<number>", 0, "Compound batch size for extraction/scoring scheduling. 0 enables automatic scheduling.", false, true);
    setMinInt_("TargetedDataExtraction:batchSize", 0);
    registerIntOption_("TargetedDataExtraction:innerBatchSize", "<number>", -1, "Inner scoring batch size for automatic scheduling (<=0 means auto).", false, true);
    setMinInt_("TargetedDataExtraction:innerBatchSize", -1);
    registerIntOption_("TargetedDataExtraction:maxConcurrentSwaths", "<number>", -1, "Maximum concurrent non-MS1 SWATH maps kept in memory for automatic scheduling (-1 means auto).", false, true);
    setMinInt_("TargetedDataExtraction:maxConcurrentSwaths", -1);
    registerIntOption_("TargetedDataExtraction:outer_loop_threads", "<number>", -1, "Legacy nested OpenMP outer-loop thread count. Leave -1 for automatic scheduling.", false, true);
    registerIntOption_("TargetedDataExtraction:ms1_isotopes", "<number>", 3, "The number of MS1 isotopes used for extraction.", false, true);
    setMinInt_("TargetedDataExtraction:ms1_isotopes", 0);

    registerSubsection_("TargetedDataExtraction:Scoring", "Scoring parameters section.");
    registerSubsection_("TargetedDataExtraction:Library", "Transition-library reader parameters section.");
    registerSubsection_("TargetedDataExtraction:Calibration", "Parameters for calibrant iRT peptides for RT normalization and mass / ion mobility correction.");
    registerSubsection_("TargetedDataExtraction:Calibration:RTNormalization", "RT normalization model and outlier handling.");
    registerSubsection_("TargetedDataExtraction:Calibration:MassIMCorrection", "m/z and ion mobility calibration.");
    registerSubsection_("TargetedDataExtraction:MRMMapping", "Parameters for mapping chromatograms to transitions in SRM/MRM data.");

    registerOutputFile_("TargetedDataExtraction:irt_mzml", "<file>", "", "Optional chromatogram mzML containing the iRT peptides.", false);
    setValidFormats_("TargetedDataExtraction:irt_mzml", {"mzML"});
    registerOutputFile_("TargetedDataExtraction:irt_trafo", "<file>", "", "Optional transformation file for RT normalization.", false);
    setValidFormats_("TargetedDataExtraction:irt_trafo", {"trafoXML"});
    registerStringList_("TargetedDataExtraction:disable_features", "<list>", StringList(), "Debug-only feature toggles: no_IM_calibration, no_IM_windowing.", false, true);
    setValidStrings_("TargetedDataExtraction:disable_features", {"no_IM_calibration", "no_IM_windowing"});
  }

  void registerPeptideQueryParameterOptions_()
  {
    registerTOPPSubsection_("PeptideQueryParameters:AssayGenerator", "Assay-generation parameters used when workflow:library_mode=empirical.");
    registerIntOption_("PeptideQueryParameters:AssayGenerator:min_transitions", "<int>", 6, "Minimal number of transitions.", false);
    registerIntOption_("PeptideQueryParameters:AssayGenerator:max_transitions", "<int>", 6, "Maximal number of transitions.", false);
    registerStringOption_("PeptideQueryParameters:AssayGenerator:allowed_fragment_types", "<type>", "b,y", "Allowed fragment types.", false);
    registerStringOption_("PeptideQueryParameters:AssayGenerator:allowed_fragment_charges", "<type>", "1,2,3,4", "Allowed fragment charge states.", false);
    registerFlag_("PeptideQueryParameters:AssayGenerator:enable_detection_specific_losses", "Allow specific neutral losses for detecting transitions.");
    registerFlag_("PeptideQueryParameters:AssayGenerator:enable_detection_unspecific_losses", "Allow unspecific neutral losses for detecting transitions.");
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_mz_threshold", "<double>", 0.025, "Precursor m/z threshold in Thomson for precursor annotation.", false);
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_lower_mz_limit", "<double>", 400, "Lower precursor m/z limit.", false);
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_upper_mz_limit", "<double>", 1200, "Upper precursor m/z limit.", false);
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:product_mz_threshold", "<double>", 0.025, "Fragment m/z threshold in Thomson for annotation.", false);
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:product_lower_mz_limit", "<double>", 350, "Lower fragment m/z limit.", false);
    registerDoubleOption_("PeptideQueryParameters:AssayGenerator:product_upper_mz_limit", "<double>", 2000, "Upper fragment m/z limit.", false);
    registerInputFile_("PeptideQueryParameters:AssayGenerator:swath_windows_file", "<file>", "", "Optional SWATH window file for exclusion of fragment ions in the precursor isolation window.", false, false);
    setValidFormats_("PeptideQueryParameters:AssayGenerator:swath_windows_file", {"txt"});
    registerInputFile_("PeptideQueryParameters:AssayGenerator:unimod_file", "<file>", "", "Optional Unimod XML file used for IPF transition generation.", false, false);
    setValidFormats_("PeptideQueryParameters:AssayGenerator:unimod_file", {"xml"});
    registerFlag_("PeptideQueryParameters:AssayGenerator:enable_ipf", "Generate identifying transitions for IPF during empirical library preparation.");
    registerIntOption_("PeptideQueryParameters:AssayGenerator:max_num_alternative_localizations", "<int>", 10000, "IPF: maximum number of site-localization permutations.", false, true);
    registerFlag_("PeptideQueryParameters:AssayGenerator:disable_identification_ms2_precursors", "IPF: disable MS2 precursor ions for identification transitions.", true);
    registerFlag_("PeptideQueryParameters:AssayGenerator:disable_identification_specific_losses", "IPF: disable specific neutral losses for identification transitions.", true);
    registerFlag_("PeptideQueryParameters:AssayGenerator:enable_identification_unspecific_losses", "IPF: allow unspecific neutral losses for identification transitions.", true);
    registerFlag_("PeptideQueryParameters:AssayGenerator:enable_swath_specifity", "IPF: use precursor-isolation windows instead of precursor m/z-specific identification transitions.", true);
    registerFlag_("PeptideQueryParameters:AssayGenerator:disable_decoy_transitions", "IPF: disable generation of decoy UIS transitions.", true);
    registerIntOption_("PeptideQueryParameters:AssayGenerator:ipf_decoy_seed", "<int>", -1, "IPF: random seed for decoy shuffle (-1 = time-based).", false, true);

    registerTOPPSubsection_("PeptideQueryParameters:DecoyGenerator", "Decoy-generation parameters used when workflow:library_mode=empirical.");
    registerStringOption_("PeptideQueryParameters:DecoyGenerator:method", "<type>", "shuffle", "Decoy generation method.", false);
    setValidStrings_("PeptideQueryParameters:DecoyGenerator:method", {"shuffle", "pseudo-reverse", "reverse", "shift"});
    registerStringOption_("PeptideQueryParameters:DecoyGenerator:decoy_tag", "<type>", "DECOY_", "Decoy tag.", false);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:min_decoy_fraction", "<double>", 0.8, "Minimum fraction of decoy / target peptides and proteins.", false, true);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:aim_decoy_fraction", "<double>", 1.0, "Requested fraction of decoys to generate.", false, true);
    registerIntOption_("PeptideQueryParameters:DecoyGenerator:shuffle_max_attempts", "<int>", 30, "Shuffle: maximum attempts to reduce target/decoy sequence identity.", false, true);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:shuffle_sequence_identity_threshold", "<double>", 0.5, "Shuffle: target-decoy amino-acid sequence identity threshold.", false, true);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:shift_precursor_mz_shift", "<double>", 0.0, "Shift: precursor ion m/z shift in Thomson.", false, true);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:shift_product_mz_shift", "<double>", 20, "Shift: fragment ion m/z shift in Thomson.", false, true);
    registerDoubleOption_("PeptideQueryParameters:DecoyGenerator:product_mz_threshold", "<double>", 0.025, "Fragment m/z threshold used during decoy annotation.", false, true);
    registerStringOption_("PeptideQueryParameters:DecoyGenerator:allowed_fragment_types", "<type>", "b,y", "Allowed fragment types.", false, true);
    registerStringOption_("PeptideQueryParameters:DecoyGenerator:allowed_fragment_charges", "<type>", "1,2,3,4", "Allowed fragment charge states.", false, true);
    registerFlag_("PeptideQueryParameters:DecoyGenerator:enable_detection_specific_losses", "Allow specific neutral losses for decoy generation.", true);
    registerFlag_("PeptideQueryParameters:DecoyGenerator:enable_detection_unspecific_losses", "Allow unspecific neutral losses for decoy generation.", true);
    registerStringOption_("PeptideQueryParameters:DecoyGenerator:switchKR", "<true/false>", "true", "Whether to switch terminal K and R to achieve different precursor masses.", false);
    setValidStrings_("PeptideQueryParameters:DecoyGenerator:switchKR", {"true", "false"});
    registerFlag_("PeptideQueryParameters:DecoyGenerator:separate", "Write only decoys instead of appending them to the targets.", true);
  }

  void registerRescoringOptions_()
  {
    registerTOPPSubsection_("Rescoring", "In-process Percolator rescoring parameters.");
    registerStringList_("Rescoring:level", "<levels>", StringList{"ms1ms2"},
                        "One or more rescoring levels to run sequentially. Transition rescoring is appended automatically when peptidoform inference is enabled.",
                        false);
    setValidStrings_("Rescoring:level", validRescoreLevels_());
  }

  void registerInferenceOptions_()
  {
    registerTOPPSubsection_("Inference", "Peptide/protein/gene/peptidoform inference parameters.");
    registerTOPPSubsection_("Inference:peptidoform", "Peptidoform inference options.");
    registerTOPPSubsection_("Inference:peptide", "Peptide inference options.");
    registerTOPPSubsection_("Inference:protein", "Protein inference options.");
    registerTOPPSubsection_("Inference:gene", "Gene inference options.");

    registerStringOption_("Inference:peptidoform:run", "<true|false>", "false", "Enable peptidoform inference.", false);
    setValidStrings_("Inference:peptidoform:run", {"true", "false"});
    const PeptidoformInferenceConfig ipf_defaults;
    registerStringOption_("Inference:peptidoform:ipf_ms1_scoring", "<true|false>", ipf_defaults.ipf_ms1_scoring ? "true" : "false", "Use MS1 precursor evidence for peptidoform inference.", false);
    setValidStrings_("Inference:peptidoform:ipf_ms1_scoring", {"true", "false"});
    registerStringOption_("Inference:peptidoform:ipf_ms2_scoring", "<true|false>", ipf_defaults.ipf_ms2_scoring ? "true" : "false", "Use MS2 precursor evidence for peptidoform inference.", false);
    setValidStrings_("Inference:peptidoform:ipf_ms2_scoring", {"true", "false"});
    registerStringOption_("Inference:peptidoform:ipf_h0", "<true|false>", ipf_defaults.ipf_h0 ? "true" : "false", "Include the null hypothesis in IPF inference.", false, true);
    setValidStrings_("Inference:peptidoform:ipf_h0", {"true", "false"});
    registerStringOption_("Inference:peptidoform:ipf_grouped_fdr", "<true|false>", ipf_defaults.ipf_grouped_fdr ? "true" : "false", "Compute model FDR separately per num_peptidoforms group.", false, true);
    setValidStrings_("Inference:peptidoform:ipf_grouped_fdr", {"true", "false"});
    registerDoubleOption_("Inference:peptidoform:ipf_max_peakgroup_pep", "<num>", ipf_defaults.ipf_max_peakgroup_pep, "Maximum peakgroup PEP for entering precursor inference.", false, true);
    registerDoubleOption_("Inference:peptidoform:ipf_max_precursor_pep", "<num>", ipf_defaults.ipf_max_precursor_pep, "Maximum precursor-level evidence PEP.", false, true);
    registerDoubleOption_("Inference:peptidoform:ipf_max_precursor_peakgroup_pep", "<num>", ipf_defaults.ipf_max_precursor_peakgroup_pep, "Maximum precursor-peakgroup PEP retained for transition inference.", false, true);
    registerDoubleOption_("Inference:peptidoform:ipf_max_transition_pep", "<num>", ipf_defaults.ipf_max_transition_pep, "Maximum transition evidence PEP retained for IPF.", false, true);
    registerStringOption_("Inference:peptidoform:propagate_signal_across_runs", "<true|false>", ipf_defaults.propagate_signal_across_runs ? "true" : "false", "Propagate confident transition evidence across aligned runs.", false);
    setValidStrings_("Inference:peptidoform:propagate_signal_across_runs", {"true", "false"});
    registerDoubleOption_("Inference:peptidoform:ipf_min_alignment_mapping_confidence", "<num>", ipf_defaults.ipf_min_alignment_mapping_confidence, "Minimum mapping confidence retained for across-run propagation.", false, true);
    registerDoubleOption_("Inference:peptidoform:across_run_confidence_threshold", "<num>", ipf_defaults.across_run_confidence_threshold, "Maximum transition PEP eligible for across-run propagation.", false, true);

    registerStringOption_("Inference:peptide:global", "<true|false>", "true", "Enable peptide inference in the global context.", false);
    setValidStrings_("Inference:peptide:global", {"true", "false"});
    registerStringOption_("Inference:peptide:experiment_wide", "<true|false>", "false", "Enable peptide inference in the experiment-wide context.", false);
    setValidStrings_("Inference:peptide:experiment_wide", {"true", "false"});
    registerStringOption_("Inference:peptide:run_specific", "<true|false>", "false", "Enable peptide inference in the run-specific context.", false);
    setValidStrings_("Inference:peptide:run_specific", {"true", "false"});

    registerStringOption_("Inference:protein:global", "<true|false>", "true", "Enable protein inference in the global context.", false);
    setValidStrings_("Inference:protein:global", {"true", "false"});
    registerStringOption_("Inference:protein:experiment_wide", "<true|false>", "false", "Enable protein inference in the experiment-wide context.", false);
    setValidStrings_("Inference:protein:experiment_wide", {"true", "false"});
    registerStringOption_("Inference:protein:run_specific", "<true|false>", "false", "Enable protein inference in the run-specific context.", false);
    setValidStrings_("Inference:protein:run_specific", {"true", "false"});

    registerStringOption_("Inference:gene:global", "<true|false>", "false", "Enable gene inference in the global context.", false);
    setValidStrings_("Inference:gene:global", {"true", "false"});
    registerStringOption_("Inference:gene:experiment_wide", "<true|false>", "false", "Enable gene inference in the experiment-wide context.", false);
    setValidStrings_("Inference:gene:experiment_wide", {"true", "false"});
    registerStringOption_("Inference:gene:run_specific", "<true|false>", "false", "Enable gene inference in the run-specific context.", false);
    setValidStrings_("Inference:gene:run_specific", {"true", "false"});

    registerTOPPSubsection_("Inference:error_estimation", "Context-level error estimation options.");
    const ErrorEstimationConfig error_defaults;
    registerStringOption_("Inference:error_estimation:parametric", "<true|false>", error_defaults.parametric ? "true" : "false", "Use parametric p-value estimation instead of empirical estimation.", false, true);
    setValidStrings_("Inference:error_estimation:parametric", {"true", "false"});
    registerStringOption_("Inference:error_estimation:pfdr", "<true|false>", error_defaults.pfdr ? "true" : "false", "Use positive FDR in q-value and s-value calculations.", false, true);
    setValidStrings_("Inference:error_estimation:pfdr", {"true", "false"});
    registerDoubleList_("Inference:error_estimation:pi0_lambda", "<l1,l2,...>", error_defaults.pi0_lambda, "Lambda grid for pi0 estimation.", false, true);
    registerStringOption_("Inference:error_estimation:pi0_method", "<choice>", error_defaults.pi0_method, "Method for pi0 estimation.", false, true);
    setValidStrings_("Inference:error_estimation:pi0_method", {"bootstrap", "smoother"});
    registerIntOption_("Inference:error_estimation:pi0_smooth_df", "<num>", error_defaults.pi0_smooth_df, "Degrees of freedom for smoother-based pi0 estimation.", false, true);
    setMinInt_("Inference:error_estimation:pi0_smooth_df", 1);
    registerStringOption_("Inference:error_estimation:pi0_smooth_log_pi0", "<true|false>", error_defaults.pi0_smooth_log_pi0 ? "true" : "false", "Smooth log(pi0) instead of pi0 directly.", false, true);
    setValidStrings_("Inference:error_estimation:pi0_smooth_log_pi0", {"true", "false"});
    registerStringOption_("Inference:error_estimation:lfdr_truncate", "<true|false>", error_defaults.lfdr_truncate ? "true" : "false", "Truncate local FDR/PEP values to [0,1].", false, true);
    setValidStrings_("Inference:error_estimation:lfdr_truncate", {"true", "false"});
    registerStringOption_("Inference:error_estimation:lfdr_monotone", "<true|false>", error_defaults.lfdr_monotone ? "true" : "false", "Enforce monotone local FDR/PEP estimates.", false, true);
    setValidStrings_("Inference:error_estimation:lfdr_monotone", {"true", "false"});
    registerStringOption_("Inference:error_estimation:lfdr_transformation", "<choice>", error_defaults.lfdr_transformation, "Transformation used by local FDR estimation.", false, true);
    setValidStrings_("Inference:error_estimation:lfdr_transformation", {"probit", "logit"});
    registerDoubleOption_("Inference:error_estimation:lfdr_adj", "<num>", error_defaults.lfdr_adj, "Adjustment factor for local FDR estimation bandwidth.", false, true);
    registerDoubleOption_("Inference:error_estimation:lfdr_eps", "<num>", error_defaults.lfdr_eps, "Epsilon used to clip p-values before local FDR estimation.", false, true);
  }

  void registerExportOptions_()
  {
    registerTOPPSubsection_("Export", "User-facing result and matrix export parameters.");
    registerTOPPSubsection_("Export:parquet", "Options for scored feature parquet export.");
    registerStringOption_("Export:parquet:run", "<true|false>", "false", "Whether to export scored feature parquet tables.", false);
    setValidStrings_("Export:parquet:run", {"true", "false"});
    registerStringOption_("Export:parquet:include_transition_data", "<true|false>", "true", "Whether to also export transition-level parquet tables.", false);
    setValidStrings_("Export:parquet:include_transition_data", {"true", "false"});
    registerStringOption_("Export:parquet:exclude_decoys", "<true|false>", "false", "Whether to exclude decoys from parquet feature exports.", false);
    setValidStrings_("Export:parquet:exclude_decoys", {"true", "false"});

    registerTOPPSubsection_("Export:results", "Options for user-facing filtered result export.");
    registerStringOption_("Export:results:run", "<true|false>", "true", "Whether to export the filtered results table.", false);
    setValidStrings_("Export:results:run", {"true", "false"});
    registerStringOption_("Export:results:format", "<choice>", "tsv", "Output format for the filtered results table.", false);
    setValidStrings_("Export:results:format", {"tsv", "parquet"});
    registerStringOption_("Export:results:ipf", "<choice>", "peptidoform", "How IPF results should be represented when SCORE_IPF is present.", false);
    setValidStrings_("Export:results:ipf", {"peptidoform", "augmented", "disable"});
    registerStringOption_("Export:results:transition_quantification", "<true|false>", "true", "Whether to aggregate transition-level quantification columns.", false);
    setValidStrings_("Export:results:transition_quantification", {"true", "false"});
    registerDoubleOption_("Export:results:max_transition_pep", "<num>", 0.7, "Maximum transition PEP used for aggregated transition quantification.", false);
    registerDoubleOption_("Export:results:ipf_max_peptidoform_pep", "<num>", 0.4, "Maximum IPF PEP used for peptidoform export filtering.", false);
    registerDoubleOption_("Export:results:max_rs_peakgroup_qvalue", "<num>", 0.05, "Maximum run-specific peakgroup q-value retained.", false);
    registerStringOption_("Export:results:use_peptide_scores", "<true|false>", "true", "Whether to append peptide-level scores if SCORE_PEPTIDE is present and to require passing the peptide-level threshold.", false);
    setValidStrings_("Export:results:use_peptide_scores", {"true", "false"});
    registerDoubleOption_("Export:results:max_global_peptide_qvalue", "<num>", 0.01, "Maximum global peptide-level q-value retained when peptide filtering is enabled.", false);
    registerStringOption_("Export:results:use_protein_scores", "<true|false>", "true", "Whether to append protein-level scores if SCORE_PROTEIN is present and to require passing the protein-level threshold.", false);
    setValidStrings_("Export:results:use_protein_scores", {"true", "false"});
    registerDoubleOption_("Export:results:max_global_protein_qvalue", "<num>", 0.01, "Maximum global protein-level q-value retained when protein filtering is enabled.", false);
    registerStringOption_("Export:results:use_gene_scores", "<true|false>", "false", "Whether to append gene-level scores if SCORE_GENE is present and to require passing the gene-level threshold.", false);
    setValidStrings_("Export:results:use_gene_scores", {"true", "false"});
    registerDoubleOption_("Export:results:max_global_gene_qvalue", "<num>", 0.01, "Maximum global gene-level q-value retained when gene filtering is enabled.", false);
    registerStringOption_("Export:results:use_alignment", "<true|false>", "false", "Whether to recover aligned features when legacy alignment tables are present.", false);
    setValidStrings_("Export:results:use_alignment", {"true", "false"});
    registerDoubleOption_("Export:results:max_alignment_pep", "<num>", 0.7, "Maximum alignment PEP retained when alignment recovery is enabled.", false);
    registerStringOption_("Export:results:exclude_decoys", "<true|false>", "true", "Whether to exclude decoys from the filtered results export.", false);
    setValidStrings_("Export:results:exclude_decoys", {"true", "false"});

    registerTOPPSubsection_("Export:matrix", "Options for quantification matrix exports.");
    registerStringOption_("Export:matrix:precursor", "<true|false>", "false", "Whether to export the precursor-level quantification matrix.", false);
    setValidStrings_("Export:matrix:precursor", {"true", "false"});
    registerStringOption_("Export:matrix:peptide", "<true|false>", "true", "Whether to export the peptide-level quantification matrix.", false);
    setValidStrings_("Export:matrix:peptide", {"true", "false"});
    registerStringOption_("Export:matrix:protein", "<true|false>", "true", "Whether to export the protein-level quantification matrix.", false);
    setValidStrings_("Export:matrix:protein", {"true", "false"});
    registerStringOption_("Export:matrix:gene", "<true|false>", "false", "Whether to export the gene-level quantification matrix.", false);
    setValidStrings_("Export:matrix:gene", {"true", "false"});
    registerStringOption_("Export:matrix:format", "<choice>", "tsv", "Output format for quantification matrices.", false);
    setValidStrings_("Export:matrix:format", {"tsv", "parquet"});
    registerStringOption_("Export:matrix:ipf", "<choice>", "peptidoform", "How IPF results should be represented when SCORE_IPF is present.", false);
    setValidStrings_("Export:matrix:ipf", {"peptidoform", "augmented", "disable"});
    registerStringOption_("Export:matrix:transition_quantification", "<true|false>", "true", "Whether to aggregate transition-level quantification columns onto the filtered feature rows.", false);
    setValidStrings_("Export:matrix:transition_quantification", {"true", "false"});
    registerDoubleOption_("Export:matrix:max_transition_pep", "<num>", 0.7, "Maximum transition PEP used for aggregated transition quantification.", false);
    registerDoubleOption_("Export:matrix:ipf_max_peptidoform_pep", "<num>", 0.4, "Maximum IPF PEP used for peptidoform export filtering.", false);
    registerDoubleOption_("Export:matrix:max_rs_peakgroup_qvalue", "<num>", 0.05, "Maximum run-specific peakgroup q-value retained.", false);
    registerStringOption_("Export:matrix:use_peptide_scores", "<true|false>", "true", "Whether peptide-level scores should be required when peptide filtering is enabled.", false);
    setValidStrings_("Export:matrix:use_peptide_scores", {"true", "false"});
    registerDoubleOption_("Export:matrix:max_global_peptide_qvalue", "<num>", 0.01, "Maximum global peptide-level q-value retained when peptide filtering is enabled.", false);
    registerStringOption_("Export:matrix:use_protein_scores", "<true|false>", "true", "Whether protein-level scores should be required when protein filtering is enabled.", false);
    setValidStrings_("Export:matrix:use_protein_scores", {"true", "false"});
    registerDoubleOption_("Export:matrix:max_global_protein_qvalue", "<num>", 0.01, "Maximum global protein-level q-value retained when protein filtering is enabled.", false);
    registerStringOption_("Export:matrix:use_gene_scores", "<true|false>", "false", "Whether gene-level scores should be required when gene filtering is enabled.", false);
    setValidStrings_("Export:matrix:use_gene_scores", {"true", "false"});
    registerDoubleOption_("Export:matrix:max_global_gene_qvalue", "<num>", 0.01, "Maximum global gene-level q-value retained when gene filtering is enabled.", false);
    registerStringOption_("Export:matrix:use_alignment", "<true|false>", "false", "Whether to recover aligned features when legacy alignment tables are present.", false);
    setValidStrings_("Export:matrix:use_alignment", {"true", "false"});
    registerDoubleOption_("Export:matrix:max_alignment_pep", "<num>", 0.7, "Maximum alignment PEP retained when alignment recovery is enabled.", false);
    registerStringOption_("Export:matrix:exclude_decoys", "<true|false>", "true", "Whether to exclude decoys from the filtered matrix input rows.", false);
    setValidStrings_("Export:matrix:exclude_decoys", {"true", "false"});
    registerIntOption_("Export:matrix:top_n", "<num>", 3, "Number of top precursors / peptides retained during matrix summarization.", false);
    registerStringOption_("Export:matrix:consistent_top", "<true|false>", "true", "Whether to use the same top precursors / peptides across all runs.", false);
    setValidStrings_("Export:matrix:consistent_top", {"true", "false"});
    registerStringOption_("Export:matrix:normalization", "<choice>", "none", "Normalization applied to matrix sample columns after summarization.", false);
    setValidStrings_("Export:matrix:normalization", {"none", "median", "medianmedian"});
  }

  void registerOptionsAndFlags_() override
  {
    registerWorkflowOptions_();
    registerPeptideQueryParameterOptions_();
    registerTargetedDataExtractionOptions_();
    registerRescoringOptions_();
    registerInferenceOptions_();
    registerExportOptions_();
  }

  Param getSubsectionDefaults_(const std::string& name) const override
  {
    if (name == "TargetedDataExtraction:Scoring")
    {
      Param feature_finder_param = MRMFeatureFinderScoring().getDefaults();
      feature_finder_param.remove("rt_extraction_window");
      feature_finder_param.setValue("stop_report_after_feature", 5);
      feature_finder_param.setValue("rt_normalization_factor", 100.0);
      feature_finder_param.setValue("Scores:use_ms1_mi", "true");
      feature_finder_param.setValue("Scores:use_mi_score", "true");
      feature_finder_param.setValue("TransitionGroupPicker:min_peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_quality", "false");
      feature_finder_param.setValue("TransitionGroupPicker:minimal_quality", -1.5);
      feature_finder_param.setValue("TransitionGroupPicker:background_subtraction", "none");
      feature_finder_param.setValue("TransitionGroupPicker:compute_peak_shape_metrics", "false");
      feature_finder_param.remove("TransitionGroupPicker:stop_after_intensity_ratio");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:use_gauss", "false");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_polynomial_order", 3);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:sgolay_frame_length", 11);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:peak_width", -1.0);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:remove_overlapping_peaks", "true");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:write_sn_log_messages", "false");
      feature_finder_param.setValue("TransitionGroupPicker:recalculate_peaks_max_z", 0.75);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:method", "corrected");
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:signal_to_noise", 0.1);
      feature_finder_param.setValue("TransitionGroupPicker:PeakPickerChromatogram:gauss_width", 30.0);
      feature_finder_param.setValue("uis_threshold_sn", -1);
      feature_finder_param.setValue("uis_threshold_peak_area", 0);
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_win_len");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:sn_bin_count");
      feature_finder_param.remove("TransitionGroupPicker:PeakPickerChromatogram:stop_after_feature");
      feature_finder_param.setValue("Scores:use_ion_mobility_scores", "auto");
      feature_finder_param.setValidStrings("Scores:use_ion_mobility_scores", {"true", "false", "auto"});
      feature_finder_param.remove("Scores:use_elution_model_score");
      feature_finder_param.setValue("EMGScoring:max_iteration", 10);
      feature_finder_param.remove("EMGScoring:interpolation_step");
      feature_finder_param.remove("EMGScoring:tolerance_stdev_bounding_box");
      feature_finder_param.remove("EMGScoring:deltaAbsError");
      feature_finder_param.remove("EMGScoring:statistics:mean");
      feature_finder_param.remove("EMGScoring:statistics:variance");
      return feature_finder_param;
    }
    else if (name == "TargetedDataExtraction:Library")
    {
      return TransitionTSVFile().getDefaults();
    }
    else if (name == "TargetedDataExtraction:Calibration")
    {
      CalibrationWorkflow cal_wf;
      Param p = cal_wf.getDefaults();
      p.setValue("tr_irt_priority_sampling", "", "Optional custom TSV transition file containing additional priority peptides for auto-iRT sampling.");
      p.setValue("rt_norm", "", "RT normalization file (if set, linear iRT files may be omitted).");
      return p;
    }
    else if (name == "TargetedDataExtraction:Calibration:RTNormalization")
    {
      Param p;
      p.setValue("alignmentMethod", "linear", "How to perform the RT alignment.");
      p.setValidStrings("alignmentMethod", {"linear", "interpolated", "lowess", "b_spline"});
      p.setValue("lowess:auto_span", "true", "Automatically select LOWESS span by cross-validation.");
      p.setValidStrings("lowess:auto_span", {"true", "false"});
      p.setValue("lowess:span", 0.05, "Span parameter for lowess.");
      p.setMinFloat("lowess:span", 0.0);
      p.setMaxFloat("lowess:span", 1.0);
      p.setValue("lowess:auto_span_min", 0.15, "Lower bound for auto-selected LOWESS span.");
      p.setMinFloat("lowess:auto_span_min", 0.001);
      p.setValue("lowess:auto_span_max", 0.80, "Upper bound for auto-selected LOWESS span.");
      p.setMaxFloat("lowess:auto_span_max", 0.99);
      p.setValue("lowess:auto_span_grid", "0.005,0.01,0.05,0.15,0.25,0.30,0.50,0.70,0.90", "Optional explicit grid of LOWESS span candidates.");
      p.setValue("b_spline:num_nodes", 5, "Number of nodes for b spline.");
      p.setMinInt("b_spline:num_nodes", 0);
      p.setValue("useIterativeChauvenet", "false", "Use Chauvenet's criterion in iterative outlier removal.");
      p.setValidStrings("useIterativeChauvenet", {"true", "false"});
      p.setValue("RANSACMaxIterations", 1000, "Maximum iterations for RANSAC outlier detection.");
      p.setValue("RANSACMaxPercentRTThreshold", 3, "Maximum RT threshold in percent of the total gradient for RANSAC.");
      p.setValue("RANSACSamplingSize", 10, "Sampling size per RANSAC iteration.");
      p.setValue("estimateBestPeptides", "false", "Whether to choose the best peptides based on peak shape for RT normalization.");
      p.setValidStrings("estimateBestPeptides", {"true", "false"});
      p.setValue("InitialQualityCutoff", 0.5, "Initial overall quality cutoff for a peak to be scored.");
      p.setValue("OverallQualityCutoff", 5.5, "Overall quality cutoff for a peak to enter RT estimation.");
      p.setValue("NrRTBins", 10, "Number of RT bins to use to compute coverage.");
      p.setValue("MinPeptidesPerBin", 1, "Minimal number of peptides required for a bin to be counted as covered.");
      p.setValue("MinBinsFilled", 8, "Minimal number of bins required to be covered.");
      return p;
    }
    else if (name == "TargetedDataExtraction:MRMMapping")
    {
      Param p;
      p.setValue("precursor_tolerance", 0.9, "Precursor tolerance when mapping (in Th).");
      p.setValue("product_tolerance", 1.2, "Product tolerance when mapping (in Th).");
      p.setValue("irt_precursor_tolerance", 1.5, "Precursor tolerance when mapping iRT transitions (in Th).");
      p.setValue("irt_product_tolerance", 1.5, "Product tolerance when mapping iRT transitions (in Th).");
      p.setValue("map_multiple_assays", "false", "Allow mapping multiple assays to chromatograms and duplicate the chromatograms in the output.");
      p.setValidStrings("map_multiple_assays", {"true", "false"});
      p.setValue("error_on_unmapped", "false", "Treat remaining unmapped chromatograms as an error.");
      p.setValidStrings("error_on_unmapped", {"true", "false"});
      return p;
    }
    else if (name == "TargetedDataExtraction:Calibration:MassIMCorrection")
    {
      return SwathMapMassCorrection().getDefaults();
    }
    else if (name == "Rescoring")
    {
      return rescoring_defaults_.getDefaults();
    }
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown subsection", name);
  }

  LibraryMode getLibraryMode_() const
  {
    const std::string mode = getStringOption_("workflow:library_mode");
    if (mode == "empirical")
    {
      return LibraryMode::EMPIRICAL;
    }
    if (mode == "prepared")
    {
      return LibraryMode::PREPARED;
    }
    return LibraryMode::AUTO;
  }

  OpenSwathLibraryPreparation::AssayGeneratorParameters getAssayGeneratorParameters_() const
  {
    OpenSwathLibraryPreparation::AssayGeneratorParameters parameters;
    parameters.min_transitions = getIntOption_("PeptideQueryParameters:AssayGenerator:min_transitions");
    parameters.max_transitions = getIntOption_("PeptideQueryParameters:AssayGenerator:max_transitions");
    StringUtils::split(getStringOption_("PeptideQueryParameters:AssayGenerator:allowed_fragment_types"), ",", parameters.allowed_fragment_types);
    std::vector<std::string> charge_strings;
    StringUtils::split(getStringOption_("PeptideQueryParameters:AssayGenerator:allowed_fragment_charges"), ",", charge_strings);
    parameters.allowed_fragment_charges.clear();
    parameters.allowed_fragment_charges.reserve(charge_strings.size());
    for (const auto& charge_string : charge_strings)
    {
      parameters.allowed_fragment_charges.push_back(static_cast<size_t>(std::atoi(charge_string.c_str())));
    }
    parameters.enable_detection_specific_losses = getFlag_("PeptideQueryParameters:AssayGenerator:enable_detection_specific_losses");
    parameters.enable_detection_unspecific_losses = getFlag_("PeptideQueryParameters:AssayGenerator:enable_detection_unspecific_losses");
    parameters.precursor_mz_threshold = getDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_mz_threshold");
    parameters.precursor_lower_mz_limit = getDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_lower_mz_limit");
    parameters.precursor_upper_mz_limit = getDoubleOption_("PeptideQueryParameters:AssayGenerator:precursor_upper_mz_limit");
    parameters.product_mz_threshold = getDoubleOption_("PeptideQueryParameters:AssayGenerator:product_mz_threshold");
    parameters.product_lower_mz_limit = getDoubleOption_("PeptideQueryParameters:AssayGenerator:product_lower_mz_limit");
    parameters.product_upper_mz_limit = getDoubleOption_("PeptideQueryParameters:AssayGenerator:product_upper_mz_limit");
    parameters.enable_ipf = getFlag_("PeptideQueryParameters:AssayGenerator:enable_ipf");
    parameters.max_num_alternative_localizations = getIntOption_("PeptideQueryParameters:AssayGenerator:max_num_alternative_localizations");
    parameters.enable_identification_ms2_precursors = !getFlag_("PeptideQueryParameters:AssayGenerator:disable_identification_ms2_precursors");
    parameters.enable_identification_specific_losses = !getFlag_("PeptideQueryParameters:AssayGenerator:disable_identification_specific_losses");
    parameters.enable_identification_unspecific_losses = getFlag_("PeptideQueryParameters:AssayGenerator:enable_identification_unspecific_losses");
    parameters.enable_swath_specifity = getFlag_("PeptideQueryParameters:AssayGenerator:enable_swath_specifity");
    parameters.disable_decoy_transitions = getFlag_("PeptideQueryParameters:AssayGenerator:disable_decoy_transitions");
    parameters.ipf_decoy_seed = getIntOption_("PeptideQueryParameters:AssayGenerator:ipf_decoy_seed");
    parameters.test_mode = getFlag_("test");
    parameters.unimod_file = getStringOption_("PeptideQueryParameters:AssayGenerator:unimod_file");

    const std::string swath_windows_file = getStringOption_("PeptideQueryParameters:AssayGenerator:swath_windows_file");
    if (!swath_windows_file.empty())
    {
      std::vector<double> lower;
      std::vector<double> upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, lower, upper);
      for (Size i = 0; i < lower.size(); ++i)
      {
        parameters.swathes.emplace_back(lower[i], upper[i]);
      }
    }
    return parameters;
  }

  OpenSwathLibraryPreparation::DecoyGeneratorParameters getDecoyGeneratorParameters_() const
  {
    OpenSwathLibraryPreparation::DecoyGeneratorParameters parameters;
    parameters.method = getStringOption_("PeptideQueryParameters:DecoyGenerator:method");
    parameters.decoy_tag = getStringOption_("PeptideQueryParameters:DecoyGenerator:decoy_tag");
    parameters.min_decoy_fraction = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:min_decoy_fraction");
    parameters.aim_decoy_fraction = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:aim_decoy_fraction");
    parameters.shuffle_max_attempts = getIntOption_("PeptideQueryParameters:DecoyGenerator:shuffle_max_attempts");
    parameters.shuffle_sequence_identity_threshold = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:shuffle_sequence_identity_threshold");
    parameters.shift_precursor_mz_shift = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:shift_precursor_mz_shift");
    parameters.shift_product_mz_shift = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:shift_product_mz_shift");
    parameters.product_mz_threshold = getDoubleOption_("PeptideQueryParameters:DecoyGenerator:product_mz_threshold");
    StringUtils::split(getStringOption_("PeptideQueryParameters:DecoyGenerator:allowed_fragment_types"), ",", parameters.allowed_fragment_types);
    std::vector<std::string> charge_strings;
    StringUtils::split(getStringOption_("PeptideQueryParameters:DecoyGenerator:allowed_fragment_charges"), ",", charge_strings);
    parameters.allowed_fragment_charges.clear();
    parameters.allowed_fragment_charges.reserve(charge_strings.size());
    for (const auto& charge_string : charge_strings)
    {
      parameters.allowed_fragment_charges.push_back(static_cast<size_t>(std::atoi(charge_string.c_str())));
    }
    parameters.enable_detection_specific_losses = getFlag_("PeptideQueryParameters:DecoyGenerator:enable_detection_specific_losses");
    parameters.enable_detection_unspecific_losses = getFlag_("PeptideQueryParameters:DecoyGenerator:enable_detection_unspecific_losses");
    parameters.switch_kr = getStringOption_("PeptideQueryParameters:DecoyGenerator:switchKR") == "true";
    parameters.separate = getFlag_("PeptideQueryParameters:DecoyGenerator:separate");
    return parameters;
  }

  PeptidoformInferenceConfig getPeptidoformConfig_() const
  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = toBool_(getStringOption_("Inference:peptidoform:ipf_ms1_scoring"));
    config.ipf_ms2_scoring = toBool_(getStringOption_("Inference:peptidoform:ipf_ms2_scoring"));
    config.ipf_h0 = toBool_(getStringOption_("Inference:peptidoform:ipf_h0"));
    config.ipf_grouped_fdr = toBool_(getStringOption_("Inference:peptidoform:ipf_grouped_fdr"));
    config.ipf_max_peakgroup_pep = getDoubleOption_("Inference:peptidoform:ipf_max_peakgroup_pep");
    config.ipf_max_precursor_pep = getDoubleOption_("Inference:peptidoform:ipf_max_precursor_pep");
    config.ipf_max_precursor_peakgroup_pep = getDoubleOption_("Inference:peptidoform:ipf_max_precursor_peakgroup_pep");
    config.ipf_max_transition_pep = getDoubleOption_("Inference:peptidoform:ipf_max_transition_pep");
    config.propagate_signal_across_runs = toBool_(getStringOption_("Inference:peptidoform:propagate_signal_across_runs"));
    config.ipf_min_alignment_mapping_confidence = getDoubleOption_("Inference:peptidoform:ipf_min_alignment_mapping_confidence");
    config.across_run_confidence_threshold = getDoubleOption_("Inference:peptidoform:across_run_confidence_threshold");
    return config;
  }

  ErrorEstimationConfig getErrorConfig_() const
  {
    ErrorEstimationConfig config;
    config.parametric = toBool_(getStringOption_("Inference:error_estimation:parametric"));
    config.pfdr = toBool_(getStringOption_("Inference:error_estimation:pfdr"));
    config.pi0_lambda = getDoubleList_("Inference:error_estimation:pi0_lambda");
    config.pi0_method = getStringOption_("Inference:error_estimation:pi0_method");
    config.pi0_smooth_df = getIntOption_("Inference:error_estimation:pi0_smooth_df");
    config.pi0_smooth_log_pi0 = toBool_(getStringOption_("Inference:error_estimation:pi0_smooth_log_pi0"));
    config.lfdr_truncate = toBool_(getStringOption_("Inference:error_estimation:lfdr_truncate"));
    config.lfdr_monotone = toBool_(getStringOption_("Inference:error_estimation:lfdr_monotone"));
    config.lfdr_transformation = getStringOption_("Inference:error_estimation:lfdr_transformation");
    config.lfdr_adj = getDoubleOption_("Inference:error_estimation:lfdr_adj");
    config.lfdr_eps = getDoubleOption_("Inference:error_estimation:lfdr_eps");
    return config;
  }

  std::vector<InferenceTask> getInferenceTasks_() const
  {
    std::vector<InferenceTask> tasks;
    if (toBool_(getStringOption_("Inference:peptidoform:run")))
    {
      appendInferenceTask_(tasks, InferenceLevel::Peptidoform, std::nullopt);
    }

    const std::array<std::pair<std::string, InferenceLevel>, 3> level_prefixes =
    {{
      {"peptide", InferenceLevel::Peptide},
      {"protein", InferenceLevel::Protein},
      {"gene", InferenceLevel::Gene}
    }};
    const std::array<std::pair<std::string, InferenceContext>, 3> contexts =
    {{
      {"global", InferenceContext::Global},
      {"experiment_wide", InferenceContext::ExperimentWide},
      {"run_specific", InferenceContext::RunSpecific}
    }};

    for (const auto& level_prefix : level_prefixes)
    {
      for (const auto& context : contexts)
      {
        if (toBool_(getStringOption_("Inference:" + level_prefix.first + ":" + context.first)))
        {
          appendInferenceTask_(tasks, level_prefix.second, context.second);
        }
      }
    }
    return tasks;
  }

  OpenSwathExportFilterConfig getFilterConfig_(const std::string& prefix) const
  {
    OpenSwathExportFilterConfig config;
    config.ipf_mode = toIPFMode_(getStringOption_(prefix + ":ipf"));
    config.transition_quantification = toBool_(getStringOption_(prefix + ":transition_quantification"));
    config.max_transition_pep = getDoubleOption_(prefix + ":max_transition_pep");
    config.ipf_max_peptidoform_pep = getDoubleOption_(prefix + ":ipf_max_peptidoform_pep");
    config.max_rs_peakgroup_qvalue = getDoubleOption_(prefix + ":max_rs_peakgroup_qvalue");
    config.peptide = toBool_(getStringOption_(prefix + ":use_peptide_scores"));
    config.max_global_peptide_qvalue = getDoubleOption_(prefix + ":max_global_peptide_qvalue");
    config.protein = toBool_(getStringOption_(prefix + ":use_protein_scores"));
    config.max_global_protein_qvalue = getDoubleOption_(prefix + ":max_global_protein_qvalue");
    config.gene = toBool_(getStringOption_(prefix + ":use_gene_scores"));
    config.max_global_gene_qvalue = getDoubleOption_(prefix + ":max_global_gene_qvalue");
    config.use_alignment = toBool_(getStringOption_(prefix + ":use_alignment"));
    config.max_alignment_pep = getDoubleOption_(prefix + ":max_alignment_pep");
    config.exclude_decoys = toBool_(getStringOption_(prefix + ":exclude_decoys"));
    return config;
  }

  OpenSwathResultsExportConfig getResultsConfig_() const
  {
    OpenSwathResultsExportConfig config;
    config.filters = getFilterConfig_("Export:results");
    config.format = toExportFormat_(getStringOption_("Export:results:format"));
    return config;
  }

  OpenSwathMatrixExportConfig getMatrixConfig_(const OpenSwathMatrixLevel level) const
  {
    OpenSwathMatrixExportConfig config;
    config.filters = getFilterConfig_("Export:matrix");
    config.level = level;
    config.normalization = toNormalization_(getStringOption_("Export:matrix:normalization"));
    const Int top_n = getIntOption_("Export:matrix:top_n");
    if (top_n < 0 || (level != OpenSwathMatrixLevel::Precursor && top_n < 1))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Parameter 'Export:matrix:top_n' must be >= 1 for peptide/protein/gene matrix export and >= 0 for precursor export.",
        StringUtils::toStr(top_n));
    }
    config.top_n = static_cast<Size>(top_n);
    config.consistent_top = toBool_(getStringOption_("Export:matrix:consistent_top"));
    config.format = toExportFormat_(getStringOption_("Export:matrix:format"));
    return config;
  }

  OpenSwathParquetExportConfig getParquetConfig_() const
  {
    OpenSwathParquetExportConfig config;
    config.include_transition_data = toBool_(getStringOption_("Export:parquet:include_transition_data"));
    config.filters.exclude_decoys = toBool_(getStringOption_("Export:parquet:exclude_decoys"));
    return config;
  }

  std::vector<ExportTask> getExportTasks_() const
  {
    std::vector<ExportTask> tasks;
    if (toBool_(getStringOption_("Export:parquet:run")))
    {
      tasks.push_back({ExportTaskType::FeatureParquet, std::nullopt});
    }
    if (toBool_(getStringOption_("Export:results:run")))
    {
      tasks.push_back({ExportTaskType::Results, std::nullopt});
    }
    if (toBool_(getStringOption_("Export:matrix:precursor")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Precursor});
    }
    if (toBool_(getStringOption_("Export:matrix:peptide")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Peptide});
    }
    if (toBool_(getStringOption_("Export:matrix:protein")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Protein});
    }
    if (toBool_(getStringOption_("Export:matrix:gene")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Gene});
    }
    return tasks;
  }

  static std::string exportTaskLabel_(const ExportTask& task)
  {
    switch (task.type)
    {
      case ExportTaskType::FeatureParquet: return "feature parquet export";
      case ExportTaskType::Results: return "results export";
      case ExportTaskType::Matrix: return toString(*task.matrix_level) + " matrix export";
    }
    return "export";
  }

  static std::string makeExportBasePath_(const StringList& input_files, const std::string& out_dir)
  {
    const std::string base_name = input_files.size() == 1 ? File::stemName(input_files.front()) : "OpenDIA";
    return out_dir + "/" + base_name;
  }

  static void logPeptidoformSummary_(const std::vector<IPFResultRow>& results)
  {
    const Size targets_at_1pct = static_cast<Size>(std::count_if(
      results.begin(), results.end(),
      [](const auto& row)
      {
        return std::isfinite(row.qvalue) && row.qvalue <= 0.01;
      }));
    OPENMS_LOG_INFO << targetLabelPlural_(InferenceLevel::Peptidoform) << " at 1% FDR: " << targets_at_1pct << std::endl;
  }

  static void logLevelContextSummary_(const std::vector<LevelContextInputRow>& input_rows,
                                      const std::vector<LevelContextResultRow>& results,
                                      const InferenceLevel level,
                                      const InferenceContext context,
                                      const std::map<Int64, std::string>& run_basenames)
  {
    constexpr double summary_fdr_threshold = 0.01;
    std::map<std::pair<Int64, Int64>, bool> is_target_by_key;
    for (const auto& row : input_rows)
    {
      is_target_by_key[{runKey_(row.run_id), row.entity_id}] = !row.decoy;
    }

    Size total_targets_at_1pct = 0;
    std::map<Int64, Size> per_run_targets_at_1pct;
    for (const auto& row : results)
    {
      if (!std::isfinite(row.qvalue) || row.qvalue > summary_fdr_threshold)
      {
        continue;
      }

      const auto target_it = is_target_by_key.find({runKey_(row.run_id), row.entity_id});
      if (target_it == is_target_by_key.end() || !target_it->second)
      {
        continue;
      }

      ++total_targets_at_1pct;
      if (context == InferenceContext::RunSpecific && row.run_id.has_value())
      {
        ++per_run_targets_at_1pct[*row.run_id];
      }
    }

    if (context != InferenceContext::RunSpecific)
    {
      OPENMS_LOG_INFO << targetLabelPlural_(level) << " at 1% FDR: " << total_targets_at_1pct << std::endl;
      return;
    }

    if (per_run_targets_at_1pct.size() != 1)
    {
      OPENMS_LOG_INFO << targetLabelPlural_(level) << " at 1% FDR across all runs: " << total_targets_at_1pct << std::endl;
    }

    for (const auto& run_count : per_run_targets_at_1pct)
    {
      const auto name_it = run_basenames.find(run_count.first);
      const std::string run_label = (name_it != run_basenames.end()) ? name_it->second : "RUN_ID " + StringUtils::toStr(run_count.first);
      OPENMS_LOG_INFO << run_label << ": " << run_count.second << " " << targetLabelPlural_(level) << std::endl;
    }
  }

  std::unordered_set<std::string> loadPriorityPeptideSequences_(const std::vector<std::string>& tsv_files,
                                                                const Param& tsv_reader_param)
  {
    std::unordered_set<std::string> priority_sequences;
    for (const auto& tsv_file : tsv_files)
    {
      if (tsv_file.empty() || !File::exists(tsv_file))
      {
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
      }
      catch (const Exception::BaseException& e)
      {
        OPENMS_LOG_WARN << "Failed to load priority peptide file " << tsv_file << ": " << e.what() << std::endl;
      }
    }
    return priority_sequences;
  }

  WorkingDirectory prepareWorkingDirectory_(const std::string& out_dir) const
  {
    WorkingDirectory working_dir;
    const bool keep_intermediate_files = toBool_(getStringOption_("workflow:keep_intermediate_files"));
    const std::string requested_dir = getOutputDirOption("workflow:intermediate_dir");

    if (!requested_dir.empty())
    {
      working_dir.path = File::absolutePath(requested_dir);
      File::makeDir(working_dir.path);
      working_dir.remove_on_success = !keep_intermediate_files;
      return working_dir;
    }

    if (keep_intermediate_files)
    {
      working_dir.path = File::absolutePath(out_dir + "/OpenDIA_intermediates");
      File::makeDir(working_dir.path);
      return working_dir;
    }

    working_dir.temp_dir = std::make_unique<File::TempDir>(true);
    working_dir.path = working_dir.temp_dir->getPath();
    working_dir.remove_on_success = true;
    return working_dir;
  }

  void cleanupWorkingDirectory_(const WorkingDirectory& working_dir) const
  {
    if (working_dir.remove_on_success && !working_dir.path.empty())
    {
      File::removeDirRecursively(working_dir.path);
    }
  }

  bool isPeptidoformInferenceRequested_() const
  {
    return toBool_(getStringOption_("Inference:peptidoform:run"));
  }

  ExitCodes runExtractionToOSW_(const std::string& prepared_library_pqp,
                                const std::string& extracted_osw,
                                const bool enable_uis_scoring)
  {
    Internal::ScopedResamplingWarningSuppression scoped_resampling_warning_suppression;

    StringList file_list = getStringList_("in");
    const std::string swath_windows_file = getStringOption_("TargetedDataExtraction:swath_windows_file");
    const std::string out_chrom = getStringOption_("TargetedDataExtraction:out_chrom");
    const std::string out_mobilogram = getStringOption_("TargetedDataExtraction:out_mobilogram");
    const bool split_file = getFlag_("TargetedDataExtraction:split_file_input");
    const bool use_emg_score = getFlag_("TargetedDataExtraction:use_elution_model_score");
    bool pasef = getFlag_("TargetedDataExtraction:pasef");
    const bool sort_swath_maps = getFlag_("TargetedDataExtraction:sort_swath_maps");
    const bool use_ms1_traces = getStringOption_("TargetedDataExtraction:enable_ms1") == "true";
    const int batchSize = static_cast<int>(getIntOption_("TargetedDataExtraction:batchSize"));
    const int outer_loop_threads = static_cast<int>(getIntOption_("TargetedDataExtraction:outer_loop_threads"));
    const int ms1_isotopes = static_cast<int>(getIntOption_("TargetedDataExtraction:ms1_isotopes"));
    const int innerBatchSize = static_cast<int>(getIntOption_("TargetedDataExtraction:innerBatchSize"));
    const int maxConcurrentSwaths = static_cast<int>(getIntOption_("TargetedDataExtraction:maxConcurrentSwaths"));
    const Size debug_level = static_cast<Size>(getIntOption_("debug"));
    Param debug_params = getParam_().copy("TargetedDataExtraction:", true);

    StringList disable_features = getStringList_("TargetedDataExtraction:disable_features");
    const bool disable_im_calibration = std::find(disable_features.begin(), disable_features.end(), "no_IM_calibration") != disable_features.end();
    const bool disable_im_windowing = std::find(disable_features.begin(), disable_features.end(), "no_IM_windowing") != disable_features.end();

    std::string readoptions = getStringOption_("TargetedDataExtraction:readOptions");
    const bool keep_cached_files = getFlag_("TargetedDataExtraction:keep_cached_files");
    const std::string tmp_dir = StringUtils::ensureLastChar(File::absolutePath(getStringOption_("TargetedDataExtraction:tempDirectory")), '/');

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

    Param irt_calibration_params = getParam_().copy("TargetedDataExtraction:Calibration:", true);
    bool auto_irt = irt_calibration_params.getValue("auto_irt:enabled").toString() == "true";
    std::string irt_tr_file = irt_calibration_params.getValue("files:linear_irt_file").toString();
    std::string priority_sampling_irt_tr_file = irt_calibration_params.getValue("tr_irt_priority_sampling").toString();
    std::string trafo_in = irt_calibration_params.getValue("rt_norm").toString();
    const UInt irt_bins_lin = irt_calibration_params.getValue("auto_irt:irt_bins");
    const UInt irt_pep_lin = irt_calibration_params.getValue("auto_irt:irt_peptides_per_bin");

    if (!irt_tr_file.empty())
    {
      if (auto_irt || !priority_sampling_irt_tr_file.empty())
      {
        OPENMS_LOG_WARN << "TargetedDataExtraction:Calibration:files:linear_irt_file provided -> disabling auto_irt and ignoring tr_irt_priority_sampling." << std::endl;
      }
      auto_irt = false;
      irt_calibration_params.setValue("auto_irt:enabled", "false");
      priority_sampling_irt_tr_file.clear();
      irt_calibration_params.setValue("tr_irt_priority_sampling", "");
    }

    if (trafo_in.empty() && irt_tr_file.empty() && !auto_irt)
    {
      OPENMS_LOG_INFO << "Neither rt_norm nor linear iRT peptides nor auto_irt are configured; a null RT transformation will be used." << std::endl;
    }

    if (auto_irt)
    {
      if (irt_bins_lin == 0)
      {
        writeLogError_("Parameter error: TargetedDataExtraction:Calibration:auto_irt:irt_bins must be > 0 when auto_irt is enabled.");
        return PARSE_ERROR;
      }
      if (irt_pep_lin == 0)
      {
        writeLogError_("Parameter error: TargetedDataExtraction:Calibration:auto_irt:irt_peptides_per_bin must be > 0 when auto_irt is enabled.");
        return PARSE_ERROR;
      }
    }

    if (!priority_sampling_irt_tr_file.empty())
    {
      if (!auto_irt)
      {
        OPENMS_LOG_WARN << "Priority iRT sampling file provided but auto_irt is disabled; ignoring: " << priority_sampling_irt_tr_file << std::endl;
      }
      else
      {
        if (!File::exists(priority_sampling_irt_tr_file))
        {
          writeLogError_("Parameter error: Priority iRT file does not exist: " + priority_sampling_irt_tr_file);
          return PARSE_ERROR;
        }
        if (FileHandler::getType(priority_sampling_irt_tr_file) != FileTypes::TSV)
        {
          writeLogError_("Parameter error: Priority iRT file must be in TSV format.");
          return PARSE_ERROR;
        }
      }
    }

    if (!swath_windows_file.empty())
    {
      std::vector<double> swath_prec_lower;
      std::vector<double> swath_prec_upper;
      SwathWindowLoader::readSwathWindows(swath_windows_file, swath_prec_lower, swath_prec_upper);
    }

    const double min_upper_edge_dist = getDoubleOption_("TargetedDataExtraction:min_upper_edge_dist");
    const bool use_ms1_im = getStringOption_("TargetedDataExtraction:use_ms1_ion_mobility") == "true";
    const bool prm = getStringOption_("TargetedDataExtraction:matching_window_only") == "true";

    ChromExtractParams cp;
    cp.min_upper_edge_dist = min_upper_edge_dist;
    cp.mz_extraction_window = getDoubleOption_("TargetedDataExtraction:mz_extraction_window");
    cp.ppm = getStringOption_("TargetedDataExtraction:mz_extraction_window_unit") == "ppm";
    cp.rt_extraction_window = getDoubleOption_("TargetedDataExtraction:rt_extraction_window");
    cp.im_extraction_window = getDoubleOption_("TargetedDataExtraction:ion_mobility_window");
    cp.extraction_function = getStringOption_("TargetedDataExtraction:extraction_function");
    cp.extra_rt_extract = getDoubleOption_("TargetedDataExtraction:extra_rt_extraction_window");

    ChromExtractParams cp_irt = cp;
    cp_irt.rt_extraction_window = -1;
    cp_irt.mz_extraction_window = getDoubleOption_("TargetedDataExtraction:irt_mz_extraction_window");
    cp_irt.im_extraction_window = getDoubleOption_("TargetedDataExtraction:irt_im_extraction_window");
    cp_irt.ppm = getStringOption_("TargetedDataExtraction:irt_mz_extraction_window_unit") == "ppm";

    if (disable_im_calibration)
    {
      cp_irt.im_extraction_window = -1;
    }
    else if ((cp_irt.im_extraction_window == -1) && (cp.im_extraction_window != -1))
    {
      OPENMS_LOG_WARN << "Warning: irt_im_extraction_window is not set, this will lead to no ion mobility calibration." << std::endl;
    }

    ChromExtractParams cp_ms1 = cp;
    cp_ms1.mz_extraction_window = getDoubleOption_("TargetedDataExtraction:mz_extraction_window_ms1");
    cp_ms1.ppm = getStringOption_("TargetedDataExtraction:mz_extraction_window_ms1_unit") == "ppm";
    cp_ms1.im_extraction_window = use_ms1_im ? getDoubleOption_("TargetedDataExtraction:im_extraction_window_ms1") : -1;

    if (disable_im_windowing)
    {
      cp.im_extraction_window = -1;
      cp_ms1.im_extraction_window = -1;
    }

    Param feature_finder_param = getParam_().copy("TargetedDataExtraction:Scoring:", true);
    feature_finder_param.setValue("use_ms1_ion_mobility", getStringOption_("TargetedDataExtraction:use_ms1_ion_mobility"));
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
    if (feature_finder_param.getValue("TransitionGroupPicker:compute_peak_shape_metrics").toBool())
    {
      feature_finder_param.setValue("Scores:use_peak_shape_metrics", "true");
    }

    Param tsv_reader_param = getParam_().copy("TargetedDataExtraction:Library:", true);

    Param tmp_mrm_map_param = getParam_().copy("TargetedDataExtraction:MRMMapping:", true);
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

    OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(FileTypes::PQP, prepared_library_pqp, tsv_reader_param);
    if (File::exists(extracted_osw) && !File::remove(extracted_osw))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, extracted_osw);
    }
    if (!File::copy(prepared_library_pqp, extracted_osw))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, extracted_osw);
    }

    std::unordered_set<std::string> priority_peptides;
    if (auto_irt)
    {
      std::vector<std::string> priority_files;
      const std::string data_path = File::getOpenMSDataPath();
      const std::string irtkit_path = data_path + "/CHEMISTRY/irtkit.tsv";
      const std::string cirtkit_path = data_path + "/CHEMISTRY/cirtkit.tsv";
      if (File::exists(irtkit_path)) priority_files.push_back(irtkit_path);
      if (File::exists(cirtkit_path)) priority_files.push_back(cirtkit_path);
      if (!priority_sampling_irt_tr_file.empty()) priority_files.push_back(priority_sampling_irt_tr_file);
      priority_peptides = loadPriorityPeptideSequences_(priority_files, TransitionTSVFile().getDefaults());
    }

    std::vector<StringList> run_groups;
    if (split_file)
    {
      run_groups.push_back(file_list);
    }
    else
    {
      for (const auto& file : file_list)
      {
        run_groups.push_back({file});
      }
    }

    const bool user_pasef = pasef;
    OpenSwathOSWWriter oswwriter(extracted_osw, enable_uis_scoring);
    oswwriter.writeHeader();

    Size run_index = 0;
    for (const auto& current_run_files : run_groups)
    {
      pasef = user_pasef;
      ChromExtractParams cp_current = cp;
      ChromExtractParams cp_ms1_current = cp_ms1;
      ChromExtractParams cp_irt_current = cp_irt;
      Param feature_finder_param_run = feature_finder_param;

      std::string per_run_tmp = tmp_dir;
      std::unique_ptr<File::TempDir> per_run_temp_dir;
      if (readoptions == "cache")
      {
        per_run_temp_dir = std::make_unique<File::TempDir>(tmp_dir, keep_cached_files);
        per_run_tmp = per_run_temp_dir->getPath();
      }

      std::shared_ptr<ExperimentalSettings> exp_meta(new ExperimentalSettings);
      std::vector<OpenSwath::SwathMap> swath_maps;
      std::vector<std::string> swath_map_sources;
      if (!loadSwathFiles(current_run_files, exp_meta, swath_maps, swath_map_sources, split_file, per_run_tmp,
                          readoptions, swath_windows_file, min_upper_edge_dist, getFlag_("force"),
                          sort_swath_maps, prm))
      {
        OPENMS_LOG_ERROR << "Failed to load SWATH files for run " << (run_index + 1) << "." << std::endl;
        return PARSE_ERROR;
      }

      if (!pasef)
      {
        const bool has_im_windows = std::any_of(swath_maps.begin(), swath_maps.end(),
          [](const OpenSwath::SwathMap& map)
          {
            return !map.ms1 && map.imLower >= 0 && map.imUpper >= 0;
          });
        if (has_im_windows)
        {
          pasef = true;
        }
      }
      if (disable_im_windowing && pasef)
      {
        pasef = false;
      }

      {
        const std::string im_score_setting = feature_finder_param_run.getValue("Scores:use_ion_mobility_scores").toString();
        if (im_score_setting == "auto")
        {
          feature_finder_param_run.setValue("Scores:use_ion_mobility_scores", pasef ? "true" : "false");
        }
      }

      if (pasef)
      {
        for (const auto& tr : transition_exp.getTransitions())
        {
          if (tr.precursor_im < 0)
          {
            throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
              "Transition " + tr.getNativeID() + " does not have a valid precursor ion mobility value required for PASEF extraction.");
          }
        }
      }

      std::string file_basename;
      if (run_groups.size() > 1)
      {
        file_basename = File::stemName(current_run_files[0]);
      }

      std::string irt_trafo_out = debug_params.getValue("irt_trafo").toString();
      if (!irt_trafo_out.empty() && run_groups.size() > 1)
      {
        irt_trafo_out = File::path(irt_trafo_out) + "/" + file_basename + "_" + File::basename(irt_trafo_out);
      }
      std::string irt_mzml_out = debug_params.getValue("irt_mzml").toString();
      if (!irt_mzml_out.empty() && run_groups.size() > 1)
      {
        irt_mzml_out = File::path(irt_mzml_out) + "/" + file_basename + "_" + File::basename(irt_mzml_out);
      }

      Param irt_detection_param = getParam_().copy("TargetedDataExtraction:Calibration:RTNormalization:", true);
      Param calibration_param = getParam_().copy("TargetedDataExtraction:Calibration:MassIMCorrection:", true);
      calibration_param.setValue("mz_extraction_window", cp_irt_current.mz_extraction_window);
      calibration_param.setValue("mz_extraction_window_ppm", cp_irt_current.ppm ? "true" : "false");
      calibration_param.setValue("im_extraction_window", cp_irt_current.im_extraction_window);

      bool mrm_mode = true;
      for (const auto& sm : swath_maps)
      {
        if (sm.ms1 || sm.sptr->getNrSpectra() > 0)
        {
          mrm_mode = false;
          break;
        }
      }

      TransformationDescription trafo_rtnorm;
      if (!trafo_in.empty())
      {
        FileHandler().loadTransformations(trafo_in, trafo_rtnorm, false, {FileTypes::TRANSFORMATIONXML});
        Param model_params;
        model_params.setValue("symmetric_regression", "false");
        model_params.setValue("span", irt_detection_param.getValue("lowess:span"));
        model_params.setValue("num_nodes", irt_detection_param.getValue("b_spline:num_nodes"));
        const std::string model_type = irt_detection_param.getValue("alignmentMethod").toString();
        trafo_rtnorm.fitModel(model_type, model_params);
      }
      else
      {
        CalibrationWorkflow calibration_wf;
        Param cal_params = irt_calibration_params.copy("", false);
        cal_params.remove("MassIMCorrection:");
        cal_params.remove("RTNormalization:");
        cal_params.remove("rt_norm");
        cal_params.remove("tr_irt_priority_sampling");
        calibration_wf.setParameters(cal_params);
        calibration_wf.setLogType(log_type_);
        IrtStrategy strategy = calibration_wf.determineIrtStrategy(transition_exp, run_groups.size());

        OpenSwath::LightTargetedExperiment prefiltered_irt_targets;
        const OpenSwath::LightTargetedExperiment* irt_sampling_transition_exp = &transition_exp;
        const bool auto_irt_prefilter_enabled = cal_params.getValue("auto_irt:prefilter:enabled").toBool();
        if (auto_irt && auto_irt_prefilter_enabled &&
            (strategy == IrtStrategy::SAMPLE_ONCE || strategy == IrtStrategy::SAMPLE_PER_RUN))
        {
          TransitionListEvidenceFilter evidence_filter;
          Param prefilter_params = cal_params.copy("auto_irt:prefilter:", true);
          evidence_filter.setParameters(prefilter_params);
          evidence_filter.setLogType(log_type_);
          try
          {
            auto prefilter_result = evidence_filter.filter(swath_maps, transition_exp, cp_ms1_current, cp_irt_current, pasef, outer_loop_threads);
            prefiltered_irt_targets = std::move(prefilter_result.filtered_targets);
            irt_sampling_transition_exp = &prefiltered_irt_targets;
            if (strategy == IrtStrategy::SAMPLE_ONCE)
            {
              strategy = IrtStrategy::SAMPLE_PER_RUN;
            }
          }
          catch (const Exception::IllegalArgument& e)
          {
            OPENMS_LOG_WARN << "Calibration:auto_irt:prefilter could not be applied (" << e.what()
                            << "); falling back to unfiltered auto-iRT sampling.\n";
          }
        }

        std::vector<std::string> priority_pep_strings(priority_peptides.begin(), priority_peptides.end());
        auto irt_experiments = calibration_wf.prepareIrtExperiments(strategy, *irt_sampling_transition_exp, priority_pep_strings, run_index);
        auto calibration_result = calibration_wf.performCalibration(
          swath_maps, transition_exp, cp_current, cp_ms1_current, irt_experiments,
          feature_finder_param_run, cp_irt_current, irt_detection_param, calibration_param,
          irt_mrm_map_param, pasef, load_into_memory, irt_trafo_out, irt_mzml_out, debug_level);
        trafo_rtnorm = calibration_result.rt_trafo;
      }

      const UInt64 cur_run = OpenMS::UniqueIdGenerator::getUniqueId();

      Interfaces::IMSDataConsumer* chromatogramConsumer = nullptr;
      std::string out_chrom_current = out_chrom;
      if (!out_chrom.empty() && run_groups.size() > 1)
      {
        out_chrom_current = File::path(out_chrom) + "/" + file_basename + "_" + File::basename(out_chrom);
      }
      prepareChromOutput(&chromatogramConsumer, exp_meta, transition_exp, out_chrom_current, cur_run, current_run_files[0]);

      std::unique_ptr<MobilogramParquetConsumer> mobilogramConsumer;
      std::string out_mobilogram_current = out_mobilogram;
      if (!out_mobilogram.empty() && run_groups.size() > 1)
      {
        out_mobilogram_current = File::path(out_mobilogram) + "/" + file_basename + "_" + File::basename(out_mobilogram);
      }
      prepareMobilogramOutput(mobilogramConsumer, exp_meta, transition_exp, out_mobilogram_current, cur_run, current_run_files[0]);

      oswwriter.addRun(cur_run, current_run_files[0]);
      if (auto* sql_cons = dynamic_cast<MSDataSqlConsumer*>(chromatogramConsumer))
      {
        sql_cons->addRun(current_run_files[0], cur_run);
        sql_cons->setRunId(cur_run);
      }
      oswwriter.setRunId(cur_run);

      FeatureMap run_feature_file;
      OpenSwathWorkflow wf(use_ms1_traces, use_ms1_im, prm, pasef, mrm_mode, outer_loop_threads);
      wf.setLogType(log_type_);
      wf.performExtraction(swath_maps, trafo_rtnorm, cp_current, cp_ms1_current, feature_finder_param_run,
                           transition_exp, run_feature_file, false, oswwriter, chromatogramConsumer,
                           batchSize, ms1_isotopes, load_into_memory, mrm_map_param,
                           mobilogramConsumer.get(), innerBatchSize, maxConcurrentSwaths);

      swath_maps.clear();
      if (mobilogramConsumer)
      {
        mobilogramConsumer->finalize();
      }
      delete chromatogramConsumer;
      ++run_index;
    }

    return EXECUTION_OK;
  }

  void runInference_(const std::string& workflow_osw)
  {
    const auto tasks = getInferenceTasks_();
    if (tasks.empty())
    {
      return;
    }

    OSWFile osw(workflow_osw);
    const PeptidoformInferenceConfig peptidoform_config = getPeptidoformConfig_();
    const ErrorEstimationConfig error_config = getErrorConfig_();
    const bool has_run_specific_task = std::any_of(tasks.begin(), tasks.end(),
      [](const auto& task)
      {
        return task.context.has_value() && *task.context == InferenceContext::RunSpecific;
      });
    const std::map<Int64, std::string> run_basenames = has_run_specific_task ? osw.readRunBasenames() : std::map<Int64, std::string>{};

    for (const auto& task : tasks)
    {
      ProgressLogger progress_logger;
      progress_logger.setLogType(ProgressLogger::CMD);
      progress_logger.startProgress(0, 1, inferenceTaskLabel_(task));

      if (task.level == InferenceLevel::Peptidoform)
      {
        const auto precursor_rows = osw.readIPFPrecursorData(peptidoform_config);
        const auto transition_rows = osw.readIPFTransitionData(peptidoform_config);
        const auto alignment_rows = osw.readIPFAlignmentData(peptidoform_config);

        if (transition_rows.empty())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Peptidoform inference was requested but no IPF transition evidence is available in the workflow OSW.");
        }

        OpenSwathPeptidoformInference inference;
        const auto results = inference.infer(precursor_rows, transition_rows, alignment_rows, peptidoform_config);
        osw.writeIPFResults("", results);
        logPeptidoformSummary_(results);
        progress_logger.endProgress();
        continue;
      }

      LevelContextInferenceConfig config;
      config.level = task.level;
      config.context = *task.context;
      config.error = error_config;
      const auto input_rows = osw.readLevelContextData(task.level, config.context);

      std::vector<LevelContextResultRow> results;
      if (task.level == InferenceLevel::Peptide)
      {
        OpenSwathPeptideInference inference;
        results = inference.infer(input_rows, config);
      }
      else if (task.level == InferenceLevel::Protein)
      {
        OpenSwathProteinInference inference;
        results = inference.infer(input_rows, config);
      }
      else
      {
        OpenSwathGeneInference inference;
        results = inference.infer(input_rows, config);
      }

      osw.writeLevelContextResults("", task.level, config.context, results);
      logLevelContextSummary_(input_rows, results, task.level, config.context, run_basenames);
      progress_logger.endProgress();
    }
  }

  void runExports_(const std::string& input_osw, const StringList& input_files, const std::string& out_dir)
  {
    const auto tasks = getExportTasks_();
    if (tasks.empty())
    {
      return;
    }

    File::makeDir(out_dir);
    const std::string base_path = makeExportBasePath_(input_files, out_dir);
    OSWFile osw(input_osw);
    std::optional<std::vector<OpenSwathExportRow>> matrix_rows_cache;

    for (const auto& task : tasks)
    {
      ProgressLogger progress_logger;
      progress_logger.setLogType(ProgressLogger::CMD);
      progress_logger.startProgress(0, 1, exportTaskLabel_(task));

      switch (task.type)
      {
        case ExportTaskType::FeatureParquet:
        {
          const auto parquet_config = getParquetConfig_();
          const auto feature_table = osw.readOpenSwathFeatureScoreTable(parquet_config);
          const std::string feature_out = base_path + ".precursor.feature.scores.parquet";
          OpenSwathParquetExporter::writeFeatureScores(feature_out, feature_table);
          if (parquet_config.include_transition_data)
          {
            const auto transition_table = osw.readOpenSwathTransitionScoreTable(parquet_config);
            if (!transition_table.rows.empty())
            {
              const std::string transition_out = base_path + ".transition.feature.scores.parquet";
              OpenSwathParquetExporter::writeTransitionScores(transition_out, transition_table);
            }
          }
          break;
        }
        case ExportTaskType::Results:
        {
          const auto results_config = getResultsConfig_();
          const auto rows = osw.readOpenSwathExportRows(results_config.filters);
          const std::string suffix = results_config.format == OpenSwathExportFileFormat::Parquet ? ".results.parquet" : ".results.tsv";
          OpenSwathResultsExporter::write(base_path + suffix, rows, results_config);
          break;
        }
        case ExportTaskType::Matrix:
        {
          const auto matrix_config = getMatrixConfig_(*task.matrix_level);
          if (!matrix_rows_cache.has_value())
          {
            matrix_rows_cache = osw.readOpenSwathExportRows(matrix_config.filters);
          }
          const auto matrix = OpenSwathMatrixExporter::buildMatrix(*matrix_rows_cache, matrix_config);
          const std::string suffix = matrix_config.format == OpenSwathExportFileFormat::Parquet ? ".matrix.parquet" : ".matrix.tsv";
          OpenSwathMatrixExporter::writeMatrix(base_path + "." + toString(*task.matrix_level) + suffix, matrix, matrix_config);
          break;
        }
      }

      progress_logger.endProgress();
    }
  }

  ExitCodes main_(int, const char**) override
  {
    const StringList input_files = getStringList_("in");
    const std::string out_dir = File::absolutePath(getOutputDirOption("out_dir"));
    const std::string input_library = getStringOption_("tr");

    FileTypes::Type tr_type = FileTypes::nameToType(getStringOption_("tr_type"));
    if (tr_type == FileTypes::UNKNOWN)
    {
      tr_type = FileHandler::getType(input_library);
    }
    if (tr_type == FileTypes::UNKNOWN)
    {
      writeLogError_("Error: Could not determine input file type for '-tr'.");
      return PARSE_ERROR;
    }

    WorkingDirectory working_dir = prepareWorkingDirectory_(out_dir);
    OPENMS_LOG_INFO << "OpenDIA working directory: " << working_dir.path << std::endl;

    try
    {
      File::makeDir(out_dir);
      const std::string prepared_library_pqp = working_dir.path + "/prepared_library.pqp";
      const std::string workflow_osw = working_dir.path + "/workflow.osw";

      OpenSwathLibraryPreparation library_preparation;
      library_preparation.setLogType(log_type_);
      const Param reader_parameters = getParam_().copy("TargetedDataExtraction:Library:", true);
      OpenSwathLibraryPreparation::LibraryStats library_stats;

      const LibraryMode requested_library_mode = getLibraryMode_();
      LibraryMode resolved_library_mode = requested_library_mode;
      bool prepared_library_ready = false;

      if (requested_library_mode == LibraryMode::AUTO)
      {
        OPENMS_LOG_INFO << "Auto-detecting OpenDIA library mode from input decoy coverage.\n";
        library_stats = library_preparation.normalizeLibraryToPQP(input_library, tr_type, prepared_library_pqp, reader_parameters);
        if (library_stats.hasDecoys())
        {
          resolved_library_mode = LibraryMode::PREPARED;
          prepared_library_ready = true;
          OPENMS_LOG_INFO << "Auto-detected prepared library input because decoy transitions are already present.\n";
        }
        else
        {
          resolved_library_mode = LibraryMode::EMPIRICAL;
          OPENMS_LOG_INFO << "Auto-detected empirical library input because no decoy transitions were found. Running peptide query preparation.\n";
          if (File::exists(prepared_library_pqp) && !File::remove(prepared_library_pqp))
          {
            throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prepared_library_pqp);
          }
        }
      }

      if (resolved_library_mode == LibraryMode::PREPARED)
      {
        if (!prepared_library_ready)
        {
          OPENMS_LOG_INFO << "Using prepared library mode.\n";
          library_stats = library_preparation.normalizeLibraryToPQP(input_library, tr_type, prepared_library_pqp, reader_parameters);
        }
      }
      else
      {
        OPENMS_LOG_INFO << "Using empirical library mode: running peptide query preparation.\n";
        const auto assay_parameters = getAssayGeneratorParameters_();
        const auto decoy_parameters = getDecoyGeneratorParameters_();
        if (isPeptidoformInferenceRequested_() && !assay_parameters.enable_ipf)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Peptidoform inference was requested, but PeptideQueryParameters:AssayGenerator:enable_ipf is false. Enable IPF-capable empirical library preparation first.");
        }
        library_stats = library_preparation.prepareEmpiricalLibraryToPQP(
          input_library, tr_type, prepared_library_pqp, assay_parameters, decoy_parameters, reader_parameters, working_dir.path);
      }

      if (!library_stats.hasDecoys())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "The prepared library does not contain decoy transitions. OpenDIA rescoring requires a target/decoy library.");
      }
      if (isPeptidoformInferenceRequested_() && !library_stats.hasIdentifyingTransitions())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Peptidoform inference was requested, but the prepared library does not contain identifying/IPF transitions.");
      }

      bool enable_uis_scoring = getStringOption_("TargetedDataExtraction:enable_ipf") == "true";
      if (isPeptidoformInferenceRequested_() && !enable_uis_scoring)
      {
        OPENMS_LOG_INFO << "Peptidoform inference requested: overriding TargetedDataExtraction:enable_ipf to true for extraction-time UIS scoring." << std::endl;
        enable_uis_scoring = true;
      }

      const ExitCodes extraction_result = runExtractionToOSW_(prepared_library_pqp, workflow_osw, enable_uis_scoring);
      if (extraction_result != EXECUTION_OK)
      {
        return extraction_result;
      }

      Rescorer rescoring;
      Param rescoring_params = getParam_().copy("Rescoring:", true);
      rescoring_params.remove("level");
      rescoring.setParameters(rescoring_params.copySubset(rescoring.getDefaults()));

      auto rescoring_levels = parseRescoreLevels_(getStringList_("Rescoring:level"));
      if (rescoring_levels.empty())
      {
        rescoring_levels.push_back(RescoreLevel::MS1MS2);
      }
      if (isPeptidoformInferenceRequested_() &&
          std::find(rescoring_levels.begin(), rescoring_levels.end(), RescoreLevel::TRANSITION) == rescoring_levels.end())
      {
        OPENMS_LOG_INFO << "Peptidoform inference requested: appending transition-level rescoring." << std::endl;
        rescoring_levels.push_back(RescoreLevel::TRANSITION);
      }

      for (const auto level : rescoring_levels)
      {
        rescoring.score(workflow_osw, level);
      }

      runInference_(workflow_osw);
      runExports_(workflow_osw, input_files, out_dir);

      cleanupWorkingDirectory_(working_dir);
      return EXECUTION_OK;
    }
    catch (const Exception::BaseException&)
    {
      OPENMS_LOG_ERROR << "OpenDIA failed. Intermediate workflow files were preserved in: " << working_dir.path << std::endl;
      throw;
    }
  }

private:
  Rescorer rescoring_defaults_;
};

int main(int argc, const char** argv)
{
  TOPPOpenDIA tool;
  return tool.main(argc, argv);
}
/// @endcond
