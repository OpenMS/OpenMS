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
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWParquetWriter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathOSWWriter.h>
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
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include <OpenMS/FORMAT/ArrowSchemaRegistry.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/FORMAT/SqliteConnector_impl.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "../openms/source/ANALYSIS/OPENSWATH/OpenSwathCanonicalLibraryMappingHelper.h"

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
1. extraction/scoring to a working `workflow.osw` or `workflow.oswpq`
2. in-process Percolator rescoring in the same working container
3. peptide/protein/gene/peptidoform inference in the same working workflow
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

  enum class WorkflowFormat
  {
    OSW,
    OSWPQ
  };

  struct InferenceTask
  {
    InferenceLevel level = InferenceLevel::Peptidoform;
    std::optional<InferenceContext> context;
  };

  using OptionalDoubleMember = std::optional<double> OpenSwathFeatureScoreRow::*;

  struct InferenceScoreColumn
  {
    const char* name = "";
    OptionalDoubleMember member = nullptr;
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

  struct OSWPQWorkspace
  {
    std::string output_path;
    std::string base_dir;
    bool archive_input = false;
    bool dirty = false;
    std::unique_ptr<File::TempDir> temp_dir;
  };

  struct PreparedLibraryPrecursor_
  {
    Int64 precursor_id = -1;
    std::string traml_id;
    std::string group_label;
    double precursor_mz = 0.0;
    Int32 charge = 0;
    std::optional<double> library_intensity;
    std::optional<double> library_rt;
    std::optional<double> library_drift_time;
    bool decoy = false;
  };

  struct PreparedLibraryPeptide_
  {
    Int64 peptide_id = -1;
    std::string unmodified_sequence;
    std::string modified_sequence;
    bool decoy = false;
  };

  struct PreparedLibraryProtein_
  {
    Int64 protein_id = -1;
    std::string accession;
    bool decoy = false;
  };

  struct PreparedLibraryGene_
  {
    Int64 gene_id = -1;
    std::string name;
    std::optional<bool> decoy;
  };

  struct PreparedLibraryTransition_
  {
    Int64 transition_id = -1;
    std::vector<Int64> precursor_ids;
    std::string traml_id;
    double product_mz = 0.0;
    Int32 charge = 0;
    std::string type;
    Int32 ordinal = 0;
    std::string annotation;
    bool detecting = false;
    std::optional<double> library_intensity;
    bool decoy = false;
    std::vector<Int64> peptide_ids;
  };

  struct PreparedLibraryLookup_
  {
    std::unordered_map<Int64, PreparedLibraryPrecursor_> precursors;
    std::unordered_map<Int64, PreparedLibraryPeptide_> peptides;
    std::unordered_map<Int64, PreparedLibraryProtein_> proteins;
    std::unordered_map<Int64, PreparedLibraryGene_> genes;
    std::unordered_map<Int64, std::vector<Int64>> precursor_to_peptides;
    std::unordered_map<Int64, std::vector<Int64>> peptide_to_proteins;
    std::unordered_map<Int64, std::vector<Int64>> peptide_to_genes;
    std::unordered_map<Int64, std::string> protein_names_by_peptide;
    std::unordered_map<Int64, std::string> gene_names_by_peptide;
    std::unordered_map<Int64, Int64> unique_protein_by_peptide;
    std::unordered_map<Int64, Int64> unique_gene_by_peptide;
    std::unordered_map<Int64, PreparedLibraryTransition_> transitions;
  };

  struct LevelContextResultMaps_
  {
    std::unordered_map<Int64, LevelContextResultRow> global;
    std::map<std::pair<Int64, Int64>, LevelContextResultRow> experiment_wide;
    std::map<std::pair<Int64, Int64>, LevelContextResultRow> run_specific;
  };

  struct FeatureTransitionObservation_
  {
    std::optional<Int64> run_id;
    std::optional<Int64> feature_id;
    std::vector<std::optional<double>> values;
    std::optional<double> score;
    std::optional<Int32> rank;
    std::optional<double> pvalue;
    std::optional<double> qvalue;
    std::optional<double> pep;
  };

  struct TransitionAggregation_
  {
    std::vector<std::string> areas;
    std::vector<std::string> apices;
    std::vector<std::string> annotations;
  };

  struct ExportQValueMaps_
  {
    std::unordered_map<Int64, double> global;
    std::map<std::pair<Int64, Int64>, double> experiment_wide;
    std::map<std::pair<Int64, Int64>, double> run_specific;
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

  static std::pair<bool, std::vector<InferenceScoreColumn>> getInferenceScoreColumns_(const std::vector<InferenceTask>& tasks)
  {
    bool include_ipf_peptide_id = false;
    std::vector<InferenceScoreColumn> columns;

    const auto append_column = [&](const char* name, OptionalDoubleMember member)
    {
      const auto duplicate = std::find_if(columns.begin(), columns.end(),
        [&](const auto& column)
        {
          return std::string_view(column.name) == name;
        });
      if (duplicate == columns.end())
      {
        columns.push_back({name, member});
      }
    };

    for (const auto& task : tasks)
    {
      switch (task.level)
      {
        case InferenceLevel::Peptidoform:
          include_ipf_peptide_id = true;
          append_column("score_ipf_precursor_peakgroup_pep", &OpenSwathFeatureScoreRow::score_ipf_precursor_peakgroup_pep);
          append_column("score_ipf_pep", &OpenSwathFeatureScoreRow::score_ipf_pep);
          append_column("score_ipf_qvalue", &OpenSwathFeatureScoreRow::score_ipf_qvalue);
          break;

        case InferenceLevel::Peptide:
          if (task.context == InferenceContext::Global)
          {
            append_column("score_peptide_global_score", &OpenSwathFeatureScoreRow::score_peptide_global_score);
            append_column("score_peptide_global_pvalue", &OpenSwathFeatureScoreRow::score_peptide_global_pvalue);
            append_column("score_peptide_global_qvalue", &OpenSwathFeatureScoreRow::score_peptide_global_qvalue);
            append_column("score_peptide_global_pep", &OpenSwathFeatureScoreRow::score_peptide_global_pep);
          }
          else if (task.context == InferenceContext::ExperimentWide)
          {
            append_column("score_peptide_experiment_wide_score", &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_score);
            append_column("score_peptide_experiment_wide_pvalue", &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pvalue);
            append_column("score_peptide_experiment_wide_qvalue", &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_qvalue);
            append_column("score_peptide_experiment_wide_pep", &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pep);
          }
          else if (task.context == InferenceContext::RunSpecific)
          {
            append_column("score_peptide_run_specific_score", &OpenSwathFeatureScoreRow::score_peptide_run_specific_score);
            append_column("score_peptide_run_specific_pvalue", &OpenSwathFeatureScoreRow::score_peptide_run_specific_pvalue);
            append_column("score_peptide_run_specific_qvalue", &OpenSwathFeatureScoreRow::score_peptide_run_specific_qvalue);
            append_column("score_peptide_run_specific_pep", &OpenSwathFeatureScoreRow::score_peptide_run_specific_pep);
          }
          break;

        case InferenceLevel::Protein:
          if (task.context == InferenceContext::Global)
          {
            append_column("score_protein_global_score", &OpenSwathFeatureScoreRow::score_protein_global_score);
            append_column("score_protein_global_pvalue", &OpenSwathFeatureScoreRow::score_protein_global_pvalue);
            append_column("score_protein_global_qvalue", &OpenSwathFeatureScoreRow::score_protein_global_qvalue);
            append_column("score_protein_global_pep", &OpenSwathFeatureScoreRow::score_protein_global_pep);
          }
          else if (task.context == InferenceContext::ExperimentWide)
          {
            append_column("score_protein_experiment_wide_score", &OpenSwathFeatureScoreRow::score_protein_experiment_wide_score);
            append_column("score_protein_experiment_wide_pvalue", &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pvalue);
            append_column("score_protein_experiment_wide_qvalue", &OpenSwathFeatureScoreRow::score_protein_experiment_wide_qvalue);
            append_column("score_protein_experiment_wide_pep", &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pep);
          }
          else if (task.context == InferenceContext::RunSpecific)
          {
            append_column("score_protein_run_specific_score", &OpenSwathFeatureScoreRow::score_protein_run_specific_score);
            append_column("score_protein_run_specific_pvalue", &OpenSwathFeatureScoreRow::score_protein_run_specific_pvalue);
            append_column("score_protein_run_specific_qvalue", &OpenSwathFeatureScoreRow::score_protein_run_specific_qvalue);
            append_column("score_protein_run_specific_pep", &OpenSwathFeatureScoreRow::score_protein_run_specific_pep);
          }
          break;

        case InferenceLevel::Gene:
          if (task.context == InferenceContext::Global)
          {
            append_column("score_gene_global_score", &OpenSwathFeatureScoreRow::score_gene_global_score);
            append_column("score_gene_global_pvalue", &OpenSwathFeatureScoreRow::score_gene_global_pvalue);
            append_column("score_gene_global_qvalue", &OpenSwathFeatureScoreRow::score_gene_global_qvalue);
            append_column("score_gene_global_pep", &OpenSwathFeatureScoreRow::score_gene_global_pep);
          }
          else if (task.context == InferenceContext::ExperimentWide)
          {
            append_column("score_gene_experiment_wide_score", &OpenSwathFeatureScoreRow::score_gene_experiment_wide_score);
            append_column("score_gene_experiment_wide_pvalue", &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pvalue);
            append_column("score_gene_experiment_wide_qvalue", &OpenSwathFeatureScoreRow::score_gene_experiment_wide_qvalue);
            append_column("score_gene_experiment_wide_pep", &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pep);
          }
          else if (task.context == InferenceContext::RunSpecific)
          {
            append_column("score_gene_run_specific_score", &OpenSwathFeatureScoreRow::score_gene_run_specific_score);
            append_column("score_gene_run_specific_pvalue", &OpenSwathFeatureScoreRow::score_gene_run_specific_pvalue);
            append_column("score_gene_run_specific_qvalue", &OpenSwathFeatureScoreRow::score_gene_run_specific_qvalue);
            append_column("score_gene_run_specific_pep", &OpenSwathFeatureScoreRow::score_gene_run_specific_pep);
          }
          break;
      }
    }

    return {include_ipf_peptide_id, columns};
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
    registerStringOption_("workflow:working_format", "<choice>", "osw", "Internal workflow container used after extraction/scoring.", false);
    setValidStrings_("workflow:working_format", {"osw", "oswpq"});
    registerStringOption_("workflow:keep_intermediate_files", "<true|false>", "false", "Whether to retain prepared_library.pqp and the single working workflow.osw or workflow.oswpq after success.", false);
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

  WorkflowFormat getWorkflowFormat_() const
  {
    return getStringOption_("workflow:working_format") == "oswpq" ? WorkflowFormat::OSWPQ : WorkflowFormat::OSW;
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

  static bool parquetValuePresent_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    return array != nullptr && !array->IsNull(row);
  }

  static OpenSwathOSWWriter::OSWValue parquetInt64Value_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    return parquetValuePresent_(array, row) ?
      OpenSwathOSWWriter::OSWValue(ParquetFile::getInt64(array, row, 0, false)) :
      OpenSwathOSWWriter::OSWValue::null();
  }

  static OpenSwathOSWWriter::OSWValue parquetDoubleValue_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    return parquetValuePresent_(array, row) ?
      OpenSwathOSWWriter::OSWValue(ParquetFile::getDouble(array, row, 0.0, false)) :
      OpenSwathOSWWriter::OSWValue::null();
  }

  static std::optional<double> parquetOptionalDouble_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    if (!parquetValuePresent_(array, row))
    {
      return std::nullopt;
    }
    return ParquetFile::getDouble(array, row, 0.0, false);
  }

  static std::optional<Int64> parquetOptionalInt64_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    if (!parquetValuePresent_(array, row))
    {
      return std::nullopt;
    }
    return ParquetFile::getInt64(array, row, 0, false);
  }

  static std::shared_ptr<arrow::Array> getOptionalParquetColumn_(
    const std::unordered_map<std::string, std::shared_ptr<arrow::Array>>& columns,
    const std::string& name)
  {
    const auto it = columns.find(name);
    return it != columns.end() ? it->second : nullptr;
  }

  static void replaceParquetColumns_(const std::string& file_path,
                                     const std::unordered_set<std::string>& columns_to_replace,
                                     const std::vector<std::shared_ptr<arrow::Field>>& extra_fields,
                                     const std::vector<std::shared_ptr<arrow::Array>>& extra_arrays)
  {
    auto table = ParquetFile::readTable(file_path);
    std::vector<std::shared_ptr<arrow::Field>> fields;
    std::vector<std::shared_ptr<arrow::Array>> arrays;
    fields.reserve(table->num_columns() + static_cast<int>(extra_fields.size()));
    arrays.reserve(table->num_columns() + static_cast<int>(extra_arrays.size()));

    for (int i = 0; i < table->num_columns(); ++i)
    {
      const auto field = table->field(i);
      if (columns_to_replace.contains(field->name()))
      {
        continue;
      }
      fields.push_back(field);
      arrays.push_back(table->column(i)->chunk(0));
    }

    for (Size i = 0; i < extra_fields.size(); ++i)
    {
      fields.push_back(extra_fields[i]);
      arrays.push_back(extra_arrays[i]);
    }

    ParquetFile::writeTable(arrow::Table::Make(arrow::schema(fields), arrays), file_path);
  }

  static void removeExistingPath_(const std::string& path)
  {
    auto remove_file = [&](const std::string& file_path)
    {
      if (!File::exists(file_path))
      {
        return;
      }
      if (File::isDirectory(file_path))
      {
        File::removeDirRecursively(file_path);
        return;
      }
      if (!File::remove(file_path))
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, file_path);
      }
    };

    if (File::exists(path) && File::isDirectory(path))
    {
      File::removeDirRecursively(path);
      return;
    }

    remove_file(path);

    for (const char* suffix : {"-journal", "-wal", "-shm"})
    {
      remove_file(path + suffix);
    }
  }

  static void checkpointSQLiteDatabase_(const std::string& sqlite_path, const std::string& label)
  {
    const std::string wal_path = sqlite_path + "-wal";
    const UInt64 wal_size_before = File::exists(wal_path) ? File::fileSize(wal_path) : 0;
    OPENMS_LOG_INFO << "Checkpointing " << label;
    if (wal_size_before > 0)
    {
      OPENMS_LOG_INFO << " (WAL " << (wal_size_before / (1024 * 1024)) << " MiB)";
    }
    OPENMS_LOG_INFO << "." << std::endl;

    try
    {
      SqliteConnector conn(sqlite_path, SqliteConnector::SqlOpenMode::READWRITE);
      conn.executeStatement("PRAGMA wal_checkpoint(TRUNCATE);");
      conn.executeStatement("PRAGMA journal_mode=DELETE;");
    }
    catch (const Exception::BaseException& e)
    {
      OPENMS_LOG_WARN << "Failed to checkpoint " << label << ": " << e.getMessage() << std::endl;
      return;
    }

    const UInt64 wal_size_after = File::exists(wal_path) ? File::fileSize(wal_path) : 0;
    OPENMS_LOG_INFO << "Finished checkpointing " << label;
    if (wal_size_after > 0)
    {
      OPENMS_LOG_INFO << " (remaining WAL " << (wal_size_after / (1024 * 1024)) << " MiB)";
    }
    OPENMS_LOG_INFO << "." << std::endl;
  }

  static void sqliteCheckResult_(const int rc, sqlite3* db, const int expected, const char* action)
  {
    if (rc == expected)
    {
      return;
    }
    throw Exception::SqlOperationFailed(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        std::string(action) + " failed: " + sqlite3_errmsg(db));
  }

  static void sqlitePrepareStatement_(sqlite3* db, sqlite3_stmt** stmt, const std::string& statement)
  {
    sqliteCheckResult_(sqlite3_prepare_v2(db, statement.c_str(), -1, stmt, nullptr),
                       db, SQLITE_OK, "sqlite3_prepare_v2");
  }

  static std::string inferenceContextValue_(const InferenceContext context)
  {
    switch (context)
    {
      case InferenceContext::Global:
        return "global";
      case InferenceContext::ExperimentWide:
        return "experiment-wide";
      case InferenceContext::RunSpecific:
        return "run-specific";
    }
    return "global";
  }


  static std::string entityIdColumnName_(const InferenceLevel level)
  {
    switch (level)
    {
      case InferenceLevel::Peptide: return "peptide_id";
      case InferenceLevel::Protein: return "protein_id";
      case InferenceLevel::Gene: return "gene_id";
      case InferenceLevel::Peptidoform: break;
    }
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Direct OSWPQ level-context helpers do not support peptidoform inference.");
  }

  static std::string inferenceParquetPath_(const OSWPQWorkspace& workspace, const InferenceLevel level)
  {
    return workspace.base_dir + "/inference/score_" + toString(level) + ".parquet";
  }

  static std::vector<std::string> sqliteTableColumns_(SqliteConnector& conn, const std::string& table_name)
  {
    sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
    sqlite3_stmt* stmt = nullptr;
    sqlitePrepareStatement_(db, &stmt, "PRAGMA table_info('" + table_name + "');");
    std::vector<std::string> columns;
    int step_result = SQLITE_DONE;
    while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW)
    {
      const auto* text = sqlite3_column_text(stmt, 1);
      columns.emplace_back(text != nullptr ? reinterpret_cast<const char*>(text) : "");
    }
    sqliteCheckResult_(step_result, db, SQLITE_DONE, "sqlite3_step");
    sqlite3_finalize(stmt);
    return columns;
  }

  static bool sqliteTableHasColumn_(SqliteConnector& conn, const std::string& table_name, const std::string& column_name)
  {
    const auto columns = sqliteTableColumns_(conn, table_name);
    return std::find(columns.begin(), columns.end(), column_name) != columns.end();
  }

  static std::string parquetGeneNameSelect_(SqliteConnector& conn)
  {
    if (!conn.tableExists("GENE"))
    {
      return "NULL";
    }
    if (sqliteTableHasColumn_(conn, "GENE", "GENE_NAME"))
    {
      return "GENE.GENE_NAME";
    }
    return "CAST(GENE.ID AS TEXT)";
  }

  static std::string parquetGeneDecoySelect_(SqliteConnector& conn)
  {
    if (!conn.tableExists("GENE"))
    {
      return "NULL";
    }
    if (sqliteTableHasColumn_(conn, "GENE", "DECOY"))
    {
      return "GENE.DECOY";
    }
    return "NULL";
  }

  static std::string parquetTransitionAnnotationSelect_(SqliteConnector& conn)
  {
    const std::string fallback_annotation =
      "TRANSITION.TYPE || CAST(TRANSITION.ORDINAL AS TEXT) || '^' || CAST(TRANSITION.CHARGE AS TEXT)";
    if (sqliteTableHasColumn_(conn, "TRANSITION", "ANNOTATION"))
    {
      return "COALESCE(TRANSITION.ANNOTATION, " + fallback_annotation + ")";
    }
    return fallback_annotation;
  }

  static bool betterLevelContextResult_(const LevelContextResultRow& candidate, const LevelContextResultRow& current)
  {
    if (candidate.pep != current.pep)
    {
      return candidate.pep < current.pep;
    }
    if (candidate.qvalue != current.qvalue)
    {
      return candidate.qvalue < current.qvalue;
    }
    if (candidate.score != current.score)
    {
      return candidate.score > current.score;
    }
    return candidate.entity_id < current.entity_id;
  }

  static void appendUniqueEntity_(std::vector<Int64>& entities, const Int64 entity_id)
  {
    if (std::find(entities.begin(), entities.end(), entity_id) == entities.end())
    {
      entities.push_back(entity_id);
    }
  }

  static std::string joinStrings_(const std::vector<std::string>& values)
  {
    std::string joined;
    for (Size i = 0; i < values.size(); ++i)
    {
      if (i != 0)
      {
        joined += ";";
      }
      joined += values[i];
    }
    return joined;
  }

  void finalizePreparedLibraryLookup_(PreparedLibraryLookup_& lookup) const
  {
    for (auto& [precursor_id, peptide_ids] : lookup.precursor_to_peptides)
    {
      std::sort(peptide_ids.begin(), peptide_ids.end());
      peptide_ids.erase(std::unique(peptide_ids.begin(), peptide_ids.end()), peptide_ids.end());
    }

    for (auto& [peptide_id, protein_ids] : lookup.peptide_to_proteins)
    {
      std::sort(protein_ids.begin(), protein_ids.end());
      protein_ids.erase(std::unique(protein_ids.begin(), protein_ids.end()), protein_ids.end());

      std::vector<std::string> names;
      names.reserve(protein_ids.size());
      for (const Int64 protein_id : protein_ids)
      {
        const auto protein_it = lookup.proteins.find(protein_id);
        if (protein_it != lookup.proteins.end())
        {
          names.push_back(protein_it->second.accession);
        }
      }
      lookup.protein_names_by_peptide[peptide_id] = joinStrings_(names);
      if (protein_ids.size() == 1)
      {
        lookup.unique_protein_by_peptide[peptide_id] = protein_ids.front();
      }
    }

    for (auto& [peptide_id, gene_ids] : lookup.peptide_to_genes)
    {
      std::sort(gene_ids.begin(), gene_ids.end());
      gene_ids.erase(std::unique(gene_ids.begin(), gene_ids.end()), gene_ids.end());

      std::vector<std::string> names;
      names.reserve(gene_ids.size());
      for (const Int64 gene_id : gene_ids)
      {
        const auto gene_it = lookup.genes.find(gene_id);
        if (gene_it != lookup.genes.end())
        {
          names.push_back(gene_it->second.name);
        }
      }
      lookup.gene_names_by_peptide[peptide_id] = joinStrings_(names);
      if (gene_ids.size() == 1)
      {
        lookup.unique_gene_by_peptide[peptide_id] = gene_ids.front();
      }
    }

    for (auto& [transition_id, transition] : lookup.transitions)
    {
      std::sort(transition.peptide_ids.begin(), transition.peptide_ids.end());
      transition.peptide_ids.erase(std::unique(transition.peptide_ids.begin(), transition.peptide_ids.end()), transition.peptide_ids.end());
    }
  }

  PreparedLibraryLookup_ buildPreparedLibraryLookupFromLightTargetedExperiment_(
    const OpenSwath::LightTargetedExperiment& targeted_exp,
    const bool load_transition_metadata) const
  {
    PreparedLibraryLookup_ lookup;
    const auto canonical_mapping = Internal::buildOpenSwathCanonicalLibraryMapping(targeted_exp);

    std::vector<std::string> peptide_sequences;
    peptide_sequences.reserve(targeted_exp.compounds.size() + targeted_exp.transitions.size());
    for (const auto& compound : targeted_exp.compounds)
    {
      if (compound.isPeptide())
      {
        peptide_sequences.push_back(compound.sequence);
      }
    }
    for (const auto& transition : targeted_exp.transitions)
    {
      for (const auto& peptidoform : transition.peptidoforms)
      {
        peptide_sequences.push_back(peptidoform);
      }
    }
    std::sort(peptide_sequences.begin(), peptide_sequences.end());
    peptide_sequences.erase(std::unique(peptide_sequences.begin(), peptide_sequences.end()), peptide_sequences.end());

    std::unordered_map<std::string, Int64> peptide_ids_by_sequence;
    peptide_ids_by_sequence.reserve(peptide_sequences.size());
    for (Size i = 0; i < peptide_sequences.size(); ++i)
    {
      const auto& modified_sequence = peptide_sequences[i];
      std::string unmodified_sequence;
      try
      {
        unmodified_sequence = AASequence::fromString(modified_sequence).toUnmodifiedString();
      }
      catch (Exception::InvalidValue&)
      {
        unmodified_sequence = modified_sequence;
      }

      const Int64 peptide_id = static_cast<Int64>(i);
      peptide_ids_by_sequence.emplace(modified_sequence, peptide_id);
      lookup.peptides.emplace(peptide_id, PreparedLibraryPeptide_{peptide_id, unmodified_sequence, modified_sequence, false});
    }

    std::vector<std::string> protein_accessions;
    protein_accessions.reserve(targeted_exp.proteins.size());
    for (const auto& protein : targeted_exp.proteins)
    {
      protein_accessions.push_back(protein.id);
    }
    std::sort(protein_accessions.begin(), protein_accessions.end());
    protein_accessions.erase(std::unique(protein_accessions.begin(), protein_accessions.end()), protein_accessions.end());

    std::unordered_map<std::string, Int64> protein_ids_by_accession;
    protein_ids_by_accession.reserve(protein_accessions.size());
    for (Size i = 0; i < protein_accessions.size(); ++i)
    {
      const auto& accession = protein_accessions[i];
      const Int64 protein_id = static_cast<Int64>(i);
      protein_ids_by_accession.emplace(accession, protein_id);
      lookup.proteins.emplace(protein_id, PreparedLibraryProtein_{protein_id, accession, false});
    }

    std::unordered_map<std::string, Int64> gene_ids_by_name;
    gene_ids_by_name.reserve(targeted_exp.compounds.size());

    for (const auto& compound : targeted_exp.compounds)
    {
      const auto precursor_id_it = canonical_mapping.compound_to_precursor.find(compound.id);
      if (precursor_id_it == canonical_mapping.compound_to_precursor.end())
      {
        continue;
      }

      const Int64 precursor_id = precursor_id_it->second;
      PreparedLibraryPrecursor_ precursor;
      precursor.precursor_id = precursor_id;
      precursor.traml_id = compound.id;
      precursor.group_label = compound.isPeptide() ? compound.peptide_group_label : "";
      const auto precursor_mz_it = canonical_mapping.precursor_mz_by_id.find(precursor_id);
      precursor.precursor_mz = precursor_mz_it != canonical_mapping.precursor_mz_by_id.end() ? precursor_mz_it->second : 0.0;
      precursor.charge = static_cast<Int32>(compound.charge);
      if (std::isfinite(compound.rt))
      {
        precursor.library_rt = compound.rt;
      }
      if (compound.drift_time != -1 && std::isfinite(compound.drift_time))
      {
        precursor.library_drift_time = compound.drift_time;
      }
      const auto precursor_decoy_it = canonical_mapping.precursor_decoy_by_id.find(precursor_id);
      precursor.decoy = precursor_decoy_it != canonical_mapping.precursor_decoy_by_id.end() ?
        precursor_decoy_it->second :
        (compound.hasDecoy() ? compound.getDecoy() : false);
      lookup.precursors[precursor_id] = std::move(precursor);

      if (!compound.isPeptide())
      {
        continue;
      }

      const auto peptide_id_it = peptide_ids_by_sequence.find(compound.sequence);
      if (peptide_id_it == peptide_ids_by_sequence.end())
      {
        continue;
      }

      const Int64 peptide_id = peptide_id_it->second;
      appendUniqueEntity_(lookup.precursor_to_peptides[precursor_id], peptide_id);
      lookup.peptides[peptide_id].decoy = lookup.peptides[peptide_id].decoy || lookup.precursors.at(precursor_id).decoy;

      auto& protein_ids = lookup.peptide_to_proteins[peptide_id];
      for (const auto& protein_ref : compound.protein_refs)
      {
        const auto protein_id_it = protein_ids_by_accession.find(protein_ref);
        if (protein_id_it != protein_ids_by_accession.end())
        {
          appendUniqueEntity_(protein_ids, protein_id_it->second);
        }
      }

      const std::string gene_name = compound.gene_name.empty() ? "NA" : compound.gene_name;
      auto [gene_it, inserted] = gene_ids_by_name.try_emplace(gene_name, static_cast<Int64>(gene_ids_by_name.size()));
      if (inserted)
      {
        lookup.genes.emplace(gene_it->second, PreparedLibraryGene_{gene_it->second, gene_name, false});
      }
      appendUniqueEntity_(lookup.peptide_to_genes[peptide_id], gene_it->second);
    }

    for (const auto& [peptide_id, protein_ids] : lookup.peptide_to_proteins)
    {
      if (lookup.peptides[peptide_id].decoy)
      {
        for (const Int64 protein_id : protein_ids)
        {
          lookup.proteins[protein_id].decoy = true;
        }
      }
    }

    for (const auto& [peptide_id, gene_ids] : lookup.peptide_to_genes)
    {
      if (lookup.peptides[peptide_id].decoy)
      {
        for (const Int64 gene_id : gene_ids)
        {
          lookup.genes[gene_id].decoy = true;
        }
      }
    }

    if (load_transition_metadata)
    {
      lookup.transitions.reserve(targeted_exp.transitions.size());
      for (Size i = 0; i < targeted_exp.transitions.size(); ++i)
      {
        const auto& transition = targeted_exp.transitions[i];
        const auto precursor_id_it = canonical_mapping.compound_to_precursor.find(transition.peptide_ref);
        if (precursor_id_it == canonical_mapping.compound_to_precursor.end())
        {
          continue;
        }

        PreparedLibraryTransition_ transition_entry;
        transition_entry.transition_id = static_cast<Int64>(i);
        transition_entry.precursor_ids.push_back(precursor_id_it->second);
        transition_entry.traml_id = transition.transition_name;
        transition_entry.product_mz = transition.product_mz;
        transition_entry.charge = static_cast<Int32>(transition.fragment_charge);
        const std::string fragment_type = transition.getFragmentType();
        transition_entry.type = fragment_type.empty() ? "" : StringUtils::substr(fragment_type, 0, 1);
        transition_entry.ordinal = static_cast<Int32>(transition.fragment_nr);
        transition_entry.annotation = transition.getAnnotation();
        transition_entry.detecting = transition.isDetectingTransition();
        transition_entry.library_intensity = transition.library_intensity;
        transition_entry.decoy = transition.getDecoy();

        for (const auto& peptidoform : transition.peptidoforms)
        {
          const auto peptide_id_it = peptide_ids_by_sequence.find(peptidoform);
          if (peptide_id_it != peptide_ids_by_sequence.end())
          {
            appendUniqueEntity_(transition_entry.peptide_ids, peptide_id_it->second);
          }
        }

        lookup.transitions.emplace(transition_entry.transition_id, std::move(transition_entry));
      }
    }

    finalizePreparedLibraryLookup_(lookup);
    return lookup;
  }

  PreparedLibraryLookup_ loadPreparedLibraryLookup_(const std::string& prepared_library_pqp,
                                                    const bool load_transition_metadata) const
  {
    PreparedLibraryLookup_ lookup;
    SqliteConnector conn(prepared_library_pqp, SqliteConnector::SqlOpenMode::READ_ONLY);
    sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
    const auto text_at = [](sqlite3_stmt* stmt, const int column) -> std::string
    {
      const auto* text = sqlite3_column_text(stmt, column);
      return text != nullptr ? reinterpret_cast<const char*>(text) : "";
    };

    {
      const bool has_library_drift_time = sqliteTableHasColumn_(conn, "PRECURSOR", "LIBRARY_DRIFT_TIME");
      sqlite3_stmt* stmt = nullptr;
      sqlitePrepareStatement_(db, &stmt,
        "SELECT ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_INTENSITY, LIBRARY_RT, "
        + std::string(has_library_drift_time ? "LIBRARY_DRIFT_TIME" : "NULL") +
        ", DECOY FROM PRECURSOR;");
      int step_result = SQLITE_DONE;
      while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW)
      {
        PreparedLibraryPrecursor_ precursor;
        precursor.precursor_id = sqlite3_column_int64(stmt, 0);
        precursor.traml_id = text_at(stmt, 1);
        precursor.group_label = text_at(stmt, 2);
        precursor.precursor_mz = sqlite3_column_double(stmt, 3);
        precursor.charge = static_cast<Int32>(sqlite3_column_int(stmt, 4));
        if (sqlite3_column_type(stmt, 5) != SQLITE_NULL) precursor.library_intensity = sqlite3_column_double(stmt, 5);
        if (sqlite3_column_type(stmt, 6) != SQLITE_NULL) precursor.library_rt = sqlite3_column_double(stmt, 6);
        if (sqlite3_column_type(stmt, 7) != SQLITE_NULL) precursor.library_drift_time = sqlite3_column_double(stmt, 7);
        precursor.decoy = sqlite3_column_int(stmt, 8) != 0;
        lookup.precursors[precursor.precursor_id] = std::move(precursor);
      }
      sqliteCheckResult_(step_result, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(stmt);
    }

    {
      sqlite3_stmt* stmt = nullptr;
      sqlitePrepareStatement_(db, &stmt,
        "SELECT ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY FROM PEPTIDE;");
      int step_result = SQLITE_DONE;
      while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW)
      {
        PreparedLibraryPeptide_ peptide;
        peptide.peptide_id = sqlite3_column_int64(stmt, 0);
        peptide.unmodified_sequence = text_at(stmt, 1);
        peptide.modified_sequence = text_at(stmt, 2);
        peptide.decoy = sqlite3_column_int(stmt, 3) != 0;
        lookup.peptides[peptide.peptide_id] = std::move(peptide);
      }
      sqliteCheckResult_(step_result, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(stmt);
    }

    {
      sqlite3_stmt* stmt = nullptr;
      sqlitePrepareStatement_(db, &stmt,
        "SELECT PRECURSOR_ID, PEPTIDE_ID FROM PRECURSOR_PEPTIDE_MAPPING ORDER BY PRECURSOR_ID, PEPTIDE_ID;");
      int step_result = SQLITE_DONE;
      while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW)
      {
        lookup.precursor_to_peptides[sqlite3_column_int64(stmt, 0)].push_back(sqlite3_column_int64(stmt, 1));
      }
      sqliteCheckResult_(step_result, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(stmt);
    }

    if (conn.tableExists("PROTEIN") && conn.tableExists("PEPTIDE_PROTEIN_MAPPING"))
    {
      sqlite3_stmt* protein_stmt = nullptr;
      sqlitePrepareStatement_(db, &protein_stmt,
        "SELECT ID, PROTEIN_ACCESSION, DECOY FROM PROTEIN;");
      int protein_step = SQLITE_DONE;
      while ((protein_step = sqlite3_step(protein_stmt)) == SQLITE_ROW)
      {
        PreparedLibraryProtein_ protein;
        protein.protein_id = sqlite3_column_int64(protein_stmt, 0);
        protein.accession = text_at(protein_stmt, 1);
        protein.decoy = sqlite3_column_int(protein_stmt, 2) != 0;
        lookup.proteins[protein.protein_id] = std::move(protein);
      }
      sqliteCheckResult_(protein_step, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(protein_stmt);

      sqlite3_stmt* mapping_stmt = nullptr;
      sqlitePrepareStatement_(db, &mapping_stmt,
        "SELECT PEPTIDE_ID, PROTEIN_ID FROM PEPTIDE_PROTEIN_MAPPING ORDER BY PEPTIDE_ID, PROTEIN_ID;");
      int mapping_step = SQLITE_DONE;
      while ((mapping_step = sqlite3_step(mapping_stmt)) == SQLITE_ROW)
      {
        lookup.peptide_to_proteins[sqlite3_column_int64(mapping_stmt, 0)].push_back(sqlite3_column_int64(mapping_stmt, 1));
      }
      sqliteCheckResult_(mapping_step, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(mapping_stmt);
    }

    if (conn.tableExists("GENE") && conn.tableExists("PEPTIDE_GENE_MAPPING"))
    {
      sqlite3_stmt* gene_stmt = nullptr;
      sqlitePrepareStatement_(db, &gene_stmt,
        "SELECT ID, " + parquetGeneNameSelect_(conn) + ", " + parquetGeneDecoySelect_(conn) + " FROM GENE;");
      int gene_step = SQLITE_DONE;
      while ((gene_step = sqlite3_step(gene_stmt)) == SQLITE_ROW)
      {
        PreparedLibraryGene_ gene;
        gene.gene_id = sqlite3_column_int64(gene_stmt, 0);
        gene.name = text_at(gene_stmt, 1);
        if (sqlite3_column_type(gene_stmt, 2) != SQLITE_NULL) gene.decoy = sqlite3_column_int(gene_stmt, 2) != 0;
        lookup.genes[gene.gene_id] = std::move(gene);
      }
      sqliteCheckResult_(gene_step, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(gene_stmt);

      sqlite3_stmt* mapping_stmt = nullptr;
      sqlitePrepareStatement_(db, &mapping_stmt,
        "SELECT PEPTIDE_ID, GENE_ID FROM PEPTIDE_GENE_MAPPING ORDER BY PEPTIDE_ID, GENE_ID;");
      int mapping_step = SQLITE_DONE;
      while ((mapping_step = sqlite3_step(mapping_stmt)) == SQLITE_ROW)
      {
        lookup.peptide_to_genes[sqlite3_column_int64(mapping_stmt, 0)].push_back(sqlite3_column_int64(mapping_stmt, 1));
      }
      sqliteCheckResult_(mapping_step, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(mapping_stmt);
    }

    if (load_transition_metadata && conn.tableExists("TRANSITION") && conn.tableExists("TRANSITION_PRECURSOR_MAPPING"))
    {
      sqlite3_stmt* transition_stmt = nullptr;
      sqlitePrepareStatement_(db, &transition_stmt,
        "SELECT TRANSITION.ID, TRANSITION_PRECURSOR_MAPPING.PRECURSOR_ID, TRANSITION.TRAML_ID, "
        "TRANSITION.PRODUCT_MZ, TRANSITION.CHARGE, TRANSITION.TYPE, TRANSITION.ORDINAL, "
        + parquetTransitionAnnotationSelect_(conn) +
        ", TRANSITION.DETECTING, TRANSITION.LIBRARY_INTENSITY, TRANSITION.DECOY "
        "FROM TRANSITION "
        "INNER JOIN TRANSITION_PRECURSOR_MAPPING ON TRANSITION.ID = TRANSITION_PRECURSOR_MAPPING.TRANSITION_ID;");
      int transition_step = SQLITE_DONE;
      while ((transition_step = sqlite3_step(transition_stmt)) == SQLITE_ROW)
      {
        const Int64 transition_id = sqlite3_column_int64(transition_stmt, 0);
        auto& transition = lookup.transitions[transition_id];
        transition.transition_id = transition_id;
        const Int64 precursor_id = sqlite3_column_int64(transition_stmt, 1);
        if (std::find(transition.precursor_ids.begin(), transition.precursor_ids.end(), precursor_id) == transition.precursor_ids.end())
        {
          transition.precursor_ids.push_back(precursor_id);
        }
        transition.traml_id = text_at(transition_stmt, 2);
        transition.product_mz = sqlite3_column_double(transition_stmt, 3);
        transition.charge = static_cast<Int32>(sqlite3_column_int(transition_stmt, 4));
        transition.type = text_at(transition_stmt, 5);
        transition.ordinal = static_cast<Int32>(sqlite3_column_int(transition_stmt, 6));
        transition.annotation = text_at(transition_stmt, 7);
        transition.detecting = sqlite3_column_int(transition_stmt, 8) != 0;
        if (sqlite3_column_type(transition_stmt, 9) != SQLITE_NULL) transition.library_intensity = sqlite3_column_double(transition_stmt, 9);
        transition.decoy = sqlite3_column_int(transition_stmt, 10) != 0;
      }
      sqliteCheckResult_(transition_step, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(transition_stmt);

      if (conn.tableExists("TRANSITION_PEPTIDE_MAPPING"))
      {
        sqlite3_stmt* peptide_mapping_stmt = nullptr;
        sqlitePrepareStatement_(db, &peptide_mapping_stmt,
          "SELECT TRANSITION_ID, PEPTIDE_ID FROM TRANSITION_PEPTIDE_MAPPING ORDER BY TRANSITION_ID, PEPTIDE_ID;");
        int peptide_mapping_step = SQLITE_DONE;
        while ((peptide_mapping_step = sqlite3_step(peptide_mapping_stmt)) == SQLITE_ROW)
        {
          const Int64 transition_id = sqlite3_column_int64(peptide_mapping_stmt, 0);
          const Int64 peptide_id = sqlite3_column_int64(peptide_mapping_stmt, 1);
          lookup.transitions[transition_id].peptide_ids.push_back(peptide_id);
        }
        sqliteCheckResult_(peptide_mapping_step, db, SQLITE_DONE, "sqlite3_step");
        sqlite3_finalize(peptide_mapping_stmt);
      }
    }

    finalizePreparedLibraryLookup_(lookup);
    return lookup;
  }

  void augmentLookupFromOSWPQLibrary_(const OSWPQWorkspace& workspace,
                                      PreparedLibraryLookup_& lookup,
                                      const bool load_transition_metadata) const
  {
    const std::string precursors_path = workspace.base_dir + "/library/precursors.parquet";
    if (!File::exists(precursors_path))
    {
      return;
    }

    const auto source_precursors = lookup.precursors;
    const auto source_precursor_to_peptides = lookup.precursor_to_peptides;
    const auto source_transitions = lookup.transitions;
    lookup.precursors.clear();
    lookup.precursor_to_peptides.clear();
    lookup.transitions.clear();

    std::unordered_map<std::string, Int64> source_precursor_ids_by_key;
    source_precursor_ids_by_key.reserve(source_precursors.size() * 2);
    for (const auto& [precursor_id, precursor] : source_precursors)
    {
      source_precursor_ids_by_key.try_emplace(StringUtils::toStr(precursor_id), precursor_id);
      if (!precursor.traml_id.empty())
      {
        source_precursor_ids_by_key.try_emplace(precursor.traml_id, precursor_id);
      }
    }

    std::unordered_map<std::string, Int64> source_transition_ids_by_key;
    if (load_transition_metadata)
    {
      source_transition_ids_by_key.reserve(source_transitions.size() * 2);
      for (const auto& [transition_id, transition] : source_transitions)
      {
        source_transition_ids_by_key.try_emplace(StringUtils::toStr(transition_id), transition_id);
        if (!transition.traml_id.empty())
        {
          source_transition_ids_by_key.try_emplace(transition.traml_id, transition_id);
        }
      }
    }

    Int64 next_peptide_id = 1;
    for (const auto& [peptide_id, peptide] : lookup.peptides)
    {
      next_peptide_id = std::max(next_peptide_id, peptide_id + 1);
    }
    Int64 next_protein_id = 1;
    std::unordered_map<std::string, Int64> protein_ids_by_accession;
    protein_ids_by_accession.reserve(lookup.proteins.size());
    for (const auto& [protein_id, protein] : lookup.proteins)
    {
      next_protein_id = std::max(next_protein_id, protein_id + 1);
      protein_ids_by_accession[protein.accession] = protein_id;
    }

    auto precursors_table = ParquetFile::readTable(precursors_path);
    const auto precursor_id_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::PRECURSOR_ID);
    const auto precursor_mz_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::PRECURSOR_MZ);
    const auto charge_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::CHARGE);
    const auto library_rt_col = ParquetFile::getColumn(precursors_table, OSWPrecursorSchema::LIBRARY_RT);
    const auto drift_time_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::LIBRARY_DRIFT_TIME);
    const auto traml_id_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::TRAML_ID);
    const auto decoy_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::DECOY);
    const auto modified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::MODIFIED_SEQUENCE);
    const auto unmodified_sequence_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::UNMODIFIED_SEQUENCE);
    const auto protein_accessions_col = ParquetFile::getOptionalColumn(precursors_table, OSWPrecursorSchema::PROTEIN_ACCESSIONS);

    for (int64_t row = 0; row < precursors_table->num_rows(); ++row)
    {
      const Int64 precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
      const std::string traml_id = ParquetFile::getString(traml_id_col, row);
      std::optional<Int64> source_precursor_id;
      if (!traml_id.empty())
      {
        const auto source_it = source_precursor_ids_by_key.find(traml_id);
        if (source_it != source_precursor_ids_by_key.end())
        {
          source_precursor_id = source_it->second;
        }
      }
      if (!source_precursor_id.has_value())
      {
        const auto source_it = source_precursor_ids_by_key.find(StringUtils::toStr(precursor_id));
        if (source_it != source_precursor_ids_by_key.end())
        {
          source_precursor_id = source_it->second;
        }
      }

      PreparedLibraryPrecursor_ precursor;
      if (source_precursor_id.has_value())
      {
        const auto source_it = source_precursors.find(*source_precursor_id);
        if (source_it != source_precursors.end())
        {
          precursor = source_it->second;
        }
      }
      precursor.precursor_id = precursor_id;
      precursor.traml_id = traml_id;
      precursor.group_label = traml_id;
      precursor.precursor_mz = ParquetFile::getDouble(precursor_mz_col, row, 0.0, false);
      precursor.charge = static_cast<Int32>(ParquetFile::getInt64(charge_col, row, 0, false));
      precursor.library_rt = ParquetFile::getDouble(library_rt_col, row, 0.0, true);
      precursor.library_drift_time.reset();
      if (drift_time_col != nullptr && !drift_time_col->IsNull(row))
      {
        precursor.library_drift_time = ParquetFile::getDouble(drift_time_col, row, 0.0, false);
      }
      precursor.decoy = ParquetFile::getBool(decoy_col, row, false, true);
      lookup.precursors[precursor_id] = std::move(precursor);

      if (source_precursor_id.has_value())
      {
        const auto peptide_mapping_it = source_precursor_to_peptides.find(*source_precursor_id);
        if (peptide_mapping_it != source_precursor_to_peptides.end())
        {
          lookup.precursor_to_peptides[precursor_id] = peptide_mapping_it->second;
          continue;
        }
      }

      auto& peptide_ids = lookup.precursor_to_peptides[precursor_id];
      if (peptide_ids.empty())
      {
        PreparedLibraryPeptide_ peptide;
        peptide.peptide_id = next_peptide_id++;
        peptide.modified_sequence = ParquetFile::getString(modified_sequence_col, row);
        peptide.unmodified_sequence = ParquetFile::getString(unmodified_sequence_col, row);
        if (peptide.unmodified_sequence.empty())
        {
          peptide.unmodified_sequence = peptide.modified_sequence;
        }
        peptide.decoy = lookup.precursors.at(precursor_id).decoy;
        peptide_ids.push_back(peptide.peptide_id);
        lookup.peptides[peptide.peptide_id] = std::move(peptide);

        const std::vector<std::string> protein_accessions = ParquetFile::getStringList(protein_accessions_col, row);
        auto& mapped_proteins = lookup.peptide_to_proteins[peptide_ids.front()];
        mapped_proteins.reserve(protein_accessions.size());
        for (const auto& accession : protein_accessions)
        {
          if (accession.empty())
          {
            continue;
          }
          Int64 protein_id = -1;
          const auto protein_id_it = protein_ids_by_accession.find(accession);
          if (protein_id_it != protein_ids_by_accession.end())
          {
            protein_id = protein_id_it->second;
          }
          else
          {
            protein_id = next_protein_id++;
            protein_ids_by_accession[accession] = protein_id;
            lookup.proteins[protein_id] = {protein_id, accession, lookup.precursors.at(precursor_id).decoy};
          }
          mapped_proteins.push_back(protein_id);
        }
        lookup.protein_names_by_peptide[peptide_ids.front()] = joinStrings_(protein_accessions);
        if (mapped_proteins.size() == 1)
        {
          lookup.unique_protein_by_peptide[peptide_ids.front()] = mapped_proteins.front();
        }
      }
    }

    if (!load_transition_metadata)
    {
      return;
    }

    const std::string transitions_path = workspace.base_dir + "/library/transitions.parquet";
    if (!File::exists(transitions_path))
    {
      return;
    }

    auto transitions_table = ParquetFile::readTable(transitions_path);
    const auto transition_id_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::TRANSITION_ID);
    const auto transition_precursor_id_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::PRECURSOR_ID);
    const auto transition_traml_id_col = ParquetFile::getOptionalColumn(transitions_table, OSWTransitionSchema::TRAML_ID);
    const auto product_mz_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::PRODUCT_MZ);
    const auto transition_charge_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::CHARGE);
    const auto transition_type_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::TYPE);
    const auto transition_annotation_col = ParquetFile::getOptionalColumn(transitions_table, OSWTransitionSchema::ANNOTATION);
    const auto transition_ordinal_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::ORDINAL);
    const auto transition_detecting_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::DETECTING);
    const auto transition_intensity_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::LIBRARY_INTENSITY);
    const auto transition_decoy_col = ParquetFile::getColumn(transitions_table, OSWTransitionSchema::DECOY);

    for (int64_t row = 0; row < transitions_table->num_rows(); ++row)
    {
      const Int64 transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
      const Int64 precursor_id = ParquetFile::getInt64(transition_precursor_id_col, row, 0, false);
      const std::string transition_traml_id = ParquetFile::getString(transition_traml_id_col, row);
      PreparedLibraryTransition_ transition;
      if (!transition_traml_id.empty())
      {
        const auto source_it = source_transition_ids_by_key.find(transition_traml_id);
        if (source_it != source_transition_ids_by_key.end())
        {
          const auto transition_source_it = source_transitions.find(source_it->second);
          if (transition_source_it != source_transitions.end())
          {
            transition = transition_source_it->second;
          }
        }
      }
      transition.transition_id = transition_id;
      transition.precursor_ids.clear();
      transition.precursor_ids.push_back(precursor_id);
      transition.traml_id = transition_traml_id;
      transition.product_mz = ParquetFile::getDouble(product_mz_col, row, 0.0, false);
      transition.charge = static_cast<Int32>(ParquetFile::getInt64(transition_charge_col, row, 0, false));
      transition.type = ParquetFile::getString(transition_type_col, row);
      transition.ordinal = static_cast<Int32>(ParquetFile::getInt64(transition_ordinal_col, row, 0, false));
      transition.annotation = ParquetFile::getString(transition_annotation_col, row);
      transition.detecting = ParquetFile::getBool(transition_detecting_col, row, true, true);
      transition.library_intensity = ParquetFile::getDouble(transition_intensity_col, row, 0.0, true);
      transition.decoy = ParquetFile::getBool(transition_decoy_col, row, false, true);
      lookup.transitions[transition_id] = std::move(transition);

      const auto peptide_mapping_it = lookup.precursor_to_peptides.find(precursor_id);
      if (peptide_mapping_it != lookup.precursor_to_peptides.end())
      {
        auto& transition_entry = lookup.transitions[transition_id];
        for (const Int64 peptide_id : peptide_mapping_it->second)
        {
          if (std::find(transition_entry.peptide_ids.begin(), transition_entry.peptide_ids.end(), peptide_id) == transition_entry.peptide_ids.end())
          {
            transition_entry.peptide_ids.push_back(peptide_id);
          }
        }
      }
    }
  }

  static std::vector<Int64> mappedEntitiesForFeature_(const PreparedLibraryLookup_& lookup,
                                                      const InferenceLevel level,
                                                      const Int64 precursor_id)
  {
    std::vector<Int64> entities;
    const auto peptide_it = lookup.precursor_to_peptides.find(precursor_id);
    if (peptide_it == lookup.precursor_to_peptides.end())
    {
      return entities;
    }

    for (const Int64 peptide_id : peptide_it->second)
    {
      if (level == InferenceLevel::Peptide)
      {
        appendUniqueEntity_(entities, peptide_id);
        continue;
      }
      if (level == InferenceLevel::Protein)
      {
        const auto protein_it = lookup.unique_protein_by_peptide.find(peptide_id);
        if (protein_it != lookup.unique_protein_by_peptide.end())
        {
          appendUniqueEntity_(entities, protein_it->second);
        }
        continue;
      }
      if (level == InferenceLevel::Gene)
      {
        const auto gene_it = lookup.unique_gene_by_peptide.find(peptide_id);
        if (gene_it != lookup.unique_gene_by_peptide.end())
        {
          appendUniqueEntity_(entities, gene_it->second);
        }
      }
    }
    return entities;
  }

  static std::optional<LevelContextResultRow> selectBestLevelContextResult_(const std::vector<Int64>& entity_ids,
                                                                            const Int64 run_id,
                                                                            const LevelContextResultMaps_& maps,
                                                                            const InferenceContext context)
  {
    std::optional<LevelContextResultRow> best;
    for (const Int64 entity_id : entity_ids)
    {
      std::optional<LevelContextResultRow> candidate;
      if (context == InferenceContext::Global)
      {
        const auto it = maps.global.find(entity_id);
        if (it != maps.global.end()) candidate = it->second;
      }
      else if (context == InferenceContext::ExperimentWide)
      {
        const auto it = maps.experiment_wide.find({run_id, entity_id});
        if (it != maps.experiment_wide.end()) candidate = it->second;
      }
      else
      {
        const auto it = maps.run_specific.find({run_id, entity_id});
        if (it != maps.run_specific.end()) candidate = it->second;
      }

      if (!candidate.has_value())
      {
        continue;
      }
      if (!best.has_value() || betterLevelContextResult_(*candidate, *best))
      {
        best = candidate;
      }
    }
    return best;
  }

  static LevelContextResultMaps_ buildLevelContextResultMaps_(const std::vector<LevelContextResultRow>& results)
  {
    LevelContextResultMaps_ maps;
    for (const auto& row : results)
    {
      switch (row.context)
      {
        case InferenceContext::Global:
          maps.global[row.entity_id] = row;
          break;
        case InferenceContext::ExperimentWide:
          if (row.run_id.has_value())
          {
            maps.experiment_wide[{*row.run_id, row.entity_id}] = row;
          }
          break;
        case InferenceContext::RunSpecific:
          if (row.run_id.has_value())
          {
            maps.run_specific[{*row.run_id, row.entity_id}] = row;
          }
          break;
      }
    }
    return maps;
  }

  static void writeLevelContextResultsParquet_(const OSWPQWorkspace& workspace,
                                               const InferenceLevel level,
                                               const std::vector<LevelContextResultRow>& results)
  {
    const std::string inference_dir = workspace.base_dir + "/inference";
    File::makeDir(inference_dir);

    arrow::StringBuilder context_builder;
    arrow::Int64Builder run_id_builder;
    arrow::Int64Builder entity_id_builder;
    arrow::DoubleBuilder score_builder;
    arrow::DoubleBuilder pvalue_builder;
    arrow::DoubleBuilder qvalue_builder;
    arrow::DoubleBuilder pep_builder;

    for (const auto& row : results)
    {
      ParquetFile::appendOrThrow(context_builder.Append(inferenceContextValue_(row.context)), "context");
      if (row.run_id.has_value())
      {
        ParquetFile::appendOrThrow(run_id_builder.Append(*row.run_id), "run_id");
      }
      else
      {
        ParquetFile::appendOrThrow(run_id_builder.AppendNull(), "run_id");
      }
      ParquetFile::appendOrThrow(entity_id_builder.Append(row.entity_id), entityIdColumnName_(level));
      ParquetFile::appendOrThrow(score_builder.Append(row.score), "score");
      ParquetFile::appendOrThrow(pvalue_builder.Append(row.pvalue), "pvalue");
      ParquetFile::appendOrThrow(qvalue_builder.Append(row.qvalue), "qvalue");
      ParquetFile::appendOrThrow(pep_builder.Append(row.pep), "pep");
    }

    std::vector<std::shared_ptr<arrow::Field>> fields =
    {
      arrow::field("context", arrow::utf8(), false),
      arrow::field("run_id", arrow::int64(), true),
      arrow::field(entityIdColumnName_(level), arrow::int64(), false),
      arrow::field("score", arrow::float64(), false),
      arrow::field("pvalue", arrow::float64(), false),
      arrow::field("qvalue", arrow::float64(), false),
      arrow::field("pep", arrow::float64(), false)
    };
    std::vector<std::shared_ptr<arrow::Array>> arrays =
    {
      ParquetFile::finishArray(context_builder, "context"),
      ParquetFile::finishArray(run_id_builder, "run_id"),
      ParquetFile::finishArray(entity_id_builder, entityIdColumnName_(level)),
      ParquetFile::finishArray(score_builder, "score"),
      ParquetFile::finishArray(pvalue_builder, "pvalue"),
      ParquetFile::finishArray(qvalue_builder, "qvalue"),
      ParquetFile::finishArray(pep_builder, "pep")
    };

    ParquetFile::writeTable(arrow::Table::Make(arrow::schema(fields), arrays), inferenceParquetPath_(workspace, level));
  }

  static std::vector<LevelContextResultRow> readLevelContextResultsParquet_(const OSWPQWorkspace& workspace,
                                                                            const InferenceLevel level)
  {
    const std::string file_path = inferenceParquetPath_(workspace, level);
    if (!File::exists(file_path))
    {
      return {};
    }

    auto table = ParquetFile::readTable(file_path);
    const auto context_col = ParquetFile::getColumn(table, "context");
    const auto run_id_col = ParquetFile::getOptionalColumn(table, "run_id");
    const auto entity_id_col = ParquetFile::getColumn(table, entityIdColumnName_(level));
    const auto score_col = ParquetFile::getColumn(table, "score");
    const auto pvalue_col = ParquetFile::getColumn(table, "pvalue");
    const auto qvalue_col = ParquetFile::getColumn(table, "qvalue");
    const auto pep_col = ParquetFile::getColumn(table, "pep");

    std::vector<LevelContextResultRow> rows;
    rows.reserve(static_cast<Size>(table->num_rows()));
    for (int64_t row = 0; row < table->num_rows(); ++row)
    {
      LevelContextResultRow result;
      const std::string context = ParquetFile::getString(context_col, row);
      if (context == "global")
      {
        result.context = InferenceContext::Global;
      }
      else if (context == "experiment-wide")
      {
        result.context = InferenceContext::ExperimentWide;
      }
      else
      {
        result.context = InferenceContext::RunSpecific;
      }
      if (run_id_col != nullptr && !run_id_col->IsNull(row))
      {
        result.run_id = ParquetFile::getInt64(run_id_col, row, 0, false);
      }
      result.entity_id = ParquetFile::getInt64(entity_id_col, row, 0, false);
      result.score = ParquetFile::getDouble(score_col, row, 0.0, false);
      result.pvalue = ParquetFile::getDouble(pvalue_col, row, 1.0, false);
      result.qvalue = ParquetFile::getDouble(qvalue_col, row, 1.0, false);
      result.pep = ParquetFile::getDouble(pep_col, row, 1.0, false);
      rows.push_back(std::move(result));
    }
    return rows;
  }

  static std::map<Int64, std::string> readOSWPQRunBasenames_(const OSWPQWorkspace& workspace)
  {
    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    std::map<Int64, std::string> basenames;
    for (int64_t row = 0; row < runs_table->num_rows(); ++row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, row, 0, false);
      std::string filename = (filename_col != nullptr && !filename_col->IsNull(row)) ? ParquetFile::getString(filename_col, row) : "";
      std::string basename = File::stemName(filename);
      if (basename.empty())
      {
        basename = File::basename(filename);
      }
      if (basename.empty())
      {
        basename = "RUN_ID " + StringUtils::toStr(run_id);
      }
      basenames[run_id] = basename;
    }
    return basenames;
  }

  template <typename RowHandler>
  static Size readRunInferenceRows_(sqlite3* db,
                                    const std::string& statement,
                                    const Int64 run_id,
                                    std::unordered_map<Int64, OpenSwathFeatureScoreRow>& row_by_feature,
                                    RowHandler&& row_handler)
  {
    sqlite3_stmt* stmt = nullptr;
    sqlitePrepareStatement_(db, &stmt, statement);
    Size row_count = 0;
    int step_result = SQLITE_DONE;
    try
    {
      while ((step_result = sqlite3_step(stmt)) == SQLITE_ROW)
      {
        const Int64 feature_id = sqlite3_column_int64(stmt, 0);
        auto& score_row = row_by_feature[feature_id];
        score_row.run_id = run_id;
        score_row.feature_id = feature_id;
        row_handler(stmt, score_row);
        ++row_count;
      }
      sqliteCheckResult_(step_result, db, SQLITE_DONE, "sqlite3_step");
      sqlite3_finalize(stmt);
      return row_count;
    }
    catch (...)
    {
      if (stmt != nullptr)
      {
        sqlite3_finalize(stmt);
      }
      throw;
    }
  }

  static const std::array<const char*, 17>& featureMS1ParquetFields_()
  {
    static const std::array<const char*, 17> fields =
    {{
      "ms1_area_intensity",
      "ms1_apex_intensity",
      "ms1_exp_im",
      "ms1_delta_im",
      "var_ms1_massdev_score",
      "var_ms1_im_ms1_delta_score",
      "var_ms1_mi_score",
      "var_ms1_mi_contrast_score",
      "var_ms1_mi_combined_score",
      "var_ms1_isotope_correlation_score",
      "var_ms1_isotope_overlap_score",
      "var_ms1_xcorr_coelution",
      "var_ms1_xcorr_coelution_contrast",
      "var_ms1_xcorr_coelution_combined",
      "var_ms1_xcorr_shape",
      "var_ms1_xcorr_shape_contrast",
      "var_ms1_xcorr_shape_combined"
    }};
    return fields;
  }

  static const std::array<const char*, 37>& featureMS2ParquetFields_()
  {
    static const std::array<const char*, 37> fields =
    {{
      "ms2_area_intensity",
      "ms2_total_area_intensity",
      "ms2_apex_intensity",
      "ms2_exp_im",
      "ms2_exp_im_leftwidth",
      "ms2_exp_im_rightwidth",
      "ms2_delta_im",
      "ms2_total_mi",
      "var_ms2_bseries_score",
      "var_ms2_dotprod_score",
      "var_ms2_intensity_score",
      "var_ms2_isotope_correlation_score",
      "var_ms2_isotope_overlap_score",
      "var_ms2_library_corr",
      "var_ms2_library_dotprod",
      "var_ms2_library_manhattan",
      "var_ms2_library_rmsd",
      "var_ms2_library_rootmeansquare",
      "var_ms2_library_sangle",
      "var_ms2_log_sn_score",
      "var_ms2_manhattan_score",
      "var_ms2_massdev_score",
      "var_ms2_massdev_score_weighted",
      "var_ms2_mi_score",
      "var_ms2_mi_weighted_score",
      "var_ms2_mi_ratio_score",
      "var_ms2_norm_rt_score",
      "var_ms2_xcorr_coelution",
      "var_ms2_xcorr_coelution_weighted",
      "var_ms2_xcorr_shape",
      "var_ms2_xcorr_shape_weighted",
      "var_ms2_yseries_score",
      "var_ms2_elution_model_fit_score",
      "var_ms2_im_xcorr_shape",
      "var_ms2_im_xcorr_coelution",
      "var_ms2_im_delta_score",
      "var_ms2_im_log_intensity"
    }};
    return fields;
  }

  static const std::array<const char*, 41>& featureTransitionParquetFields_()
  {
    static const std::array<const char*, 41> fields =
    {{
      "area_intensity",
      "total_area_intensity",
      "apex_rt",
      "apex_intensity",
      "rt_fwhm",
      "masserror_ppm",
      "total_mi",
      "var_intensity_score",
      "var_intensity_ratio_score",
      "var_log_intensity",
      "var_xcorr_coelution",
      "var_xcorr_shape",
      "var_log_sn_score",
      "var_massdev_score",
      "var_mi_score",
      "var_mi_ratio_score",
      "var_isotope_correlation_score",
      "var_isotope_overlap_score",
      "exp_im",
      "exp_im_leftwidth",
      "exp_im_rightwidth",
      "delta_im",
      "var_im_delta_score",
      "var_im_log_intensity",
      "var_im_xcorr_coelution_contrast",
      "var_im_xcorr_shape_contrast",
      "var_im_xcorr_coelution_combined",
      "var_im_xcorr_shape_combined",
      "start_position_at_5",
      "end_position_at_5",
      "start_position_at_10",
      "end_position_at_10",
      "start_position_at_50",
      "end_position_at_50",
      "total_width",
      "tailing_factor",
      "asymmetry_factor",
      "slope_of_baseline",
      "baseline_delta_2_height",
      "points_across_baseline",
      "points_across_half_height"
    }};
    return fields;
  }

  OSWPQWorkspace prepareOSWPQWorkspace_(const std::string& workflow_oswpq) const
  {
    OSWPQWorkspace workspace;
    workspace.output_path = workflow_oswpq;
    workspace.archive_input = !File::isDirectory(workflow_oswpq);
    if (workspace.archive_input)
    {
      workspace.base_dir = ZipArchiveFile::unzipDirectory(workflow_oswpq, workspace.temp_dir);
    }
    else
    {
      workspace.base_dir = workflow_oswpq;
    }
    return workspace;
  }

  void commitOSWPQWorkspace_(OSWPQWorkspace& workspace) const
  {
    if (!workspace.dirty)
    {
      return;
    }
    if (workspace.archive_input)
    {
      OPENMS_LOG_INFO << "Repacking workflow.oswpq archive." << std::endl;
      removeExistingPath_(workspace.output_path);
      ZipArchiveFile::zipDirectory(workspace.base_dir, workspace.output_path);
      ZipArchiveFile::writeSidecarIndex(workspace.output_path);
      OPENMS_LOG_INFO << "Finished repacking workflow.oswpq archive." << std::endl;
    }
    workspace.dirty = false;
  }

  struct FeatureScoreInsertRow_
  {
    Int64 feature_id = -1;
    double score = 0.0;
    Int32 rank = 0;
    std::optional<double> pvalue;
    double qvalue = 1.0;
    double pep = 1.0;
  };

  struct TransitionScoreInsertRow_
  {
    Int64 feature_id = -1;
    Int64 transition_id = -1;
    double score = 0.0;
    Int32 rank = 0;
    std::optional<double> pvalue;
    double qvalue = 1.0;
    double pep = 1.0;
  };

  static void writeFeatureScoreTable_(const std::string& osw_path,
                                      const std::string& table_name,
                                      const std::vector<FeatureScoreInsertRow_>& rows)
  {
    if (rows.empty())
    {
      return;
    }

    SqliteConnector conn(osw_path, SqliteConnector::SqlOpenMode::READWRITE);
    sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
    sqlite3_stmt* stmt = nullptr;
    conn.executeStatement("DROP TABLE IF EXISTS " + table_name + ";");
    conn.executeStatement(
      "CREATE TABLE " + table_name + " ("
      "FEATURE_ID INTEGER NOT NULL, "
      "SCORE REAL NOT NULL, "
      "RANK INTEGER NOT NULL, "
      "PVALUE REAL, "
      "QVALUE REAL NOT NULL, "
      "PEP REAL NOT NULL);");
    conn.executeStatement("BEGIN IMMEDIATE TRANSACTION;");
    try
    {
      if (table_name == "SCORE_MS2")
      {
        conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_ms2_feature_id ON SCORE_MS2 (FEATURE_ID);");
      }
      else if (table_name == "SCORE_MS1")
      {
        conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_ms1_feature_id ON SCORE_MS1 (FEATURE_ID);");
      }

      sqlitePrepareStatement_(db, &stmt,
        "INSERT INTO " + table_name + " (FEATURE_ID, SCORE, RANK, PVALUE, QVALUE, PEP) VALUES (?, ?, ?, ?, ?, ?);");

      for (const auto& row : rows)
      {
        sqliteCheckResult_(sqlite3_bind_int64(stmt, 1, row.feature_id), db, SQLITE_OK, "sqlite3_bind_int64");
        sqliteCheckResult_(sqlite3_bind_double(stmt, 2, row.score), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_bind_int(stmt, 3, row.rank), db, SQLITE_OK, "sqlite3_bind_int");
        if (row.pvalue.has_value())
        {
          sqliteCheckResult_(sqlite3_bind_double(stmt, 4, *row.pvalue), db, SQLITE_OK, "sqlite3_bind_double");
        }
        else
        {
          sqliteCheckResult_(sqlite3_bind_null(stmt, 4), db, SQLITE_OK, "sqlite3_bind_null");
        }
        sqliteCheckResult_(sqlite3_bind_double(stmt, 5, row.qvalue), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_bind_double(stmt, 6, row.pep), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_step(stmt), db, SQLITE_DONE, "sqlite3_step");
        sqliteCheckResult_(sqlite3_reset(stmt), db, SQLITE_OK, "sqlite3_reset");
        sqliteCheckResult_(sqlite3_clear_bindings(stmt), db, SQLITE_OK, "sqlite3_clear_bindings");
      }

      sqlite3_finalize(stmt);
      stmt = nullptr;
      conn.executeStatement("COMMIT;");
    }
    catch (...)
    {
      if (stmt != nullptr)
      {
        sqlite3_finalize(stmt);
      }
      try
      {
        conn.executeStatement("ROLLBACK;");
      }
      catch (...)
      {
      }
      throw;
    }
  }

  static void writeTransitionScoreTable_(const std::string& osw_path,
                                         const std::vector<TransitionScoreInsertRow_>& rows)
  {
    if (rows.empty())
    {
      return;
    }

    SqliteConnector conn(osw_path, SqliteConnector::SqlOpenMode::READWRITE);
    sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
    sqlite3_stmt* stmt = nullptr;
    conn.executeStatement("DROP TABLE IF EXISTS SCORE_TRANSITION;");
    conn.executeStatement(
      "CREATE TABLE SCORE_TRANSITION ("
      "FEATURE_ID INTEGER NOT NULL, "
      "TRANSITION_ID INTEGER NOT NULL, "
      "SCORE REAL NOT NULL, "
      "RANK INTEGER NOT NULL, "
      "PVALUE REAL, "
      "QVALUE REAL NOT NULL, "
      "PEP REAL NOT NULL);");
    conn.executeStatement("BEGIN IMMEDIATE TRANSACTION;");
    try
    {
      conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_transition_feature_id ON SCORE_TRANSITION (FEATURE_ID);");
      conn.executeStatement("CREATE INDEX IF NOT EXISTS idx_score_transition_transition_id ON SCORE_TRANSITION (TRANSITION_ID);");
      sqlitePrepareStatement_(db, &stmt,
        "INSERT INTO SCORE_TRANSITION (FEATURE_ID, TRANSITION_ID, SCORE, RANK, PVALUE, QVALUE, PEP) VALUES (?, ?, ?, ?, ?, ?, ?);");

      for (const auto& row : rows)
      {
        sqliteCheckResult_(sqlite3_bind_int64(stmt, 1, row.feature_id), db, SQLITE_OK, "sqlite3_bind_int64");
        sqliteCheckResult_(sqlite3_bind_int64(stmt, 2, row.transition_id), db, SQLITE_OK, "sqlite3_bind_int64");
        sqliteCheckResult_(sqlite3_bind_double(stmt, 3, row.score), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_bind_int(stmt, 4, row.rank), db, SQLITE_OK, "sqlite3_bind_int");
        if (row.pvalue.has_value())
        {
          sqliteCheckResult_(sqlite3_bind_double(stmt, 5, *row.pvalue), db, SQLITE_OK, "sqlite3_bind_double");
        }
        else
        {
          sqliteCheckResult_(sqlite3_bind_null(stmt, 5), db, SQLITE_OK, "sqlite3_bind_null");
        }
        sqliteCheckResult_(sqlite3_bind_double(stmt, 6, row.qvalue), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_bind_double(stmt, 7, row.pep), db, SQLITE_OK, "sqlite3_bind_double");
        sqliteCheckResult_(sqlite3_step(stmt), db, SQLITE_DONE, "sqlite3_step");
        sqliteCheckResult_(sqlite3_reset(stmt), db, SQLITE_OK, "sqlite3_reset");
        sqliteCheckResult_(sqlite3_clear_bindings(stmt), db, SQLITE_OK, "sqlite3_clear_bindings");
      }

      sqlite3_finalize(stmt);
      stmt = nullptr;
      conn.executeStatement("COMMIT;");
    }
    catch (...)
    {
      if (stmt != nullptr)
      {
        sqlite3_finalize(stmt);
      }
      try
      {
        conn.executeStatement("ROLLBACK;");
      }
      catch (...)
      {
      }
      throw;
    }
  }

  void buildOSWPQSQLiteBridge_(const OSWPQWorkspace& workspace,
                               const std::string& bridge_osw,
                               const std::string& prepared_library_pqp,
                               const bool enable_uis_scoring)
  {
    OPENMS_LOG_INFO << "Building temporary SQLite bridge for OSWPQ inference and export." << std::endl;
    removeExistingPath_(bridge_osw);
    if (!File::copy(prepared_library_pqp, bridge_osw))
    {
      throw Exception::FileNotReadable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, prepared_library_pqp);
    }

    OpenSwathOSWWriter osw_writer(bridge_osw, enable_uis_scoring);
    osw_writer.writeHeader();

    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_array = ParquetFile::getOptionalColumn(runs_table, "filename");

    std::vector<FeatureScoreInsertRow_> ms1_scores;
    std::vector<FeatureScoreInsertRow_> ms2_scores;
    std::vector<TransitionScoreInsertRow_> transition_scores;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
      osw_writer.addRun(run_id, ParquetFile::getString(filename_array, run_row));

      const std::string run_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id);
      const std::string features_path = run_path + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
      const auto feature_run_id_col = ParquetFile::getColumn(features_table, "run_id");
      const auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
      std::unordered_map<std::string, std::shared_ptr<arrow::Array>> feature_columns;
      const auto add_feature_column = [&](const std::string& name)
      {
        feature_columns.emplace(name, ParquetFile::getOptionalColumn(features_table, name));
      };
      for (const auto& name : std::array<std::string, 8>{
             "exp_rt", "exp_im", "norm_rt", "delta_rt", "left_width", "right_width", "exp_im_leftwidth", "exp_im_rightwidth"})
      {
        add_feature_column(name);
      }
      for (const auto& name : featureMS1ParquetFields_()) add_feature_column(name);
      for (const auto& name : featureMS2ParquetFields_()) add_feature_column(name);
      for (const auto& name : std::array<std::string, 10>{
             "score_ms1_score", "score_ms1_peak_group_rank", "score_ms1_pvalue", "score_ms1_qvalue", "score_ms1_pep",
             "score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_pvalue", "score_ms2_qvalue", "score_ms2_pep"})
      {
        add_feature_column(name);
      }

      OpenSwathOSWWriter::OSWData osw_rows;
      const std::string feature_precursor_path = run_path + "/feature_precursor.parquet";
      const std::string feature_transition_path = run_path + "/feature_transition.parquet";
      Size feature_transition_row_count = 0;
      if (File::exists(feature_transition_path))
      {
        feature_transition_row_count = static_cast<Size>(ParquetFile::readTable(feature_transition_path)->num_rows());
      }
      osw_rows.reserve(static_cast<Size>(features_table->num_rows()), feature_transition_row_count);

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        OpenSwathOSWWriter::FeatureRow feature_row{};
        feature_row[0] = parquetInt64Value_(feature_id_col, row);
        feature_row[1] = parquetInt64Value_(feature_run_id_col, row);
        feature_row[2] = parquetInt64Value_(precursor_id_col, row);
        feature_row[3] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "exp_rt"), row);
        feature_row[4] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "exp_im"), row);
        feature_row[5] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "norm_rt"), row);
        feature_row[6] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "delta_rt"), row);
        feature_row[7] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "left_width"), row);
        feature_row[8] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "right_width"), row);
        feature_row[9] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "exp_im_leftwidth"), row);
        feature_row[10] = parquetDoubleValue_(getOptionalParquetColumn_(feature_columns, "exp_im_rightwidth"), row);
        osw_rows.feature_rows.push_back(std::move(feature_row));

        OpenSwathOSWWriter::FeatureMS1Row feature_ms1_row{};
        feature_ms1_row[0] = parquetInt64Value_(feature_id_col, row);
        for (Size column = 0; column < featureMS1ParquetFields_().size(); ++column)
        {
          feature_ms1_row[column + 1] = parquetDoubleValue_(
            getOptionalParquetColumn_(feature_columns, featureMS1ParquetFields_()[column]), row);
        }
        osw_rows.feature_ms1_rows.push_back(std::move(feature_ms1_row));

        OpenSwathOSWWriter::FeatureMS2Row feature_ms2_row{};
        feature_ms2_row[0] = parquetInt64Value_(feature_id_col, row);
        for (Size column = 0; column < featureMS2ParquetFields_().size(); ++column)
        {
          feature_ms2_row[column + 1] = parquetDoubleValue_(
            getOptionalParquetColumn_(feature_columns, featureMS2ParquetFields_()[column]), row);
        }
        osw_rows.feature_ms2_rows.push_back(std::move(feature_ms2_row));

        const auto score_ms1_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_score"), row);
        const auto score_ms1_rank = parquetOptionalInt64_(getOptionalParquetColumn_(feature_columns, "score_ms1_peak_group_rank"), row);
        const auto score_ms1_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_qvalue"), row);
        const auto score_ms1_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_pep"), row);
        if (score_ms1_score.has_value() && score_ms1_rank.has_value() && score_ms1_qvalue.has_value() && score_ms1_pep.has_value())
        {
          ms1_scores.push_back({
            ParquetFile::getInt64(feature_id_col, row, 0, false),
            *score_ms1_score,
            static_cast<Int32>(*score_ms1_rank),
            parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_pvalue"), row),
            *score_ms1_qvalue,
            *score_ms1_pep
          });
        }

        const auto score_ms2_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_score"), row);
        const auto score_ms2_rank = parquetOptionalInt64_(getOptionalParquetColumn_(feature_columns, "score_ms2_peak_group_rank"), row);
        const auto score_ms2_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_qvalue"), row);
        const auto score_ms2_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_pep"), row);
        if (score_ms2_score.has_value() && score_ms2_rank.has_value() && score_ms2_qvalue.has_value() && score_ms2_pep.has_value())
        {
          ms2_scores.push_back({
            ParquetFile::getInt64(feature_id_col, row, 0, false),
            *score_ms2_score,
            static_cast<Int32>(*score_ms2_rank),
            parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_pvalue"), row),
            *score_ms2_qvalue,
            *score_ms2_pep
          });
        }
      }

      if (File::exists(feature_precursor_path))
      {
        auto feature_precursor_table = ParquetFile::readTable(feature_precursor_path);
        const auto fp_feature_id_col = ParquetFile::getColumn(feature_precursor_table, "feature_id");
        const auto fp_isotope_col = ParquetFile::getOptionalColumn(feature_precursor_table, "precursor_isotope");
        const auto fp_area_col = ParquetFile::getOptionalColumn(feature_precursor_table, "precursor_area_intensity");
        const auto fp_apex_col = ParquetFile::getOptionalColumn(feature_precursor_table, "precursor_apex_intensity");
        for (int64_t row = 0; row < feature_precursor_table->num_rows(); ++row)
        {
          OpenSwathOSWWriter::FeaturePrecursorRow feature_precursor_row{};
          feature_precursor_row[0] = parquetInt64Value_(fp_feature_id_col, row);
          feature_precursor_row[1] = parquetInt64Value_(fp_isotope_col, row);
          feature_precursor_row[2] = parquetDoubleValue_(fp_area_col, row);
          feature_precursor_row[3] = parquetDoubleValue_(fp_apex_col, row);
          osw_rows.feature_precursor_rows.push_back(std::move(feature_precursor_row));
        }
      }

      if (File::exists(feature_transition_path))
      {
        auto feature_transition_table = ParquetFile::readTable(feature_transition_path);
        const auto ft_feature_id_col = ParquetFile::getColumn(feature_transition_table, "feature_id");
        const auto ft_transition_id_col = ParquetFile::getColumn(feature_transition_table, "transition_id");
        std::unordered_map<std::string, std::shared_ptr<arrow::Array>> ft_columns;
        const auto add_ft_column = [&](const std::string& name)
        {
          ft_columns.emplace(name, ParquetFile::getOptionalColumn(feature_transition_table, name));
        };
        for (const auto& name : featureTransitionParquetFields_()) add_ft_column(name);
        for (const auto& name : std::array<std::string, 5>{
               "score_transition_score", "score_transition_rank", "score_transition_pvalue", "score_transition_qvalue", "score_transition_pep"})
        {
          add_ft_column(name);
        }

        for (int64_t row = 0; row < feature_transition_table->num_rows(); ++row)
        {
          OpenSwathOSWWriter::FeatureTransitionRow feature_transition_row{};
          feature_transition_row[0] = parquetInt64Value_(ft_feature_id_col, row);
          feature_transition_row[1] = parquetInt64Value_(ft_transition_id_col, row);
          for (Size column = 0; column < featureTransitionParquetFields_().size(); ++column)
          {
            feature_transition_row[column + 2] = parquetDoubleValue_(
              getOptionalParquetColumn_(ft_columns, featureTransitionParquetFields_()[column]), row);
          }
          osw_rows.feature_transition_rows.push_back(std::move(feature_transition_row));

          const auto transition_score = parquetOptionalDouble_(getOptionalParquetColumn_(ft_columns, "score_transition_score"), row);
          const auto transition_rank = parquetOptionalInt64_(getOptionalParquetColumn_(ft_columns, "score_transition_rank"), row);
          const auto transition_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(ft_columns, "score_transition_qvalue"), row);
          const auto transition_pep = parquetOptionalDouble_(getOptionalParquetColumn_(ft_columns, "score_transition_pep"), row);
          if (transition_score.has_value() && transition_rank.has_value() && transition_qvalue.has_value() && transition_pep.has_value())
          {
            transition_scores.push_back({
              ParquetFile::getInt64(ft_feature_id_col, row, 0, false),
              ParquetFile::getInt64(ft_transition_id_col, row, 0, false),
              *transition_score,
              static_cast<Int32>(*transition_rank),
              parquetOptionalDouble_(getOptionalParquetColumn_(ft_columns, "score_transition_pvalue"), row),
              *transition_qvalue,
              *transition_pep
            });
          }
        }
      }

      osw_writer.writeRows(osw_rows);
    }

    writeFeatureScoreTable_(bridge_osw, "SCORE_MS1", ms1_scores);
    writeFeatureScoreTable_(bridge_osw, "SCORE_MS2", ms2_scores);
    writeTransitionScoreTable_(bridge_osw, transition_scores);
  }

  void syncInferenceScoresToOSWPQ_(const std::string& bridge_osw, OSWPQWorkspace& workspace)
  {
    const auto tasks = getInferenceTasks_();
    if (tasks.empty())
    {
      return;
    }

    const auto [include_ipf_peptide_id, score_members] = getInferenceScoreColumns_(tasks);
    if (!include_ipf_peptide_id && score_members.empty())
    {
      return;
    }

    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");
    std::unordered_set<std::string> replace_columns;
    replace_columns.reserve(score_members.size() + 1);
    if (include_ipf_peptide_id)
    {
      replace_columns.insert("ipf_peptide_id");
    }
    for (const auto& score_member : score_members)
    {
      replace_columns.insert(score_member.name);
    }

    SqliteConnector conn(bridge_osw, SqliteConnector::SqlOpenMode::READ_ONLY);
    sqlite3* db = Internal::SqliteHelper::getNativeHandle(conn);
    ProgressLogger bridge_read_progress;
    bridge_read_progress.setLogType(ProgressLogger::CMD);
    const auto query_steps = std::max<int64_t>(1, runs_table->num_rows() * static_cast<int64_t>(tasks.size()));
    bridge_read_progress.startProgress(0, query_steps, "reading feature scores from temporary SQLite bridge for OSWPQ sync");
    int64_t completed_query_steps = 0;

    ProgressLogger progress_logger;
    progress_logger.setLogType(ProgressLogger::CMD);
    progress_logger.startProgress(0, runs_table->num_rows(), "syncing inference scores back to workflow.oswpq");

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
      std::unordered_map<Int64, OpenSwathFeatureScoreRow> row_by_feature;
      OPENMS_LOG_INFO << "Reading inference scores from the temporary SQLite bridge for run_id="
                      << run_id << "." << std::endl;

      for (const auto& task : tasks)
      {
        const std::string run_filter = " FEATURE.RUN_ID = " + StringUtils::toStr(run_id) + " ";
        Size loaded_rows = 0;
        if (task.level == InferenceLevel::Peptidoform)
        {
          const std::string statement =
            "SELECT FEATURE.ID, SCORE_IPF_BEST.PEPTIDE_ID, SCORE_IPF_BEST.PRECURSOR_PEAKGROUP_PEP, SCORE_IPF_BEST.PEP, SCORE_IPF_BEST.QVALUE "
            "FROM FEATURE "
            "INNER JOIN ("
            "  SELECT FEATURE_ID, PEPTIDE_ID, PRECURSOR_PEAKGROUP_PEP, PEP, QVALUE "
            "  FROM ("
            "    SELECT FEATURE_ID, PEPTIDE_ID, PRECURSOR_PEAKGROUP_PEP, PEP, QVALUE, "
            "           ROW_NUMBER() OVER (PARTITION BY FEATURE_ID ORDER BY PEP, PEPTIDE_ID) AS RN "
            "    FROM SCORE_IPF"
            "  ) "
            "  WHERE RN = 1"
            ") AS SCORE_IPF_BEST ON SCORE_IPF_BEST.FEATURE_ID = FEATURE.ID "
            "WHERE " + run_filter +
            "ORDER BY FEATURE.ID;";
          loaded_rows = readRunInferenceRows_(db, statement, run_id, row_by_feature,
            [&](sqlite3_stmt* stmt, OpenSwathFeatureScoreRow& row)
            {
              if (sqlite3_column_type(stmt, 1) != SQLITE_NULL) row.ipf_peptide_id = sqlite3_column_int64(stmt, 1);
              if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) row.score_ipf_precursor_peakgroup_pep = sqlite3_column_double(stmt, 2);
              if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) row.score_ipf_pep = sqlite3_column_double(stmt, 3);
              if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) row.score_ipf_qvalue = sqlite3_column_double(stmt, 4);
            });
        }
        else if (task.level == InferenceLevel::Peptide && task.context.has_value())
        {
          const std::string context = inferenceContextValue_(*task.context);
          std::string statement =
            "SELECT FEATURE.ID, SCORE_PEPTIDE.SCORE, SCORE_PEPTIDE.PVALUE, SCORE_PEPTIDE.QVALUE, SCORE_PEPTIDE.PEP "
            "FROM FEATURE "
            "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON FEATURE.PRECURSOR_ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "
            "INNER JOIN SCORE_PEPTIDE ON SCORE_PEPTIDE.PEPTIDE_ID = PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID ";
          if (*task.context == InferenceContext::Global)
          {
            statement += "AND SCORE_PEPTIDE.CONTEXT = '" + context + "' ";
          }
          else
          {
            statement += "AND SCORE_PEPTIDE.RUN_ID = FEATURE.RUN_ID "
                         "AND SCORE_PEPTIDE.CONTEXT = '" + context + "' ";
          }
          statement += "WHERE " + run_filter + "ORDER BY FEATURE.ID;";
          loaded_rows = readRunInferenceRows_(db, statement, run_id, row_by_feature,
            [&](sqlite3_stmt* stmt, OpenSwathFeatureScoreRow& row)
            {
              OptionalDoubleMember score_member = nullptr;
              OptionalDoubleMember pvalue_member = nullptr;
              OptionalDoubleMember qvalue_member = nullptr;
              OptionalDoubleMember pep_member = nullptr;
              switch (*task.context)
              {
                case InferenceContext::Global:
                  score_member = &OpenSwathFeatureScoreRow::score_peptide_global_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_peptide_global_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_peptide_global_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_peptide_global_pep;
                  break;
                case InferenceContext::ExperimentWide:
                  score_member = &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pep;
                  break;
                case InferenceContext::RunSpecific:
                  score_member = &OpenSwathFeatureScoreRow::score_peptide_run_specific_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_peptide_run_specific_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_peptide_run_specific_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_peptide_run_specific_pep;
                  break;
              }
              if (sqlite3_column_type(stmt, 1) != SQLITE_NULL) row.*score_member = sqlite3_column_double(stmt, 1);
              if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) row.*pvalue_member = sqlite3_column_double(stmt, 2);
              if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) row.*qvalue_member = sqlite3_column_double(stmt, 3);
              if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) row.*pep_member = sqlite3_column_double(stmt, 4);
            });
        }
        else if (task.level == InferenceLevel::Protein && task.context.has_value())
        {
          const std::string context = inferenceContextValue_(*task.context);
          std::string statement =
            "WITH UNIQUE_PEPTIDE_PROTEIN_MAPPING AS ("
            "  SELECT PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID AS PEPTIDE_ID, PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID "
            "  FROM ("
            "    SELECT PEPTIDE_ID, COUNT(*) AS NUM_PROTEINS "
            "    FROM PEPTIDE_PROTEIN_MAPPING "
            "    GROUP BY PEPTIDE_ID"
            "  ) AS PROTEINS_PER_PEPTIDE "
            "  INNER JOIN PEPTIDE_PROTEIN_MAPPING ON PROTEINS_PER_PEPTIDE.PEPTIDE_ID = PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID "
            "  WHERE NUM_PROTEINS = 1"
            ") "
            "SELECT FEATURE.ID, SCORE_PROTEIN.SCORE, SCORE_PROTEIN.PVALUE, SCORE_PROTEIN.QVALUE, SCORE_PROTEIN.PEP "
            "FROM FEATURE "
            "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON FEATURE.PRECURSOR_ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "
            "INNER JOIN UNIQUE_PEPTIDE_PROTEIN_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = UNIQUE_PEPTIDE_PROTEIN_MAPPING.PEPTIDE_ID "
            "INNER JOIN SCORE_PROTEIN ON SCORE_PROTEIN.PROTEIN_ID = UNIQUE_PEPTIDE_PROTEIN_MAPPING.PROTEIN_ID ";
          if (*task.context == InferenceContext::Global)
          {
            statement += "AND SCORE_PROTEIN.CONTEXT = '" + context + "' ";
          }
          else
          {
            statement += "AND SCORE_PROTEIN.RUN_ID = FEATURE.RUN_ID "
                         "AND SCORE_PROTEIN.CONTEXT = '" + context + "' ";
          }
          statement += "WHERE " + run_filter + "ORDER BY FEATURE.ID;";
          loaded_rows = readRunInferenceRows_(db, statement, run_id, row_by_feature,
            [&](sqlite3_stmt* stmt, OpenSwathFeatureScoreRow& row)
            {
              OptionalDoubleMember score_member = nullptr;
              OptionalDoubleMember pvalue_member = nullptr;
              OptionalDoubleMember qvalue_member = nullptr;
              OptionalDoubleMember pep_member = nullptr;
              switch (*task.context)
              {
                case InferenceContext::Global:
                  score_member = &OpenSwathFeatureScoreRow::score_protein_global_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_protein_global_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_protein_global_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_protein_global_pep;
                  break;
                case InferenceContext::ExperimentWide:
                  score_member = &OpenSwathFeatureScoreRow::score_protein_experiment_wide_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_protein_experiment_wide_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pep;
                  break;
                case InferenceContext::RunSpecific:
                  score_member = &OpenSwathFeatureScoreRow::score_protein_run_specific_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_protein_run_specific_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_protein_run_specific_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_protein_run_specific_pep;
                  break;
              }
              if (sqlite3_column_type(stmt, 1) != SQLITE_NULL) row.*score_member = sqlite3_column_double(stmt, 1);
              if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) row.*pvalue_member = sqlite3_column_double(stmt, 2);
              if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) row.*qvalue_member = sqlite3_column_double(stmt, 3);
              if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) row.*pep_member = sqlite3_column_double(stmt, 4);
            });
        }
        else if (task.level == InferenceLevel::Gene && task.context.has_value())
        {
          const std::string context = inferenceContextValue_(*task.context);
          std::string statement =
            "WITH UNIQUE_PEPTIDE_GENE_MAPPING AS ("
            "  SELECT PEPTIDE_GENE_MAPPING.PEPTIDE_ID AS PEPTIDE_ID, PEPTIDE_GENE_MAPPING.GENE_ID "
            "  FROM ("
            "    SELECT PEPTIDE_ID, COUNT(*) AS NUM_GENES "
            "    FROM PEPTIDE_GENE_MAPPING "
            "    GROUP BY PEPTIDE_ID"
            "  ) AS GENES_PER_PEPTIDE "
            "  INNER JOIN PEPTIDE_GENE_MAPPING ON GENES_PER_PEPTIDE.PEPTIDE_ID = PEPTIDE_GENE_MAPPING.PEPTIDE_ID "
            "  WHERE NUM_GENES = 1"
            ") "
            "SELECT FEATURE.ID, SCORE_GENE.SCORE, SCORE_GENE.PVALUE, SCORE_GENE.QVALUE, SCORE_GENE.PEP "
            "FROM FEATURE "
            "INNER JOIN PRECURSOR_PEPTIDE_MAPPING ON FEATURE.PRECURSOR_ID = PRECURSOR_PEPTIDE_MAPPING.PRECURSOR_ID "
            "INNER JOIN UNIQUE_PEPTIDE_GENE_MAPPING ON PRECURSOR_PEPTIDE_MAPPING.PEPTIDE_ID = UNIQUE_PEPTIDE_GENE_MAPPING.PEPTIDE_ID "
            "INNER JOIN SCORE_GENE ON SCORE_GENE.GENE_ID = UNIQUE_PEPTIDE_GENE_MAPPING.GENE_ID ";
          if (*task.context == InferenceContext::Global)
          {
            statement += "AND SCORE_GENE.CONTEXT = '" + context + "' ";
          }
          else
          {
            statement += "AND SCORE_GENE.RUN_ID = FEATURE.RUN_ID "
                         "AND SCORE_GENE.CONTEXT = '" + context + "' ";
          }
          statement += "WHERE " + run_filter + "ORDER BY FEATURE.ID;";
          loaded_rows = readRunInferenceRows_(db, statement, run_id, row_by_feature,
            [&](sqlite3_stmt* stmt, OpenSwathFeatureScoreRow& row)
            {
              OptionalDoubleMember score_member = nullptr;
              OptionalDoubleMember pvalue_member = nullptr;
              OptionalDoubleMember qvalue_member = nullptr;
              OptionalDoubleMember pep_member = nullptr;
              switch (*task.context)
              {
                case InferenceContext::Global:
                  score_member = &OpenSwathFeatureScoreRow::score_gene_global_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_gene_global_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_gene_global_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_gene_global_pep;
                  break;
                case InferenceContext::ExperimentWide:
                  score_member = &OpenSwathFeatureScoreRow::score_gene_experiment_wide_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_gene_experiment_wide_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pep;
                  break;
                case InferenceContext::RunSpecific:
                  score_member = &OpenSwathFeatureScoreRow::score_gene_run_specific_score;
                  pvalue_member = &OpenSwathFeatureScoreRow::score_gene_run_specific_pvalue;
                  qvalue_member = &OpenSwathFeatureScoreRow::score_gene_run_specific_qvalue;
                  pep_member = &OpenSwathFeatureScoreRow::score_gene_run_specific_pep;
                  break;
              }
              if (sqlite3_column_type(stmt, 1) != SQLITE_NULL) row.*score_member = sqlite3_column_double(stmt, 1);
              if (sqlite3_column_type(stmt, 2) != SQLITE_NULL) row.*pvalue_member = sqlite3_column_double(stmt, 2);
              if (sqlite3_column_type(stmt, 3) != SQLITE_NULL) row.*qvalue_member = sqlite3_column_double(stmt, 3);
              if (sqlite3_column_type(stmt, 4) != SQLITE_NULL) row.*pep_member = sqlite3_column_double(stmt, 4);
            });
        }

        ++completed_query_steps;
        bridge_read_progress.setProgress(completed_query_steps);
        OPENMS_LOG_INFO << "Loaded " << loaded_rows << " feature-score rows for "
                        << inferenceTaskLabel_(task) << " sync in run_id=" << run_id << "." << std::endl;
      }

      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
      OPENMS_LOG_INFO << "Updating OSWPQ feature parquet for run_id=" << run_id
                      << " (" << features_table->num_rows() << " rows)." << std::endl;

      arrow::Int64Builder ipf_peptide_id_builder;
      std::vector<std::unique_ptr<arrow::DoubleBuilder>> double_builders;
      double_builders.reserve(score_members.size());
      for (Size i = 0; i < score_members.size(); ++i)
      {
        double_builders.push_back(std::make_unique<arrow::DoubleBuilder>());
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        const Int64 feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
        const OpenSwathFeatureScoreRow* score_row = nullptr;
        const auto feature_it = row_by_feature.find(feature_id);
        if (feature_it != row_by_feature.end())
        {
          score_row = &feature_it->second;
        }

        if (include_ipf_peptide_id && score_row != nullptr && score_row->ipf_peptide_id.has_value())
        {
          ParquetFile::appendOrThrow(ipf_peptide_id_builder.Append(*score_row->ipf_peptide_id), "ipf_peptide_id");
        }
        else if (include_ipf_peptide_id)
        {
          ParquetFile::appendOrThrow(ipf_peptide_id_builder.AppendNull(), "ipf_peptide_id");
        }

        for (Size column = 0; column < score_members.size(); ++column)
        {
          const auto member_value = (score_row != nullptr) ? score_row->*(score_members[column].member) : std::optional<double>{};
          if (member_value.has_value())
          {
            ParquetFile::appendOrThrow(double_builders[column]->Append(*member_value), score_members[column].name);
          }
          else
          {
            ParquetFile::appendOrThrow(double_builders[column]->AppendNull(), score_members[column].name);
          }
        }
      }

      std::vector<std::shared_ptr<arrow::Field>> extra_fields;
      std::vector<std::shared_ptr<arrow::Array>> extra_arrays;
      extra_fields.reserve(score_members.size() + (include_ipf_peptide_id ? 1 : 0));
      extra_arrays.reserve(score_members.size() + (include_ipf_peptide_id ? 1 : 0));
      if (include_ipf_peptide_id)
      {
        extra_fields.push_back(arrow::field("ipf_peptide_id", arrow::int64(), true));
        extra_arrays.push_back(ParquetFile::finishArray(ipf_peptide_id_builder, "ipf_peptide_id"));
      }
      for (Size column = 0; column < score_members.size(); ++column)
      {
        extra_fields.push_back(arrow::field(score_members[column].name, arrow::float64(), true));
        extra_arrays.push_back(ParquetFile::finishArray(*double_builders[column], score_members[column].name));
      }

      replaceParquetColumns_(features_path, replace_columns, extra_fields, extra_arrays);
      progress_logger.setProgress(run_row + 1);
    }

    bridge_read_progress.endProgress();
    progress_logger.endProgress();
    OPENMS_LOG_INFO << "Finished syncing inference scores back to workflow.oswpq." << std::endl;
    workspace.dirty = true;
  }

  std::vector<LevelContextInputRow> buildOSWPQLevelContextInputRows_(const OSWPQWorkspace& workspace,
                                                                     const PreparedLibraryLookup_& lookup,
                                                                     const InferenceLevel level,
                                                                     const InferenceContext context) const
  {
    if (level == InferenceLevel::Peptidoform)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ level-context inference does not support peptidoform rows.");
    }

    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    std::map<std::pair<Int64, Int64>, LevelContextInputRow> best_rows;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto precursor_id_array = ParquetFile::getColumn(features_table, "precursor_id");
      const auto score_ms2_array = ParquetFile::getOptionalColumn(features_table, "score_ms2_score");
      if (score_ms2_array == nullptr)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Level-context inference on OSWPQ requires score_ms2_score in features.parquet.");
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        if (score_ms2_array->IsNull(row))
        {
          continue;
        }

        const Int64 precursor_id = ParquetFile::getInt64(precursor_id_array, row, 0, false);
        const auto precursor_it = lookup.precursors.find(precursor_id);
        if (precursor_it == lookup.precursors.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Missing prepared-library precursor metadata for precursor_id=" + StringUtils::toStr(precursor_id));
        }

        const double score = ParquetFile::getDouble(score_ms2_array, row, 0.0, false);
        const auto entity_ids = mappedEntitiesForFeature_(lookup, level, precursor_id);
        for (const Int64 entity_id : entity_ids)
        {
          const Int64 run_key = context == InferenceContext::Global ? std::numeric_limits<Int64>::min() : run_id;
          const auto key = std::make_pair(run_key, entity_id);
          LevelContextInputRow candidate;
          if (context != InferenceContext::Global)
          {
            candidate.run_id = run_id;
          }
          candidate.group_id = context == InferenceContext::Global ?
            StringUtils::toStr(entity_id) :
            StringUtils::toStr(run_id) + "_" + StringUtils::toStr(entity_id);
          candidate.entity_id = entity_id;
          candidate.decoy = precursor_it->second.decoy;
          candidate.score = score;
          candidate.context = context;

          const auto existing = best_rows.find(key);
          if (existing == best_rows.end() || candidate.score > existing->second.score)
          {
            best_rows[key] = std::move(candidate);
          }
        }
      }
    }

    std::vector<LevelContextInputRow> rows;
    rows.reserve(best_rows.size());
    for (auto& [key, row] : best_rows)
    {
      rows.push_back(std::move(row));
    }

    OPENMS_LOG_INFO << "Read " << rows.size() << " best-score rows for "
                    << toString(level) << " inference in '" << toString(context)
                    << "' context." << std::endl;
    return rows;
  }

  void applyLevelContextResultsToOSWPQ_(OSWPQWorkspace& workspace,
                                        const PreparedLibraryLookup_& lookup,
                                        const std::map<InferenceLevel, std::vector<LevelContextResultRow>>& results_by_level) const
  {
    const auto tasks = getInferenceTasks_();
    if (tasks.empty())
    {
      return;
    }

    const auto [include_ipf_peptide_id, score_members] = getInferenceScoreColumns_(tasks);
    if (include_ipf_peptide_id)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ score sync currently supports peptide, protein, and gene inference only.");
    }
    if (score_members.empty())
    {
      return;
    }

    std::map<InferenceLevel, LevelContextResultMaps_> result_maps;
    for (const auto& [level, results] : results_by_level)
    {
      result_maps[level] = buildLevelContextResultMaps_(results);
    }

    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");

    std::unordered_set<std::string> replace_columns;
    replace_columns.reserve(score_members.size());
    for (const auto& score_member : score_members)
    {
      replace_columns.insert(score_member.name);
    }

    ProgressLogger progress_logger;
    progress_logger.setLogType(ProgressLogger::CMD);
    progress_logger.startProgress(0, runs_table->num_rows(), "syncing level-context scores back to workflow.oswpq");

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto precursor_id_array = ParquetFile::getColumn(features_table, "precursor_id");

      std::vector<std::unique_ptr<arrow::DoubleBuilder>> double_builders;
      double_builders.reserve(score_members.size());
      for (Size i = 0; i < score_members.size(); ++i)
      {
        double_builders.push_back(std::make_unique<arrow::DoubleBuilder>());
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        const Int64 precursor_id = ParquetFile::getInt64(precursor_id_array, row, 0, false);
        OpenSwathFeatureScoreRow score_row;

        const auto assign_level = [&](const InferenceLevel level,
                                      const InferenceContext context,
                                      OptionalDoubleMember score_member,
                                      OptionalDoubleMember pvalue_member,
                                      OptionalDoubleMember qvalue_member,
                                      OptionalDoubleMember pep_member)
        {
          const auto level_it = result_maps.find(level);
          if (level_it == result_maps.end())
          {
            return;
          }
          const auto entity_ids = mappedEntitiesForFeature_(lookup, level, precursor_id);
          const auto best_result = selectBestLevelContextResult_(entity_ids, run_id, level_it->second, context);
          if (!best_result.has_value())
          {
            return;
          }
          score_row.*score_member = best_result->score;
          score_row.*pvalue_member = best_result->pvalue;
          score_row.*qvalue_member = best_result->qvalue;
          score_row.*pep_member = best_result->pep;
        };

        assign_level(InferenceLevel::Peptide, InferenceContext::Global,
                     &OpenSwathFeatureScoreRow::score_peptide_global_score,
                     &OpenSwathFeatureScoreRow::score_peptide_global_pvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_global_pep);
        assign_level(InferenceLevel::Peptide, InferenceContext::ExperimentWide,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pep);
        assign_level(InferenceLevel::Peptide, InferenceContext::RunSpecific,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_pvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_pep);

        assign_level(InferenceLevel::Protein, InferenceContext::Global,
                     &OpenSwathFeatureScoreRow::score_protein_global_score,
                     &OpenSwathFeatureScoreRow::score_protein_global_pvalue,
                     &OpenSwathFeatureScoreRow::score_protein_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_global_pep);
        assign_level(InferenceLevel::Protein, InferenceContext::ExperimentWide,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pvalue,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pep);
        assign_level(InferenceLevel::Protein, InferenceContext::RunSpecific,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_pvalue,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_pep);

        assign_level(InferenceLevel::Gene, InferenceContext::Global,
                     &OpenSwathFeatureScoreRow::score_gene_global_score,
                     &OpenSwathFeatureScoreRow::score_gene_global_pvalue,
                     &OpenSwathFeatureScoreRow::score_gene_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_global_pep);
        assign_level(InferenceLevel::Gene, InferenceContext::ExperimentWide,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pvalue,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pep);
        assign_level(InferenceLevel::Gene, InferenceContext::RunSpecific,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_pvalue,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_pep);

        for (Size column = 0; column < score_members.size(); ++column)
        {
          const auto value = score_row.*(score_members[column].member);
          if (value.has_value())
          {
            ParquetFile::appendOrThrow(double_builders[column]->Append(*value), score_members[column].name);
          }
          else
          {
            ParquetFile::appendOrThrow(double_builders[column]->AppendNull(), score_members[column].name);
          }
        }
      }

      std::vector<std::shared_ptr<arrow::Field>> extra_fields;
      std::vector<std::shared_ptr<arrow::Array>> extra_arrays;
      extra_fields.reserve(score_members.size());
      extra_arrays.reserve(score_members.size());
      for (Size column = 0; column < score_members.size(); ++column)
      {
        extra_fields.push_back(arrow::field(score_members[column].name, arrow::float64(), true));
        extra_arrays.push_back(ParquetFile::finishArray(*double_builders[column], score_members[column].name));
      }

      replaceParquetColumns_(features_path, replace_columns, extra_fields, extra_arrays);
      progress_logger.setProgress(run_row + 1);
    }

    progress_logger.endProgress();
    workspace.dirty = true;
  }

  void runInferenceOSWPQ_(OSWPQWorkspace& workspace, const PreparedLibraryLookup_& lookup)
  {
    const auto tasks = getInferenceTasks_();
    if (tasks.empty())
    {
      return;
    }

    if (std::any_of(tasks.begin(), tasks.end(),
                    [](const auto& task) { return task.level == InferenceLevel::Peptidoform; }))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ inference currently supports peptide, protein, and gene inference only. Use workflow:working_format osw for peptidoform/IPF inference.");
    }

    const ErrorEstimationConfig error_config = getErrorConfig_();
    const bool has_run_specific_task = std::any_of(tasks.begin(), tasks.end(),
      [](const auto& task)
      {
        return task.context.has_value() && *task.context == InferenceContext::RunSpecific;
      });
    const std::map<Int64, std::string> run_basenames = has_run_specific_task ? readOSWPQRunBasenames_(workspace) : std::map<Int64, std::string>{};

    std::map<InferenceLevel, std::vector<LevelContextResultRow>> results_by_level;
    for (const auto& task : tasks)
    {
      ProgressLogger progress_logger;
      progress_logger.setLogType(ProgressLogger::CMD);
      progress_logger.startProgress(0, 1, inferenceTaskLabel_(task));

      LevelContextInferenceConfig config;
      config.level = task.level;
      config.context = *task.context;
      config.error = error_config;
      const auto input_rows = buildOSWPQLevelContextInputRows_(workspace, lookup, task.level, config.context);

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

      auto& level_results = results_by_level[task.level];
      level_results.insert(level_results.end(), results.begin(), results.end());
      logLevelContextSummary_(input_rows, results, task.level, config.context, run_basenames);
      progress_logger.endProgress();
    }

    for (const auto& [level, results] : results_by_level)
    {
      writeLevelContextResultsParquet_(workspace, level, results);
    }
    applyLevelContextResultsToOSWPQ_(workspace, lookup, results_by_level);
    workspace.dirty = true;
  }

  static ExportQValueMaps_ buildPeptideQValueMaps_(const std::vector<LevelContextResultRow>& results)
  {
    ExportQValueMaps_ maps;
    for (const auto& row : results)
    {
      switch (row.context)
      {
        case InferenceContext::Global:
          maps.global[row.entity_id] = row.qvalue;
          break;
        case InferenceContext::ExperimentWide:
          if (row.run_id.has_value()) maps.experiment_wide[{*row.run_id, row.entity_id}] = row.qvalue;
          break;
        case InferenceContext::RunSpecific:
          if (row.run_id.has_value()) maps.run_specific[{*row.run_id, row.entity_id}] = row.qvalue;
          break;
      }
    }
    return maps;
  }

  static ExportQValueMaps_ buildAggregatedEntityQValueMaps_(const std::vector<LevelContextResultRow>& results,
                                                            const std::unordered_map<Int64, std::vector<Int64>>& peptide_to_entities)
  {
    std::unordered_map<Int64, std::vector<Int64>> entity_to_peptides;
    for (const auto& [peptide_id, entity_ids] : peptide_to_entities)
    {
      for (const Int64 entity_id : entity_ids)
      {
        entity_to_peptides[entity_id].push_back(peptide_id);
      }
    }

    ExportQValueMaps_ maps;
    auto update_min = [](auto& target, const auto& key, const double qvalue)
    {
      const auto existing = target.find(key);
      if (existing == target.end() || qvalue < existing->second)
      {
        target[key] = qvalue;
      }
    };

    for (const auto& row : results)
    {
      const auto peptide_it = entity_to_peptides.find(row.entity_id);
      if (peptide_it == entity_to_peptides.end())
      {
        continue;
      }
      for (const Int64 peptide_id : peptide_it->second)
      {
        switch (row.context)
        {
          case InferenceContext::Global:
            update_min(maps.global, peptide_id, row.qvalue);
            break;
          case InferenceContext::ExperimentWide:
            if (row.run_id.has_value()) update_min(maps.experiment_wide, std::make_pair(*row.run_id, peptide_id), row.qvalue);
            break;
          case InferenceContext::RunSpecific:
            if (row.run_id.has_value()) update_min(maps.run_specific, std::make_pair(*row.run_id, peptide_id), row.qvalue);
            break;
        }
      }
    }
    return maps;
  }

  std::unordered_map<Int64, TransitionAggregation_> buildTransitionAggregations_(const OSWPQWorkspace& workspace,
                                                                                  const PreparedLibraryLookup_& lookup,
                                                                                  const double max_transition_pep) const
  {
    std::unordered_map<Int64, TransitionAggregation_> aggregations;
    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string feature_transition_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet";
      if (!File::exists(feature_transition_path))
      {
        continue;
      }

      auto table = ParquetFile::readTable(feature_transition_path);
      const auto feature_id_col = ParquetFile::getColumn(table, "feature_id");
      const auto transition_id_col = ParquetFile::getColumn(table, "transition_id");
      const auto area_col = ParquetFile::getOptionalColumn(table, "area_intensity");
      const auto apex_col = ParquetFile::getOptionalColumn(table, "apex_intensity");
      const auto pep_col = ParquetFile::getOptionalColumn(table, "score_transition_pep");
      const bool filter_by_transition_pep = pep_col != nullptr;

      for (int64_t row = 0; row < table->num_rows(); ++row)
      {
        if (filter_by_transition_pep)
        {
          if (pep_col->IsNull(row))
          {
            continue;
          }
          if (ParquetFile::getDouble(pep_col, row, 1.0, false) >= max_transition_pep)
          {
            continue;
          }
        }

        const Int64 transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
        const auto transition_it = lookup.transitions.find(transition_id);
        if (transition_it == lookup.transitions.end() || transition_it->second.decoy)
        {
          continue;
        }

        const Int64 feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
        auto& aggregation = aggregations[feature_id];
        aggregation.areas.push_back(StringUtils::toStr(ParquetFile::getDouble(area_col, row, 0.0, true)));
        aggregation.apices.push_back(StringUtils::toStr(ParquetFile::getDouble(apex_col, row, 0.0, true)));
        aggregation.annotations.push_back(
          StringUtils::toStr(transition_id) + "_" + transition_it->second.type +
          StringUtils::toStr(transition_it->second.ordinal) + "_" + StringUtils::toStr(transition_it->second.charge));
      }
    }

    return aggregations;
  }

  OpenSwathFeatureScoreTable readOSWPQFeatureScoreTable_(const OSWPQWorkspace& workspace,
                                                         const PreparedLibraryLookup_& lookup,
                                                         const OpenSwathParquetExportConfig& config) const
  {
    OpenSwathFeatureScoreTable table;
    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    bool discovered_dynamic_columns = false;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string filename = (filename_col != nullptr && !filename_col->IsNull(run_row)) ? ParquetFile::getString(filename_col, run_row) : "";
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
      const auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");

      std::unordered_map<std::string, std::shared_ptr<arrow::Array>> feature_columns;
      const auto add_feature_column = [&](const std::string& name)
      {
        feature_columns.emplace(name, ParquetFile::getOptionalColumn(features_table, name));
      };
      for (const auto& name : std::array<std::string, 27>{
             "exp_rt", "exp_im", "norm_rt", "delta_rt", "left_width", "right_width", "exp_im_leftwidth", "exp_im_rightwidth",
             "score_ms1_score", "score_ms1_peak_group_rank", "score_ms1_pvalue", "score_ms1_qvalue", "score_ms1_pep",
             "score_ms2_score", "score_ms2_peak_group_rank", "score_ms2_pvalue", "score_ms2_qvalue", "score_ms2_pep",
             "ipf_peptide_id", "score_ipf_precursor_peakgroup_pep", "score_ipf_pep", "score_ipf_qvalue",
             "score_peptide_global_score", "score_peptide_global_pvalue", "score_peptide_global_qvalue", "score_peptide_global_pep",
             "score_peptide_experiment_wide_score"})
      {
        add_feature_column(name);
      }
      for (const auto& name : std::array<std::string, 20>{
             "score_peptide_experiment_wide_pvalue", "score_peptide_experiment_wide_qvalue", "score_peptide_experiment_wide_pep",
             "score_peptide_run_specific_score", "score_peptide_run_specific_pvalue", "score_peptide_run_specific_qvalue", "score_peptide_run_specific_pep",
             "score_protein_global_score", "score_protein_global_pvalue", "score_protein_global_qvalue", "score_protein_global_pep",
             "score_protein_experiment_wide_score", "score_protein_experiment_wide_pvalue", "score_protein_experiment_wide_qvalue", "score_protein_experiment_wide_pep",
             "score_protein_run_specific_score", "score_protein_run_specific_pvalue", "score_protein_run_specific_qvalue", "score_protein_run_specific_pep",
             "score_gene_global_score"})
      {
        add_feature_column(name);
      }
      for (const auto& name : std::array<std::string, 11>{
             "score_gene_global_pvalue", "score_gene_global_qvalue", "score_gene_global_pep",
             "score_gene_experiment_wide_score", "score_gene_experiment_wide_pvalue", "score_gene_experiment_wide_qvalue", "score_gene_experiment_wide_pep",
             "score_gene_run_specific_score", "score_gene_run_specific_pvalue", "score_gene_run_specific_qvalue", "score_gene_run_specific_pep"})
      {
        add_feature_column(name);
      }

      if (!discovered_dynamic_columns)
      {
        for (const auto& name : featureMS1ParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(features_table, name) != nullptr)
          {
            table.feature_ms1_column_names.emplace_back(name);
          }
        }
        for (const auto& name : featureMS2ParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(features_table, name) != nullptr)
          {
            table.feature_ms2_column_names.emplace_back(name);
          }
        }
        discovered_dynamic_columns = true;
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        const Int64 precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
        const auto precursor_it = lookup.precursors.find(precursor_id);
        if (precursor_it == lookup.precursors.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Missing prepared-library precursor metadata for precursor_id=" + StringUtils::toStr(precursor_id));
        }
        if (config.filters.exclude_decoys && precursor_it->second.decoy)
        {
          continue;
        }

        const auto peptide_mapping_it = lookup.precursor_to_peptides.find(precursor_id);
        if (peptide_mapping_it == lookup.precursor_to_peptides.end())
        {
          continue;
        }

        for (const Int64 peptide_id : peptide_mapping_it->second)
        {
          const auto peptide_it = lookup.peptides.find(peptide_id);
          if (peptide_it == lookup.peptides.end())
          {
            continue;
          }

          std::vector<std::optional<Int64>> protein_ids{std::nullopt};
          const auto protein_mapping_it = lookup.peptide_to_proteins.find(peptide_id);
          if (protein_mapping_it != lookup.peptide_to_proteins.end() && !protein_mapping_it->second.empty())
          {
            protein_ids.clear();
            protein_ids.reserve(protein_mapping_it->second.size());
            for (const Int64 protein_id : protein_mapping_it->second)
            {
              protein_ids.push_back(protein_id);
            }
          }

          std::vector<std::optional<Int64>> gene_ids{std::nullopt};
          const auto gene_mapping_it = lookup.peptide_to_genes.find(peptide_id);
          if (gene_mapping_it != lookup.peptide_to_genes.end() && !gene_mapping_it->second.empty())
          {
            gene_ids.clear();
            gene_ids.reserve(gene_mapping_it->second.size());
            for (const Int64 gene_id : gene_mapping_it->second)
            {
              gene_ids.push_back(gene_id);
            }
          }

          for (const auto protein_id : protein_ids)
          {
            for (const auto gene_id : gene_ids)
            {
              OpenSwathFeatureScoreRow score_row;
              score_row.protein_id = protein_id.value_or(-1);
              score_row.peptide_id = peptide_id;
              score_row.ipf_peptide_id = parquetOptionalInt64_(getOptionalParquetColumn_(feature_columns, "ipf_peptide_id"), row);
              score_row.precursor_id = precursor_id;
              score_row.unmodified_sequence = peptide_it->second.unmodified_sequence;
              score_row.modified_sequence = peptide_it->second.modified_sequence;
              score_row.precursor_traml_id = precursor_it->second.traml_id;
              score_row.precursor_group_label = precursor_it->second.group_label;
              score_row.precursor_mz = precursor_it->second.precursor_mz;
              score_row.precursor_charge = precursor_it->second.charge;
              score_row.precursor_library_intensity = precursor_it->second.library_intensity;
              score_row.precursor_library_rt = precursor_it->second.library_rt;
              score_row.precursor_library_drift_time = precursor_it->second.library_drift_time;
              score_row.peptide_decoy = peptide_it->second.decoy;
              score_row.precursor_decoy = precursor_it->second.decoy;
              score_row.run_id = run_id;
              score_row.filename = filename;
              score_row.feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
              score_row.exp_rt = ParquetFile::getDouble(getOptionalParquetColumn_(feature_columns, "exp_rt"), row, 0.0, true);
              score_row.exp_im = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "exp_im"), row);
              score_row.norm_rt = ParquetFile::getDouble(getOptionalParquetColumn_(feature_columns, "norm_rt"), row, 0.0, true);
              score_row.delta_rt = ParquetFile::getDouble(getOptionalParquetColumn_(feature_columns, "delta_rt"), row, 0.0, true);
              score_row.left_width = ParquetFile::getDouble(getOptionalParquetColumn_(feature_columns, "left_width"), row, 0.0, true);
              score_row.right_width = ParquetFile::getDouble(getOptionalParquetColumn_(feature_columns, "right_width"), row, 0.0, true);
              score_row.im_left_width = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "exp_im_leftwidth"), row);
              score_row.im_right_width = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "exp_im_rightwidth"), row);

              for (const auto& name : table.feature_ms1_column_names)
              {
                score_row.feature_ms1_values.push_back(ParquetFile::getDouble(ParquetFile::getOptionalColumn(features_table, name), row, 0.0, true));
              }
              for (const auto& name : table.feature_ms2_column_names)
              {
                score_row.feature_ms2_values.push_back(ParquetFile::getDouble(ParquetFile::getOptionalColumn(features_table, name), row, 0.0, true));
              }

              score_row.score_ms1_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_score"), row);
              const auto score_ms1_rank = parquetOptionalInt64_(getOptionalParquetColumn_(feature_columns, "score_ms1_peak_group_rank"), row);
              if (score_ms1_rank.has_value()) score_row.score_ms1_rank = static_cast<Int32>(*score_ms1_rank);
              score_row.score_ms1_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_pvalue"), row);
              score_row.score_ms1_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_qvalue"), row);
              score_row.score_ms1_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms1_pep"), row);
              score_row.score_ms2_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_score"), row);
              const auto score_ms2_rank = parquetOptionalInt64_(getOptionalParquetColumn_(feature_columns, "score_ms2_peak_group_rank"), row);
              if (score_ms2_rank.has_value()) score_row.score_ms2_peak_group_rank = static_cast<Int32>(*score_ms2_rank);
              score_row.score_ms2_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_pvalue"), row);
              score_row.score_ms2_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_qvalue"), row);
              score_row.score_ms2_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ms2_pep"), row);
              score_row.score_ipf_precursor_peakgroup_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ipf_precursor_peakgroup_pep"), row);
              score_row.score_ipf_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ipf_pep"), row);
              score_row.score_ipf_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_ipf_qvalue"), row);
              score_row.score_peptide_global_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_global_score"), row);
              score_row.score_peptide_global_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_global_pvalue"), row);
              score_row.score_peptide_global_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_global_qvalue"), row);
              score_row.score_peptide_global_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_global_pep"), row);
              score_row.score_peptide_experiment_wide_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_score"), row);
              score_row.score_peptide_experiment_wide_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_pvalue"), row);
              score_row.score_peptide_experiment_wide_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_qvalue"), row);
              score_row.score_peptide_experiment_wide_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_pep"), row);
              score_row.score_peptide_run_specific_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_score"), row);
              score_row.score_peptide_run_specific_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_pvalue"), row);
              score_row.score_peptide_run_specific_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_qvalue"), row);
              score_row.score_peptide_run_specific_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_pep"), row);
              score_row.score_protein_global_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_global_score"), row);
              score_row.score_protein_global_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_global_pvalue"), row);
              score_row.score_protein_global_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_global_qvalue"), row);
              score_row.score_protein_global_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_global_pep"), row);
              score_row.score_protein_experiment_wide_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_score"), row);
              score_row.score_protein_experiment_wide_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_pvalue"), row);
              score_row.score_protein_experiment_wide_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_qvalue"), row);
              score_row.score_protein_experiment_wide_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_pep"), row);
              score_row.score_protein_run_specific_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_score"), row);
              score_row.score_protein_run_specific_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_pvalue"), row);
              score_row.score_protein_run_specific_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_qvalue"), row);
              score_row.score_protein_run_specific_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_pep"), row);
              score_row.score_gene_global_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_global_score"), row);
              score_row.score_gene_global_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_global_pvalue"), row);
              score_row.score_gene_global_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_global_qvalue"), row);
              score_row.score_gene_global_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_global_pep"), row);
              score_row.score_gene_experiment_wide_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_score"), row);
              score_row.score_gene_experiment_wide_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_pvalue"), row);
              score_row.score_gene_experiment_wide_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_qvalue"), row);
              score_row.score_gene_experiment_wide_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_pep"), row);
              score_row.score_gene_run_specific_score = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_score"), row);
              score_row.score_gene_run_specific_pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_pvalue"), row);
              score_row.score_gene_run_specific_qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_qvalue"), row);
              score_row.score_gene_run_specific_pep = parquetOptionalDouble_(getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_pep"), row);

              if (protein_id.has_value())
              {
                const auto protein_it = lookup.proteins.find(*protein_id);
                if (protein_it != lookup.proteins.end())
                {
                  score_row.protein_accession = protein_it->second.accession;
                  score_row.protein_decoy = protein_it->second.decoy;
                }
              }
              if (gene_id.has_value())
              {
                const auto gene_it = lookup.genes.find(*gene_id);
                if (gene_it != lookup.genes.end())
                {
                  score_row.gene_id = *gene_id;
                  score_row.gene_name = gene_it->second.name;
                  score_row.gene_decoy = gene_it->second.decoy;
                }
              }

              table.rows.push_back(std::move(score_row));
            }
          }
        }
      }
    }

    std::stable_sort(table.rows.begin(), table.rows.end(),
      [](const OpenSwathFeatureScoreRow& lhs, const OpenSwathFeatureScoreRow& rhs)
      {
        if (lhs.precursor_id != rhs.precursor_id) return lhs.precursor_id < rhs.precursor_id;
        if (lhs.feature_id != rhs.feature_id) return lhs.feature_id < rhs.feature_id;
        if (lhs.peptide_id != rhs.peptide_id) return lhs.peptide_id < rhs.peptide_id;
        if (lhs.protein_id != rhs.protein_id) return lhs.protein_id < rhs.protein_id;
        return lhs.gene_id.value_or(-1) < rhs.gene_id.value_or(-1);
      });

    OPENMS_LOG_INFO << "Read " << table.rows.size() << " precursor feature score rows." << std::endl;
    return table;
  }

  OpenSwathTransitionScoreTable readOSWPQTransitionScoreTable_(const OSWPQWorkspace& workspace,
                                                               const PreparedLibraryLookup_& lookup,
                                                               const OpenSwathParquetExportConfig& config) const
  {
    OpenSwathTransitionScoreTable table;
    if (!config.include_transition_data)
    {
      return table;
    }

    std::unordered_map<Int64, std::vector<FeatureTransitionObservation_>> observations_by_transition;
    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    bool discovered_dynamic_columns = false;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string feature_transition_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet";
      if (!File::exists(feature_transition_path))
      {
        continue;
      }

      auto feature_transition_table = ParquetFile::readTable(feature_transition_path);
      if (!discovered_dynamic_columns)
      {
        for (const auto& name : featureTransitionParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(feature_transition_table, name) != nullptr)
          {
            table.feature_transition_column_names.emplace_back(name);
          }
        }
        discovered_dynamic_columns = true;
      }

      const auto feature_id_col = ParquetFile::getColumn(feature_transition_table, "feature_id");
      const auto transition_id_col = ParquetFile::getColumn(feature_transition_table, "transition_id");
      std::unordered_map<std::string, std::shared_ptr<arrow::Array>> transition_columns;
      const auto add_transition_column = [&](const std::string& name)
      {
        transition_columns.emplace(name, ParquetFile::getOptionalColumn(feature_transition_table, name));
      };
      for (const auto& name : table.feature_transition_column_names)
      {
        add_transition_column(name);
      }
      for (const auto& name : std::array<std::string, 5>{
             "score_transition_score", "score_transition_rank", "score_transition_pvalue", "score_transition_qvalue", "score_transition_pep"})
      {
        add_transition_column(name);
      }

      for (int64_t row = 0; row < feature_transition_table->num_rows(); ++row)
      {
        FeatureTransitionObservation_ observation;
        observation.run_id = run_id;
        observation.feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
        observation.values.reserve(table.feature_transition_column_names.size());
        for (const auto& name : table.feature_transition_column_names)
        {
          const auto column = getOptionalParquetColumn_(transition_columns, name);
          if (column != nullptr && !column->IsNull(row))
          {
            observation.values.push_back(ParquetFile::getDouble(column, row, 0.0, false));
          }
          else
          {
            observation.values.push_back(std::nullopt);
          }
        }
        observation.score = parquetOptionalDouble_(getOptionalParquetColumn_(transition_columns, "score_transition_score"), row);
        const auto rank = parquetOptionalInt64_(getOptionalParquetColumn_(transition_columns, "score_transition_rank"), row);
        if (rank.has_value()) observation.rank = static_cast<Int32>(*rank);
        observation.pvalue = parquetOptionalDouble_(getOptionalParquetColumn_(transition_columns, "score_transition_pvalue"), row);
        observation.qvalue = parquetOptionalDouble_(getOptionalParquetColumn_(transition_columns, "score_transition_qvalue"), row);
        observation.pep = parquetOptionalDouble_(getOptionalParquetColumn_(transition_columns, "score_transition_pep"), row);
        const Int64 transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
        observations_by_transition[transition_id].push_back(std::move(observation));
      }
    }

    for (const auto& [transition_id, transition] : lookup.transitions)
    {
      if (config.filters.exclude_decoys && transition.decoy)
      {
        continue;
      }

      std::vector<std::optional<Int64>> peptide_ids{std::nullopt};
      if (!transition.peptide_ids.empty())
      {
        peptide_ids.clear();
        peptide_ids.reserve(transition.peptide_ids.size());
        for (const Int64 peptide_id : transition.peptide_ids)
        {
          peptide_ids.push_back(peptide_id);
        }
      }

      std::vector<Int64> precursor_ids = transition.precursor_ids;
      if (precursor_ids.empty())
      {
        precursor_ids.push_back(-1);
      }

      std::vector<FeatureTransitionObservation_> empty_observations(1);
      const auto obs_it = observations_by_transition.find(transition_id);
      const auto& observations = obs_it != observations_by_transition.end() ? obs_it->second : empty_observations;

      for (const Int64 precursor_id : precursor_ids)
      {
        for (const auto peptide_id : peptide_ids)
        {
          for (const auto& observation : observations)
          {
            OpenSwathTransitionScoreRow score_row;
            score_row.run_id = observation.run_id;
            score_row.ipf_peptide_id = peptide_id;
            score_row.precursor_id = precursor_id;
            score_row.transition_id = transition_id;
            score_row.transition_traml_id = transition.traml_id;
            score_row.product_mz = transition.product_mz;
            score_row.transition_charge = transition.charge;
            score_row.transition_type = transition.type;
            score_row.transition_ordinal = transition.ordinal;
            score_row.annotation = transition.annotation;
            score_row.transition_detecting = transition.detecting;
            score_row.transition_library_intensity = transition.library_intensity;
            score_row.transition_decoy = transition.decoy;
            score_row.feature_id = observation.feature_id;
            score_row.feature_transition_values = observation.values;
            score_row.score_transition_score = observation.score;
            score_row.score_transition_rank = observation.rank;
            score_row.score_transition_pvalue = observation.pvalue;
            score_row.score_transition_qvalue = observation.qvalue;
            score_row.score_transition_pep = observation.pep;
            table.rows.push_back(std::move(score_row));
          }
        }
      }
    }

    std::stable_sort(table.rows.begin(), table.rows.end(),
      [](const OpenSwathTransitionScoreRow& lhs, const OpenSwathTransitionScoreRow& rhs)
      {
        if (lhs.precursor_id != rhs.precursor_id) return lhs.precursor_id < rhs.precursor_id;
        if (lhs.transition_id != rhs.transition_id) return lhs.transition_id < rhs.transition_id;
        return lhs.feature_id.value_or(std::numeric_limits<Int64>::max()) < rhs.feature_id.value_or(std::numeric_limits<Int64>::max());
      });

    OPENMS_LOG_INFO << "Read " << table.rows.size() << " transition score rows." << std::endl;
    return table;
  }

  static bool oswpqHasIPFColumns_(const OSWPQWorkspace& workspace)
  {
    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      auto feature_table = ParquetFile::readTable(workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet");
      if (ParquetFile::getOptionalColumn(feature_table, "score_ipf_qvalue") != nullptr ||
          ParquetFile::getOptionalColumn(feature_table, "score_ipf_pep") != nullptr)
      {
        return true;
      }
    }
    return false;
  }

  std::vector<OpenSwathExportRow> readOSWPQExportRows_(const OSWPQWorkspace& workspace,
                                                       const PreparedLibraryLookup_& lookup,
                                                       const OpenSwathExportFilterConfig& config) const
  {
    if (config.use_alignment)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ export does not support alignment recovery because no alignment parquet tables are written yet.");
    }
    if (config.ipf_mode != OpenSwathIPFExportMode::Disable && oswpqHasIPFColumns_(workspace))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ export currently supports standard OpenSWATH-style exports only. Use Export:*:ipf disable or workflow:working_format osw for IPF-aware export.");
    }

    const auto peptide_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Peptide);
    const auto protein_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Protein);
    const auto gene_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Gene);
    const ExportQValueMaps_ peptide_qvalues = buildPeptideQValueMaps_(peptide_results);
    const ExportQValueMaps_ protein_qvalues = buildAggregatedEntityQValueMaps_(protein_results, lookup.peptide_to_proteins);
    const ExportQValueMaps_ gene_qvalues = buildAggregatedEntityQValueMaps_(gene_results, lookup.peptide_to_genes);
    const auto transition_aggregations = config.transition_quantification ?
      buildTransitionAggregations_(workspace, lookup, config.max_transition_pep) :
      std::unordered_map<Int64, TransitionAggregation_>{};

    auto runs_table = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet");
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    std::vector<OpenSwathExportRow> rows;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string filename = (filename_col != nullptr && !filename_col->IsNull(run_row)) ? ParquetFile::getString(filename_col, run_row) : "";
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = ParquetFile::readTable(features_path);
      const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
      const auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");
      const auto exp_rt_col = ParquetFile::getOptionalColumn(features_table, "exp_rt");
      const auto norm_rt_col = ParquetFile::getOptionalColumn(features_table, "norm_rt");
      const auto delta_rt_col = ParquetFile::getOptionalColumn(features_table, "delta_rt");
      const auto left_width_col = ParquetFile::getOptionalColumn(features_table, "left_width");
      const auto right_width_col = ParquetFile::getOptionalColumn(features_table, "right_width");
      const auto exp_im_col = ParquetFile::getOptionalColumn(features_table, "exp_im");
      const auto exp_im_left_col = ParquetFile::getOptionalColumn(features_table, "exp_im_leftwidth");
      const auto exp_im_right_col = ParquetFile::getOptionalColumn(features_table, "exp_im_rightwidth");
      const auto ms2_area_col = ParquetFile::getOptionalColumn(features_table, "ms2_area_intensity");
      const auto ms1_area_col = ParquetFile::getOptionalColumn(features_table, "ms1_area_intensity");
      const auto ms1_apex_col = ParquetFile::getOptionalColumn(features_table, "ms1_apex_intensity");
      const auto score_ms1_pep_col = ParquetFile::getOptionalColumn(features_table, "score_ms1_pep");
      const auto score_ms2_score_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_score");
      const auto score_ms2_qvalue_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_qvalue");
      const auto score_ms2_pep_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_pep");
      const auto score_ms2_rank_col = ParquetFile::getOptionalColumn(features_table, "score_ms2_peak_group_rank");
      if (score_ms2_qvalue_col == nullptr || score_ms2_score_col == nullptr)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Direct OSWPQ export requires score_ms2_score and score_ms2_qvalue in features.parquet.");
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        if (score_ms2_qvalue_col->IsNull(row))
        {
          continue;
        }
        const double ms2_qvalue = ParquetFile::getDouble(score_ms2_qvalue_col, row, 1.0, false);
        if (!(ms2_qvalue < config.max_rs_peakgroup_qvalue))
        {
          continue;
        }

        const Int64 precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
        const auto precursor_it = lookup.precursors.find(precursor_id);
        if (precursor_it == lookup.precursors.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Missing prepared-library precursor metadata for precursor_id=" + StringUtils::toStr(precursor_id));
        }
        const auto peptide_mapping_it = lookup.precursor_to_peptides.find(precursor_id);
        if (peptide_mapping_it == lookup.precursor_to_peptides.end())
        {
          continue;
        }

        for (const Int64 peptide_id : peptide_mapping_it->second)
        {
          const auto peptide_it = lookup.peptides.find(peptide_id);
          if (peptide_it == lookup.peptides.end())
          {
            continue;
          }

          OpenSwathExportRow export_row;
          export_row.run_id = run_id;
          export_row.filename = filename;
          export_row.run_name = File::stemName(filename);
          if (export_row.run_name.empty())
          {
            export_row.run_name = File::basename(filename);
          }
          if (export_row.run_name.empty())
          {
            export_row.run_name = "RUN_ID " + StringUtils::toStr(run_id);
          }

          export_row.feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
          export_row.peptide_id = peptide_id;
          export_row.precursor_id = precursor_id;
          export_row.transition_group_id = StringUtils::toStr(precursor_id);
          export_row.decoy = precursor_it->second.decoy;
          export_row.sequence = peptide_it->second.unmodified_sequence;
          export_row.full_peptide_name = peptide_it->second.modified_sequence;
          export_row.protein_name = lookup.protein_names_by_peptide.count(peptide_id) ? lookup.protein_names_by_peptide.at(peptide_id) : "";
          export_row.gene_name = lookup.gene_names_by_peptide.count(peptide_id) ? lookup.gene_names_by_peptide.at(peptide_id) : "";
          export_row.charge = precursor_it->second.charge;
          export_row.mz = precursor_it->second.precursor_mz;
          export_row.rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);
          export_row.assay_rt = export_row.rt - ParquetFile::getDouble(delta_rt_col, row, 0.0, true);
          export_row.delta_rt = ParquetFile::getDouble(delta_rt_col, row, 0.0, true);
          export_row.irt = ParquetFile::getDouble(norm_rt_col, row, 0.0, true);
          export_row.assay_irt = precursor_it->second.library_rt.value_or(std::numeric_limits<double>::quiet_NaN());
          export_row.delta_irt = export_row.irt - export_row.assay_irt;
          export_row.intensity = ParquetFile::getDouble(ms2_area_col, row, 0.0, true);
          export_row.aggr_prec_peak_area = parquetOptionalDouble_(ms1_area_col, row);
          export_row.aggr_prec_peak_apex = parquetOptionalDouble_(ms1_apex_col, row);
          export_row.left_width = ParquetFile::getDouble(left_width_col, row, 0.0, true);
          export_row.right_width = ParquetFile::getDouble(right_width_col, row, 0.0, true);
          export_row.exp_im = parquetOptionalDouble_(exp_im_col, row);
          export_row.im_left_width = parquetOptionalDouble_(exp_im_left_col, row);
          export_row.im_right_width = parquetOptionalDouble_(exp_im_right_col, row);
          export_row.ms1_pep = parquetOptionalDouble_(score_ms1_pep_col, row);
          export_row.ms2_pep = parquetOptionalDouble_(score_ms2_pep_col, row);
          export_row.peak_group_rank = static_cast<Int32>(parquetOptionalInt64_(score_ms2_rank_col, row).value_or(0));
          export_row.d_score = ParquetFile::getDouble(score_ms2_score_col, row, 0.0, false);
          export_row.m_score = ms2_qvalue;
          export_row.pep = parquetOptionalDouble_(score_ms2_pep_col, row);

          const auto peptide_global_it = peptide_qvalues.global.find(peptide_id);
          if (peptide_global_it != peptide_qvalues.global.end()) export_row.peptide_global_qvalue = peptide_global_it->second;
          const auto peptide_experiment_it = peptide_qvalues.experiment_wide.find({run_id, peptide_id});
          if (peptide_experiment_it != peptide_qvalues.experiment_wide.end()) export_row.peptide_experiment_wide_qvalue = peptide_experiment_it->second;
          const auto peptide_run_it = peptide_qvalues.run_specific.find({run_id, peptide_id});
          if (peptide_run_it != peptide_qvalues.run_specific.end()) export_row.peptide_run_specific_qvalue = peptide_run_it->second;

          const auto protein_global_it = protein_qvalues.global.find(peptide_id);
          if (protein_global_it != protein_qvalues.global.end()) export_row.protein_global_qvalue = protein_global_it->second;
          const auto protein_experiment_it = protein_qvalues.experiment_wide.find({run_id, peptide_id});
          if (protein_experiment_it != protein_qvalues.experiment_wide.end()) export_row.protein_experiment_wide_qvalue = protein_experiment_it->second;
          const auto protein_run_it = protein_qvalues.run_specific.find({run_id, peptide_id});
          if (protein_run_it != protein_qvalues.run_specific.end()) export_row.protein_run_specific_qvalue = protein_run_it->second;

          const auto gene_global_it = gene_qvalues.global.find(peptide_id);
          if (gene_global_it != gene_qvalues.global.end()) export_row.gene_global_qvalue = gene_global_it->second;
          const auto gene_experiment_it = gene_qvalues.experiment_wide.find({run_id, peptide_id});
          if (gene_experiment_it != gene_qvalues.experiment_wide.end()) export_row.gene_experiment_wide_qvalue = gene_experiment_it->second;
          const auto gene_run_it = gene_qvalues.run_specific.find({run_id, peptide_id});
          if (gene_run_it != gene_qvalues.run_specific.end()) export_row.gene_run_specific_qvalue = gene_run_it->second;

          const auto transition_it = transition_aggregations.find(export_row.feature_id);
          if (transition_it != transition_aggregations.end())
          {
            export_row.aggr_peak_area = joinStrings_(transition_it->second.areas);
            export_row.aggr_peak_apex = joinStrings_(transition_it->second.apices);
            export_row.aggr_fragment_annotation = joinStrings_(transition_it->second.annotations);
          }

          rows.push_back(std::move(export_row));
        }
      }
    }

    if (config.exclude_decoys)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(), [](const auto& row) { return row.decoy; }), rows.end());
    }
    if (config.peptide)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(),
                                [&](const auto& row)
                                {
                                  return !row.peptide_global_qvalue.has_value() ||
                                         *row.peptide_global_qvalue >= config.max_global_peptide_qvalue;
                                }),
                 rows.end());
    }
    if (config.protein)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(),
                                [&](const auto& row)
                                {
                                  return !row.protein_global_qvalue.has_value() ||
                                         *row.protein_global_qvalue >= config.max_global_protein_qvalue;
                                }),
                 rows.end());
    }
    if (config.gene)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(),
                                [&](const auto& row)
                                {
                                  return !row.gene_global_qvalue.has_value() ||
                                         *row.gene_global_qvalue >= config.max_global_gene_qvalue;
                                }),
                 rows.end());
    }

    std::stable_sort(rows.begin(), rows.end(),
      [](const OpenSwathExportRow& lhs, const OpenSwathExportRow& rhs)
      {
        if (lhs.precursor_id != rhs.precursor_id) return lhs.precursor_id < rhs.precursor_id;
        if (lhs.feature_id != rhs.feature_id) return lhs.feature_id < rhs.feature_id;
        return lhs.peptide_id < rhs.peptide_id;
      });

    OPENMS_LOG_INFO << "Read " << rows.size() << " filtered export rows." << std::endl;
    return rows;
  }

  void runExportsOSWPQ_(const OSWPQWorkspace& workspace,
                        const PreparedLibraryLookup_& lookup,
                        const StringList& input_files,
                        const std::string& out_dir) const
  {
    const auto tasks = getExportTasks_();
    if (tasks.empty())
    {
      return;
    }

    File::makeDir(out_dir);
    const std::string base_path = makeExportBasePath_(input_files, out_dir);
    std::optional<std::vector<OpenSwathExportRow>> matrix_rows_cache;

    for (const auto& task : tasks)
    {
      ProgressLogger progress_logger;
      progress_logger.setLogType(ProgressLogger::CMD);
      progress_logger.startProgress(0, 1, exportTaskLabel_(task));

      if (task.type == ExportTaskType::Results)
      {
        const auto results_config = getResultsConfig_();
        const auto rows = readOSWPQExportRows_(workspace, lookup, results_config.filters);
        const std::string suffix = results_config.format == OpenSwathExportFileFormat::Parquet ? ".results.parquet" : ".results.tsv";
        OpenSwathResultsExporter::write(base_path + suffix, rows, results_config);
      }
      else if (task.type == ExportTaskType::FeatureParquet)
      {
        const auto parquet_config = getParquetConfig_();
        const auto feature_table = readOSWPQFeatureScoreTable_(workspace, lookup, parquet_config);
        const std::string feature_out = base_path + ".precursor.feature.scores.parquet";
        OpenSwathParquetExporter::writeFeatureScores(feature_out, feature_table);
        if (parquet_config.include_transition_data)
        {
          const auto transition_table = readOSWPQTransitionScoreTable_(workspace, lookup, parquet_config);
          if (!transition_table.rows.empty())
          {
            const std::string transition_out = base_path + ".transition.feature.scores.parquet";
            OpenSwathParquetExporter::writeTransitionScores(transition_out, transition_table);
          }
        }
      }
      else if (task.type == ExportTaskType::Matrix)
      {
        const auto matrix_config = getMatrixConfig_(*task.matrix_level);
        if (!matrix_rows_cache.has_value())
        {
          matrix_rows_cache = readOSWPQExportRows_(workspace, lookup, matrix_config.filters);
        }
        const auto matrix = OpenSwathMatrixExporter::buildMatrix(*matrix_rows_cache, matrix_config);
        const std::string suffix = matrix_config.format == OpenSwathExportFileFormat::Parquet ? ".matrix.parquet" : ".matrix.tsv";
        OpenSwathMatrixExporter::writeMatrix(base_path + "." + toString(*task.matrix_level) + suffix, matrix, matrix_config);
      }

      progress_logger.endProgress();
    }
  }

  ExitCodes runExtractionToWorkflow_(const std::string& prepared_library_pqp,
                                     const std::string& workflow_output,
                                     const WorkflowFormat workflow_format,
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
    if (workflow_format == WorkflowFormat::OSW)
    {
      removeExistingPath_(workflow_output);
      if (!File::copy(prepared_library_pqp, workflow_output))
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, workflow_output);
      }
    }
    else
    {
      removeExistingPath_(workflow_output);
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
    const bool write_osw = workflow_format == WorkflowFormat::OSW;
    const bool write_parquet = workflow_format == WorkflowFormat::OSWPQ;
    OpenSwathOSWWriter oswwriter(write_osw ? workflow_output : "", enable_uis_scoring);
    OpenSwathOSWParquetWriter parquet_writer;
    parquet_writer.setPreserveExisting(true);
    if (write_osw)
    {
      oswwriter.writeHeader();
    }

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
      OpenSwathOSWWriter::OSWData run_osw_rows;
      const bool store_features_in_feature_file = false;
      OpenSwathWorkflow wf(use_ms1_traces, use_ms1_im, prm, pasef, mrm_mode, outer_loop_threads);
      wf.setLogType(log_type_);
      wf.performExtraction(swath_maps, trafo_rtnorm, cp_current, cp_ms1_current, feature_finder_param_run,
                           transition_exp, run_feature_file, store_features_in_feature_file, oswwriter, chromatogramConsumer,
                           batchSize, ms1_isotopes, load_into_memory, mrm_map_param,
                           mobilogramConsumer.get(), innerBatchSize, maxConcurrentSwaths,
                           write_parquet ? &run_osw_rows : nullptr);

      swath_maps.clear();
      if (mobilogramConsumer)
      {
        mobilogramConsumer->finalize();
      }
      delete chromatogramConsumer;
      if (write_parquet)
      {
        parquet_writer.write(workflow_output, transition_exp, run_osw_rows, cur_run, current_run_files[0], enable_uis_scoring);
      }
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
      const WorkflowFormat workflow_format = getWorkflowFormat_();
      const std::string workflow_output = workflow_format == WorkflowFormat::OSWPQ ?
        working_dir.path + "/workflow.oswpq" :
        working_dir.path + "/workflow.osw";

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

      const ExitCodes extraction_result = runExtractionToWorkflow_(prepared_library_pqp, workflow_output, workflow_format, enable_uis_scoring);
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
        rescoring.score(workflow_output, level);
      }

      if (workflow_format == WorkflowFormat::OSW)
      {
        runInference_(workflow_output);
        runExports_(workflow_output, input_files, out_dir);
      }
      else
      {
        OSWPQWorkspace workflow_workspace = prepareOSWPQWorkspace_(workflow_output);
        const auto export_tasks = getExportTasks_();
        const bool needs_transition_metadata =
          std::any_of(export_tasks.begin(), export_tasks.end(),
            [&](const auto& task)
            {
              if (task.type == ExportTaskType::FeatureParquet)
              {
                return true;
              }
              if (task.type == ExportTaskType::Results)
              {
                return getResultsConfig_().filters.transition_quantification;
              }
              if (task.type == ExportTaskType::Matrix)
              {
                return getMatrixConfig_(*task.matrix_level).filters.transition_quantification;
              }
              return false;
            });
        Param tsv_reader_param = getParam_().copy("TargetedDataExtraction:Library:", true);
        OpenSwath::LightTargetedExperiment transition_exp = loadTransitionList(FileTypes::PQP, prepared_library_pqp, tsv_reader_param);
        PreparedLibraryLookup_ prepared_library_lookup =
          buildPreparedLibraryLookupFromLightTargetedExperiment_(transition_exp, needs_transition_metadata);
        try
        {
          runInferenceOSWPQ_(workflow_workspace, prepared_library_lookup);
          runExportsOSWPQ_(workflow_workspace, prepared_library_lookup, input_files, out_dir);
          commitOSWPQWorkspace_(workflow_workspace);
        }
        catch (...)
        {
          try
          {
            commitOSWPQWorkspace_(workflow_workspace);
          }
          catch (const Exception::BaseException& commit_error)
          {
            OPENMS_LOG_WARN << "Failed to preserve updated workflow.oswpq after an error: "
                            << commit_error.what() << std::endl;
          }
          throw;
        }
      }

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
