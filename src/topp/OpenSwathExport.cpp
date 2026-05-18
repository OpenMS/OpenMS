// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathMatrixExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathParquetExporter.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathResultsExporter.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/SYSTEM/File.h>

#include <optional>
#include <vector>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_OpenSwathExport OpenSwathExport

@brief Export scored OpenSWATH OSW files to user-facing TSV / Parquet results tables and quantification matrices.

OpenSwathExport operates on scored OpenSWATH OSW files and can generate multiple
exports in one invocation:
- scored feature Parquet exports (`*.precursor.feature.scores.parquet`, optional `*.transition.feature.scores.parquet`)
- filtered user-facing results exports (`*.results.tsv` or `*.results.parquet`)
- quantification matrices at precursor, peptide, protein, or gene level

The tool currently supports OSW input. OpenSWATH Parquet input (`.oswpq`) is
reserved for future work and is rejected with a clear error for now.

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathExport.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathExport.html
*/

/// @cond TOPPCLASSES
class TOPPOpenSwathExport :
  public TOPPBase
{
public:
  TOPPOpenSwathExport() :
    TOPPBase("OpenSwathExport", "Export scored OpenSWATH results to TSV or Parquet tables.")
  {
  }

protected:
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

  static bool toBool_(const String& value)
  {
    return value == "true";
  }

  static OpenSwathExportFileFormat toExportFormat_(const String& value)
  {
    return value == "parquet" ? OpenSwathExportFileFormat::Parquet : OpenSwathExportFileFormat::TSV;
  }

  static OpenSwathIPFExportMode toIPFMode_(const String& value)
  {
    if (value == "augmented") return OpenSwathIPFExportMode::Augmented;
    if (value == "disable") return OpenSwathIPFExportMode::Disable;
    return OpenSwathIPFExportMode::Peptidoform;
  }

  static OpenSwathMatrixNormalization toNormalization_(const String& value)
  {
    if (value == "median") return OpenSwathMatrixNormalization::Median;
    if (value == "medianmedian") return OpenSwathMatrixNormalization::MedianMedian;
    return OpenSwathMatrixNormalization::None;
  }

  OpenSwathExportFilterConfig getFilterConfig_(const String& prefix) const
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
    config.filters = getFilterConfig_("results");
    config.format = toExportFormat_(getStringOption_("results:format"));
    return config;
  }

  OpenSwathMatrixExportConfig getMatrixConfig_(const OpenSwathMatrixLevel level) const
  {
    OpenSwathMatrixExportConfig config;
    config.filters = getFilterConfig_("matrix");
    config.level = level;
    config.normalization = toNormalization_(getStringOption_("matrix:normalization"));
    const Int top_n = getIntOption_("matrix:top_n");
    if (top_n < 0 || (level != OpenSwathMatrixLevel::Precursor && top_n < 1))
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
        "Parameter 'matrix:top_n' must be >= 1 for peptide/protein/gene matrix export and >= 0 for precursor export.",
        String(top_n));
    }
    config.top_n = static_cast<Size>(top_n);
    config.consistent_top = toBool_(getStringOption_("matrix:consistent_top"));
    config.format = toExportFormat_(getStringOption_("matrix:format"));
    return config;
  }

  OpenSwathParquetExportConfig getParquetConfig_() const
  {
    OpenSwathParquetExportConfig config;
    config.include_transition_data = toBool_(getStringOption_("parquet:include_transition_data"));
    config.filters.exclude_decoys = toBool_(getStringOption_("parquet:exclude_decoys"));
    return config;
  }

  static String makeOutputDir_(const String& input_file, const String& requested_out_dir)
  {
    if (!requested_out_dir.empty())
    {
      return requested_out_dir;
    }
    const String input_dir = File::path(input_file);
    return input_dir.empty() ? "." : input_dir;
  }

  static String makeBasePath_(const String& input_file, const String& out_dir)
  {
    const String dir = makeOutputDir_(input_file, out_dir);
    return dir + "/" + File::stemName(input_file);
  }

  std::vector<ExportTask> getTasks_() const
  {
    std::vector<ExportTask> tasks;
    if (toBool_(getStringOption_("parquet:run")))
    {
      tasks.push_back({ExportTaskType::FeatureParquet, std::nullopt});
    }
    if (toBool_(getStringOption_("results:run")))
    {
      tasks.push_back({ExportTaskType::Results, std::nullopt});
    }
    if (toBool_(getStringOption_("matrix:precursor")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Precursor});
    }
    if (toBool_(getStringOption_("matrix:peptide")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Peptide});
    }
    if (toBool_(getStringOption_("matrix:protein")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Protein});
    }
    if (toBool_(getStringOption_("matrix:gene")))
    {
      tasks.push_back({ExportTaskType::Matrix, OpenSwathMatrixLevel::Gene});
    }
    if (tasks.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No export types were enabled. Set at least one parquet/results/matrix selector to 'true'.");
    }
    return tasks;
  }

  static String taskLabel_(const ExportTask& task)
  {
    switch (task.type)
    {
      case ExportTaskType::FeatureParquet: return "feature parquet export";
      case ExportTaskType::Results: return "results export";
      case ExportTaskType::Matrix: return toString(*task.matrix_level) + " matrix export";
    }
    return "export";
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input scored OpenSWATH result file in OSW format.", true);
    setValidFormats_("in", {"osw", "oswpq"});
    registerOutputDir_("out_dir", "<dir>", "", "Optional output directory. If omitted, files are written next to the input OSW.", false, false);

    registerTOPPSubsection_("parquet", "Options for scored feature parquet export");
    registerStringOption_("parquet:run", "<true|false>", "false", "Whether to export scored feature parquet tables.", false);
    setValidStrings_("parquet:run", {"true", "false"});
    registerStringOption_("parquet:include_transition_data", "<true|false>", "true", "Whether to also export transition-level parquet tables.", false);
    setValidStrings_("parquet:include_transition_data", {"true", "false"});
    registerStringOption_("parquet:exclude_decoys", "<true|false>", "false", "Whether to exclude decoys from parquet feature exports.", false);
    setValidStrings_("parquet:exclude_decoys", {"true", "false"});

    registerTOPPSubsection_("results", "Options for user-facing results export");
    registerStringOption_("results:run", "<true|false>", "true", "Whether to export the filtered results table.", false);
    setValidStrings_("results:run", {"true", "false"});
    registerStringOption_("results:format", "<choice>", "parquet", "Output format for the filtered results table.", false);
    setValidStrings_("results:format", {"tsv", "parquet"});
    registerStringOption_("results:ipf", "<choice>", "peptidoform", "How IPF results should be represented when SCORE_IPF is present.", false);
    setValidStrings_("results:ipf", {"peptidoform", "augmented", "disable"});
    registerStringOption_("results:transition_quantification", "<true|false>", "true", "Whether to aggregate transition-level quantification columns.", false);
    setValidStrings_("results:transition_quantification", {"true", "false"});
    registerDoubleOption_("results:max_transition_pep", "<num>", 0.7, "Maximum transition PEP used for aggregated transition quantification.", false);
    registerDoubleOption_("results:ipf_max_peptidoform_pep", "<num>", 0.4, "Maximum IPF PEP used for peptidoform export filtering.", false);
    registerDoubleOption_("results:max_rs_peakgroup_qvalue", "<num>", 0.05, "Maximum run-specific peakgroup q-value retained.", false);
    registerStringOption_("results:use_peptide_scores", "<true|false>", "true", "Whether to append peptide-level scores if SCORE_PEPTIDE is present and to require passing the peptide-level threshold.", false);
    setValidStrings_("results:use_peptide_scores", {"true", "false"});
    registerDoubleOption_("results:max_global_peptide_qvalue", "<num>", 0.01, "Maximum global peptide-level q-value retained when peptide filtering is enabled.", false);
    registerStringOption_("results:use_protein_scores", "<true|false>", "true", "Whether to append protein-level scores if SCORE_PROTEIN is present and to require passing the protein-level threshold.", false);
    setValidStrings_("results:use_protein_scores", {"true", "false"});
    registerDoubleOption_("results:max_global_protein_qvalue", "<num>", 0.01, "Maximum global protein-level q-value retained when protein filtering is enabled.", false);
    registerStringOption_("results:use_gene_scores", "<true|false>", "false", "Whether to append gene-level scores if SCORE_GENE is present and to require passing the gene-level threshold.", false);
    setValidStrings_("results:use_gene_scores", {"true", "false"});
    registerDoubleOption_("results:max_global_gene_qvalue", "<num>", 0.01, "Maximum global gene-level q-value retained when gene filtering is enabled.", false);
    registerStringOption_("results:use_alignment", "<true|false>", "false", "Whether to recover aligned features when legacy alignment tables are present.", false);
    setValidStrings_("results:use_alignment", {"true", "false"});
    registerDoubleOption_("results:max_alignment_pep", "<num>", 0.7, "Maximum alignment PEP retained when alignment recovery is enabled.", false);
    registerStringOption_("results:exclude_decoys", "<true|false>", "true", "Whether to exclude decoys from the filtered results export.", false);
    setValidStrings_("results:exclude_decoys", {"true", "false"});

    registerTOPPSubsection_("matrix", "Options for quantification matrix exports");
    registerStringOption_("matrix:precursor", "<true|false>", "false", "Whether to export the precursor-level quantification matrix.", false);
    setValidStrings_("matrix:precursor", {"true", "false"});
    registerStringOption_("matrix:peptide", "<true|false>", "true", "Whether to export the peptide-level quantification matrix.", false);
    setValidStrings_("matrix:peptide", {"true", "false"});
    registerStringOption_("matrix:protein", "<true|false>", "true", "Whether to export the protein-level quantification matrix.", false);
    setValidStrings_("matrix:protein", {"true", "false"});
    registerStringOption_("matrix:gene", "<true|false>", "false", "Whether to export the gene-level quantification matrix.", false);
    setValidStrings_("matrix:gene", {"true", "false"});
    registerStringOption_("matrix:format", "<choice>", "tsv", "Output format for quantification matrices.", false);
    setValidStrings_("matrix:format", {"tsv", "parquet"});
    registerStringOption_("matrix:ipf", "<choice>", "peptidoform", "How IPF results should be represented when SCORE_IPF is present.", false);
    setValidStrings_("matrix:ipf", {"peptidoform", "augmented", "disable"});
    registerStringOption_("matrix:transition_quantification", "<true|false>", "true", "Whether to aggregate transition-level quantification columns onto the filtered feature rows.", false);
    setValidStrings_("matrix:transition_quantification", {"true", "false"});
    registerDoubleOption_("matrix:max_transition_pep", "<num>", 0.7, "Maximum transition PEP used for aggregated transition quantification.", false);
    registerDoubleOption_("matrix:ipf_max_peptidoform_pep", "<num>", 0.4, "Maximum IPF PEP used for peptidoform export filtering.", false);
    registerDoubleOption_("matrix:max_rs_peakgroup_qvalue", "<num>", 0.05, "Maximum run-specific peakgroup q-value retained.", false);
    registerStringOption_("matrix:use_peptide_scores", "<true|false>", "true", "Whether peptide-level scores should be required when peptide filtering is enabled.", false);
    setValidStrings_("matrix:use_peptide_scores", {"true", "false"});
    registerDoubleOption_("matrix:max_global_peptide_qvalue", "<num>", 0.01, "Maximum global peptide-level q-value retained when peptide filtering is enabled.", false);
    registerStringOption_("matrix:use_protein_scores", "<true|false>", "true", "Whether protein-level scores should be required when protein filtering is enabled.", false);
    setValidStrings_("matrix:use_protein_scores", {"true", "false"});
    registerDoubleOption_("matrix:max_global_protein_qvalue", "<num>", 0.01, "Maximum global protein-level q-value retained when protein filtering is enabled.", false);
    registerStringOption_("matrix:use_gene_scores", "<true|false>", "false", "Whether gene-level scores should be required when gene filtering is enabled.", false);
    setValidStrings_("matrix:use_gene_scores", {"true", "false"});
    registerDoubleOption_("matrix:max_global_gene_qvalue", "<num>", 0.01, "Maximum global gene-level q-value retained when gene filtering is enabled.", false);
    registerStringOption_("matrix:use_alignment", "<true|false>", "false", "Whether to recover aligned features when legacy alignment tables are present.", false);
    setValidStrings_("matrix:use_alignment", {"true", "false"});
    registerDoubleOption_("matrix:max_alignment_pep", "<num>", 0.7, "Maximum alignment PEP retained when alignment recovery is enabled.", false);
    registerStringOption_("matrix:exclude_decoys", "<true|false>", "true", "Whether to exclude decoys from the filtered matrix input rows.", false);
    setValidStrings_("matrix:exclude_decoys", {"true", "false"});
    registerIntOption_("matrix:top_n", "<num>", 3, "Number of top precursors / peptides retained during matrix summarization.", false);
    registerStringOption_("matrix:consistent_top", "<true|false>", "true", "Whether to use the same top precursors / peptides across all runs.", false);
    setValidStrings_("matrix:consistent_top", {"true", "false"});
    registerStringOption_("matrix:normalization", "<choice>", "none", "Normalization applied to matrix sample columns after summarization.", false);
    setValidStrings_("matrix:normalization", {"none", "median", "medianmedian"});
  }

  ExitCodes main_(int, const char**) override
  {
    const String input_file = getStringOption_("in");
    if (input_file.hasSuffix(".oswpq") || input_file.hasSuffix(".oswpq.zip"))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "OpenSwathExport currently supports OSW input only. .oswpq support is reserved for future work.");
    }

    const String out_dir = makeOutputDir_(input_file, getOutputDirOption("out_dir"));
    File::makeDir(out_dir);
    const String base_path = makeBasePath_(input_file, out_dir);
    const auto tasks = getTasks_();

    OSWFile osw(input_file);
    std::optional<std::vector<OpenSwathExportRow>> matrix_rows_cache;
    for (const auto& task : tasks)
    {
      ProgressLogger progress_logger;
      progress_logger.setLogType(ProgressLogger::CMD);
      progress_logger.startProgress(0, 1, taskLabel_(task));

      switch (task.type)
      {
        case ExportTaskType::FeatureParquet:
        {
          const auto parquet_config = getParquetConfig_();
          const auto feature_table = osw.readOpenSwathFeatureScoreTable(parquet_config);
          const String feature_out = base_path + ".precursor.feature.scores.parquet";
          OpenSwathParquetExporter::writeFeatureScores(feature_out, feature_table);
          OPENMS_LOG_INFO << "Wrote " << feature_table.rows.size() << " precursor feature score rows to '" << feature_out << "'." << std::endl;
          if (parquet_config.include_transition_data)
          {
            const auto transition_table = osw.readOpenSwathTransitionScoreTable(parquet_config);
            if (!transition_table.rows.empty())
            {
              const String transition_out = base_path + ".transition.feature.scores.parquet";
              OpenSwathParquetExporter::writeTransitionScores(transition_out, transition_table);
              OPENMS_LOG_INFO << "Wrote " << transition_table.rows.size() << " transition feature score rows to '" << transition_out << "'." << std::endl;
            }
          }
          break;
        }
        case ExportTaskType::Results:
        {
          const auto results_config = getResultsConfig_();
          const auto rows = osw.readOpenSwathExportRows(results_config.filters);
          const String suffix = results_config.format == OpenSwathExportFileFormat::Parquet ? ".results.parquet" : ".results.tsv";
          const String output = base_path + suffix;
          OpenSwathResultsExporter::write(output, rows, results_config);
          OPENMS_LOG_INFO << "Wrote " << rows.size() << " filtered result rows to '" << output << "'." << std::endl;
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
          const String suffix = matrix_config.format == OpenSwathExportFileFormat::Parquet ? ".matrix.parquet" : ".matrix.tsv";
          const String output = base_path + "." + toString(*task.matrix_level) + suffix;
          OpenSwathMatrixExporter::writeMatrix(output, matrix, matrix_config);
          OPENMS_LOG_INFO << "Wrote " << matrix.identifier_rows.size() << " " << toString(*task.matrix_level)
                          << " matrix rows to '" << output << "'." << std::endl;
          break;
        }
      }
      progress_logger.endProgress();
    }

    return EXECUTION_OK;
  }
};
/// @endcond

int main(int argc, const char** argv)
{
  TOPPOpenSwathExport tool;
  return tool.main(argc, argv);
}
