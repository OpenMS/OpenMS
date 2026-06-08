// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathGeneInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathProteinInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptideInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <map>
#include <optional>
#include <utility>
#include <vector>

using namespace OpenMS;
using namespace std;

/**
@page TOPP_OpenSwathInfer OpenSwathInfer

@brief Runs peptidoform, peptide, protein, or gene inference on OpenSWATH OSW files.

The tool operates directly on OSW files. Peptidoform inference writes @c SCORE_IPF.
Peptide, protein, and gene inference write @c SCORE_PEPTIDE, @c SCORE_PROTEIN, and
@c SCORE_GENE for the selected context. Multiple inference passes can be selected
in one invocation via per-level @c true|false selectors. By default,
@c peptide:global and @c protein:global are enabled.

For more information, see the original publications describing the underlying algorithms:
- Rosenberger, G., Bludau, I., Schmitt, U. et al. Statistical control of peptide and protein error rates in large-scale targeted data-independent acquisition analyses. Nat Methods 14, 921–927 (2017). https://doi.org/10.1038/nmeth.4398
- Rosenberger, G., Liu, Y., Röst, H. et al. Inference and quantification of peptidoforms in large sample cohorts by SWATH-MS. Nat Biotechnol 35, 781–788 (2017). https://doi.org/10.1038/nbt.3908

<B>The command line parameters of this tool are:</B>
@verbinclude TOPP_OpenSwathInfer.cli
<B>INI file documentation of this tool:</B>
@htmlinclude TOPP_OpenSwathInfer.html
*/

/// @cond TOPPCLASSES
class TOPPOpenSwathInfer :
  public TOPPBase
{
public:
  TOPPOpenSwathInfer() :
    TOPPBase("OpenSwathInfer", "Runs inference on OpenSWATH results.")
  {
  }

protected:
  struct InferenceTask
  {
    InferenceLevel level = InferenceLevel::Peptidoform;
    std::optional<InferenceContext> context;
  };

  static constexpr double summary_fdr_threshold_ = 0.01;

  static bool toBool_(const String& value)
  {
    return value == "true";
  }

  PeptidoformInferenceConfig getPeptidoformConfig_() const
  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = toBool_(getStringOption_("peptidoform:ipf_ms1_scoring"));
    config.ipf_ms2_scoring = toBool_(getStringOption_("peptidoform:ipf_ms2_scoring"));
    config.ipf_h0 = toBool_(getStringOption_("peptidoform:ipf_h0"));
    config.ipf_grouped_fdr = toBool_(getStringOption_("peptidoform:ipf_grouped_fdr"));
    config.ipf_max_peakgroup_pep = getDoubleOption_("peptidoform:ipf_max_peakgroup_pep");
    config.ipf_max_precursor_pep = getDoubleOption_("peptidoform:ipf_max_precursor_pep");
    config.ipf_max_precursor_peakgroup_pep = getDoubleOption_("peptidoform:ipf_max_precursor_peakgroup_pep");
    config.ipf_max_transition_pep = getDoubleOption_("peptidoform:ipf_max_transition_pep");
    config.propagate_signal_across_runs = toBool_(getStringOption_("peptidoform:propagate_signal_across_runs"));
    config.ipf_min_alignment_mapping_confidence = getDoubleOption_("peptidoform:ipf_min_alignment_mapping_confidence");
    config.across_run_confidence_threshold = getDoubleOption_("peptidoform:across_run_confidence_threshold");
    return config;
  }

  ErrorEstimationConfig getErrorConfig_() const
  {
    ErrorEstimationConfig config;
    config.parametric = toBool_(getStringOption_("error_estimation:parametric"));
    config.pfdr = toBool_(getStringOption_("error_estimation:pfdr"));
    config.pi0_lambda = getDoubleList_("error_estimation:pi0_lambda");
    config.pi0_method = getStringOption_("error_estimation:pi0_method");
    config.pi0_smooth_df = getIntOption_("error_estimation:pi0_smooth_df");
    config.pi0_smooth_log_pi0 = toBool_(getStringOption_("error_estimation:pi0_smooth_log_pi0"));
    config.lfdr_truncate = toBool_(getStringOption_("error_estimation:lfdr_truncate"));
    config.lfdr_monotone = toBool_(getStringOption_("error_estimation:lfdr_monotone"));
    config.lfdr_transformation = getStringOption_("error_estimation:lfdr_transformation");
    config.lfdr_adj = getDoubleOption_("error_estimation:lfdr_adj");
    config.lfdr_eps = getDoubleOption_("error_estimation:lfdr_eps");
    return config;
  }

  static void appendTask_(std::vector<InferenceTask>& tasks, const InferenceLevel level, const std::optional<InferenceContext>& context)
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

  static String taskLabel_(const InferenceTask& task)
  {
    if (task.level == InferenceLevel::Peptidoform)
    {
      return "peptidoform inference";
    }
    return toString(task.level) + " inference in '" + toString(*task.context) + "' context";
  }

  static String targetLabelPlural_(const InferenceLevel level)
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

  static void logPeptidoformSummary_(const std::vector<IPFResultRow>& results)
  {
    const Size targets_at_1pct = static_cast<Size>(std::count_if(
      results.begin(),
      results.end(),
      [](const auto& row)
      {
        return std::isfinite(row.qvalue) && row.qvalue <= summary_fdr_threshold_;
      }));
    OPENMS_LOG_INFO << targetLabelPlural_(InferenceLevel::Peptidoform) << " at 1% FDR: " << targets_at_1pct << std::endl;
  }

  static void logLevelContextSummary_(const std::vector<LevelContextInputRow>& input_rows,
                                      const std::vector<LevelContextResultRow>& results,
                                      const InferenceLevel level,
                                      const InferenceContext context,
                                      const std::map<Int64, String>& run_basenames)
  {
    std::map<std::pair<Int64, Int64>, bool> is_target_by_key;
    for (const auto& row : input_rows)
    {
      is_target_by_key[{runKey_(row.run_id), row.entity_id}] = !row.decoy;
    }

    Size total_targets_at_1pct = 0;
    std::map<Int64, Size> per_run_targets_at_1pct;
    for (const auto& row : results)
    {
      if (!std::isfinite(row.qvalue) || row.qvalue > summary_fdr_threshold_)
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
      OPENMS_LOG_INFO <<  targetLabelPlural_(level) << " at 1% FDR across all runs: " << total_targets_at_1pct << std::endl;
    }

    for (const auto& run_count : per_run_targets_at_1pct)
    {
      const auto name_it = run_basenames.find(run_count.first);
      const String run_label = (name_it != run_basenames.end()) ? name_it->second : String("RUN_ID ") + String(run_count.first);
      OPENMS_LOG_INFO << run_label << ": " << run_count.second << " " << targetLabelPlural_(level) << std::endl;
    }
  }

  std::vector<InferenceTask> getInferenceTasks_() const
  {
    std::vector<InferenceTask> tasks;
    if (toBool_(getStringOption_("peptidoform:run")))
    {
      appendTask_(tasks, InferenceLevel::Peptidoform, std::nullopt);
    }

    const std::array<std::pair<String, InferenceLevel>, 3> level_prefixes =
    {{
      {"peptide", InferenceLevel::Peptide},
      {"protein", InferenceLevel::Protein},
      {"gene", InferenceLevel::Gene}
    }};
    const std::array<std::pair<String, InferenceContext>, 3> contexts =
    {{
      {"global", InferenceContext::Global},
      {"experiment_wide", InferenceContext::ExperimentWide},
      {"run_specific", InferenceContext::RunSpecific}
    }};

    for (const auto& level_prefix : level_prefixes)
    {
      for (const auto& context : contexts)
      {
        if (toBool_(getStringOption_(level_prefix.first + ":" + context.first)))
        {
          appendTask_(tasks, level_prefix.second, context.second);
        }
      }
    }

    if (tasks.empty())
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "No inference levels were enabled. Set at least one of the peptidoform/peptide/protein/gene selectors to 'true'.");
    }
    return tasks;
  }

  static String prepareWorkingOutput_(const String& input_file, const String& output_file)
  {
    if (output_file.empty() || output_file == input_file)
    {
      return input_file;
    }

    if (File::exists(output_file) && !File::remove(output_file))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, output_file);
    }
    if (!File::copy(input_file, output_file))
    {
      throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, output_file);
    }
    return output_file;
  }

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input OpenSWATH OSW file.");
    setValidFormats_("in", {"osw"});
    registerOutputFile_("out", "<file>", "", "Optional output OSW file. If empty, the input file is updated in place.", false);
    setValidFormats_("out", {"osw"});

    registerTOPPSubsection_("peptidoform", "Peptidoform/IPF inference options.");
    registerStringOption_("peptidoform:run", "<true|false>", "false", "Enable peptidoform inference in this invocation.", false);
    setValidStrings_("peptidoform:run", {"true", "false"});
    const PeptidoformInferenceConfig ipf_defaults;
    registerStringOption_("peptidoform:ipf_ms1_scoring", "<true|false>", ipf_defaults.ipf_ms1_scoring ? "true" : "false", "Use MS1 precursor evidence.", false);
    setValidStrings_("peptidoform:ipf_ms1_scoring", {"true", "false"});
    registerStringOption_("peptidoform:ipf_ms2_scoring", "<true|false>", ipf_defaults.ipf_ms2_scoring ? "true" : "false", "Use MS2 precursor evidence from transition-level precursor rows.", false);
    setValidStrings_("peptidoform:ipf_ms2_scoring", {"true", "false"});
    registerStringOption_("peptidoform:ipf_h0", "<true|false>", ipf_defaults.ipf_h0 ? "true" : "false", "Include the null hypothesis in IPF transition inference.", false, true);
    setValidStrings_("peptidoform:ipf_h0", {"true", "false"});
    registerStringOption_("peptidoform:ipf_grouped_fdr", "<true|false>", ipf_defaults.ipf_grouped_fdr ? "true" : "false", "Compute model FDR separately per num_peptidoforms group.", false, true);
    setValidStrings_("peptidoform:ipf_grouped_fdr", {"true", "false"});
    registerDoubleOption_("peptidoform:ipf_max_peakgroup_pep", "<num>", ipf_defaults.ipf_max_peakgroup_pep, "Maximum peakgroup PEP for entering precursor inference.", false, true);
    registerDoubleOption_("peptidoform:ipf_max_precursor_pep", "<num>", ipf_defaults.ipf_max_precursor_pep, "Maximum precursor-level evidence PEP for precursor inference.", false, true);
    registerDoubleOption_("peptidoform:ipf_max_precursor_peakgroup_pep", "<num>", ipf_defaults.ipf_max_precursor_peakgroup_pep, "Maximum precursor-peakgroup PEP retained for transition inference.", false, true);
    registerDoubleOption_("peptidoform:ipf_max_transition_pep", "<num>", ipf_defaults.ipf_max_transition_pep, "Maximum transition evidence PEP retained for IPF.", false, true);
    registerStringOption_("peptidoform:propagate_signal_across_runs", "<true|false>", ipf_defaults.propagate_signal_across_runs ? "true" : "false", "Propagate confident transition evidence across aligned runs using FEATURE_MS2_ALIGNMENT_CANDIDATE.", false);
    setValidStrings_("peptidoform:propagate_signal_across_runs", {"true", "false"});
    registerDoubleOption_("peptidoform:ipf_min_alignment_mapping_confidence", "<num>", ipf_defaults.ipf_min_alignment_mapping_confidence, "Minimum mapping confidence retained for across-run propagation from FEATURE_MS2_ALIGNMENT_CANDIDATE.", false, true);
    registerDoubleOption_("peptidoform:across_run_confidence_threshold", "<num>", ipf_defaults.across_run_confidence_threshold, "Maximum transition PEP eligible for across-run propagation.", false, true);

    registerTOPPSubsection_("peptide", "Peptide-level inference selection options.");
    registerStringOption_("peptide:global", "<true|false>", "true", "Enable peptide inference in the global context.", false);
    setValidStrings_("peptide:global", {"true", "false"});
    registerStringOption_("peptide:experiment_wide", "<true|false>", "false", "Enable peptide inference in the experiment-wide context.", false);
    setValidStrings_("peptide:experiment_wide", {"true", "false"});
    registerStringOption_("peptide:run_specific", "<true|false>", "false", "Enable peptide inference in the run-specific context.", false);
    setValidStrings_("peptide:run_specific", {"true", "false"});

    registerTOPPSubsection_("protein", "Protein-level inference selection options.");
    registerStringOption_("protein:global", "<true|false>", "true", "Enable protein inference in the global context.", false);
    setValidStrings_("protein:global", {"true", "false"});
    registerStringOption_("protein:experiment_wide", "<true|false>", "false", "Enable protein inference in the experiment-wide context.", false);
    setValidStrings_("protein:experiment_wide", {"true", "false"});
    registerStringOption_("protein:run_specific", "<true|false>", "false", "Enable protein inference in the run-specific context.", false);
    setValidStrings_("protein:run_specific", {"true", "false"});

    registerTOPPSubsection_("gene", "Gene-level inference selection options.");
    registerStringOption_("gene:global", "<true|false>", "false", "Enable gene inference in the global context.", false);
    setValidStrings_("gene:global", {"true", "false"});
    registerStringOption_("gene:experiment_wide", "<true|false>", "false", "Enable gene inference in the experiment-wide context.", false);
    setValidStrings_("gene:experiment_wide", {"true", "false"});
    registerStringOption_("gene:run_specific", "<true|false>", "false", "Enable gene inference in the run-specific context.", false);
    setValidStrings_("gene:run_specific", {"true", "false"});

    registerTOPPSubsection_("error_estimation", "Context-level error estimation options.");
    const ErrorEstimationConfig error_defaults;
    registerStringOption_("error_estimation:parametric", "<true|false>", error_defaults.parametric ? "true" : "false", "Use parametric p-value estimation instead of empirical estimation.", false, true);
    setValidStrings_("error_estimation:parametric", {"true", "false"});
    registerStringOption_("error_estimation:pfdr", "<true|false>", error_defaults.pfdr ? "true" : "false", "Use positive FDR in q-value and s-value calculations.", false, true);
    setValidStrings_("error_estimation:pfdr", {"true", "false"});
    registerDoubleList_("error_estimation:pi0_lambda", "<l1,l2,...>", error_defaults.pi0_lambda, "Lambda grid for pi0 estimation.", false, true);
    registerStringOption_("error_estimation:pi0_method", "<choice>", error_defaults.pi0_method, "Method for pi0 estimation.", false, true);
    setValidStrings_("error_estimation:pi0_method", {"bootstrap", "smoother"});
    registerIntOption_("error_estimation:pi0_smooth_df", "<num>", error_defaults.pi0_smooth_df, "Degrees of freedom for smoother-based pi0 estimation.", false, true);
    setMinInt_("error_estimation:pi0_smooth_df", 1);
    registerStringOption_("error_estimation:pi0_smooth_log_pi0", "<true|false>", error_defaults.pi0_smooth_log_pi0 ? "true" : "false", "Smooth log(pi0) instead of pi0 directly.", false, true);
    setValidStrings_("error_estimation:pi0_smooth_log_pi0", {"true", "false"});
    registerStringOption_("error_estimation:lfdr_truncate", "<true|false>", error_defaults.lfdr_truncate ? "true" : "false", "Truncate local FDR/PEP values to [0,1].", false, true);
    setValidStrings_("error_estimation:lfdr_truncate", {"true", "false"});
    registerStringOption_("error_estimation:lfdr_monotone", "<true|false>", error_defaults.lfdr_monotone ? "true" : "false", "Enforce monotone local FDR/PEP estimates.", false, true);
    setValidStrings_("error_estimation:lfdr_monotone", {"true", "false"});
    registerStringOption_("error_estimation:lfdr_transformation", "<choice>", error_defaults.lfdr_transformation, "Transformation used by local FDR estimation.", false, true);
    setValidStrings_("error_estimation:lfdr_transformation", {"probit", "logit"});
    registerDoubleOption_("error_estimation:lfdr_adj", "<num>", error_defaults.lfdr_adj, "Adjustment factor for local FDR estimation bandwidth.", false, true);
    registerDoubleOption_("error_estimation:lfdr_eps", "<num>", error_defaults.lfdr_eps, "Epsilon used to clip p-values before local FDR estimation.", false, true);
  }

  ExitCodes main_(int, const char**) override
  {
    const String input_file = getStringOption_("in");
    const String output_file = getStringOption_("out");
    const std::vector<InferenceTask> tasks = getInferenceTasks_();
    const String working_file = prepareWorkingOutput_(input_file, output_file);

    OSWFile osw(working_file);
    const PeptidoformInferenceConfig peptidoform_config = getPeptidoformConfig_();
    const ErrorEstimationConfig error_config = getErrorConfig_();
    const bool has_run_specific_task = std::any_of(tasks.begin(), tasks.end(),
      [](const auto& task)
      {
        return task.context.has_value() && *task.context == InferenceContext::RunSpecific;
      });
    const std::map<Int64, String> run_basenames = has_run_specific_task ? osw.readRunBasenames() : std::map<Int64, String>{};
    for (const auto& task : tasks)
    {
      const String task_label = taskLabel_(task);
      ProgressLogger progresslogger;
      progresslogger.setLogType(ProgressLogger::CMD);
      progresslogger.startProgress(0, 1, task_label);

      if (task.level == InferenceLevel::Peptidoform)
      {
        const auto precursor_rows = osw.readIPFPrecursorData(peptidoform_config);
        const auto transition_rows = osw.readIPFTransitionData(peptidoform_config);
        const auto alignment_rows = osw.readIPFAlignmentData(peptidoform_config);

        OpenSwathPeptidoformInference inference;
        const auto results = inference.infer(precursor_rows, transition_rows, alignment_rows, peptidoform_config);
        osw.writeIPFResults("", results);
        logPeptidoformSummary_(results);
        progresslogger.endProgress();
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
      progresslogger.endProgress();
    }

    return EXECUTION_OK;
  }
};
/// @endcond

int main(int argc, const char** argv)
{
  TOPPOpenSwathInfer tool;
  return tool.main(argc, argv);
}
