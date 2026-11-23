// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing, Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/OSWParquetFile.h>
#include <OpenMS/FORMAT/ParquetFile.h>
#include <OpenMS/FORMAT/ZipArchiveFile.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <iterator>
#include <limits>
#include <optional>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace OpenMS
{
namespace
{
  struct OpenSwathCanonicalLibraryMapping
  {
    std::unordered_map<std::string, int64_t> compound_to_precursor;
    std::unordered_map<int64_t, double> precursor_mz_by_id;
    std::unordered_map<int64_t, bool> precursor_decoy_by_id;
    std::unordered_map<std::string, int64_t> transition_to_id;
  };

  inline OpenSwathCanonicalLibraryMapping buildOpenSwathCanonicalLibraryMapping(const OpenSwath::LightTargetedExperiment& targeted_exp)
  {
    OpenSwathCanonicalLibraryMapping mapping;

    // IDs are already canonical dense integer strings. Preserve that exact ID
    // domain because Parquet writers persist compound.id directly.
    mapping.compound_to_precursor.reserve(targeted_exp.compounds.size());
    for (const auto& compound : targeted_exp.compounds)
    {
      mapping.compound_to_precursor.emplace(compound.id, StringUtils::toInt64(compound.id));
    }

    mapping.precursor_mz_by_id.reserve(targeted_exp.compounds.size());
    mapping.precursor_decoy_by_id.reserve(targeted_exp.compounds.size());
    mapping.transition_to_id.reserve(targeted_exp.transitions.size());
    for (const auto& transition : targeted_exp.transitions)
    {
      const auto precursor_it = mapping.compound_to_precursor.find(transition.peptide_ref);
      if (precursor_it == mapping.compound_to_precursor.end())
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Transition references unknown peptide_ref '" + std::string(transition.peptide_ref) + "'");
      }

      const int64_t precursor_id = precursor_it->second;
      if (! mapping.precursor_mz_by_id.contains(precursor_id)) { mapping.precursor_mz_by_id.emplace(precursor_id, transition.precursor_mz); }
      if (! mapping.precursor_decoy_by_id.contains(precursor_id) && transition.isDetectingTransition())
      {
        mapping.precursor_decoy_by_id.emplace(precursor_id, transition.getDecoy());
      }

      // transition_name is canonicalized to the integer transition ID.
      mapping.transition_to_id.try_emplace(transition.transition_name, StringUtils::toInt64(transition.transition_name));
    }

    return mapping;
  }
} // namespace

class OSWParquetFile::Impl
{
public:
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

  struct OSWPQWorkspace
  {
    std::string output_path;
    std::string base_dir;
    bool archive_input = false;
    bool dirty = false;
    std::unique_ptr<TempDir> temp_dir;
    mutable std::shared_ptr<arrow::Table> runs_table_cache;
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

  static std::pair<bool, std::vector<InferenceScoreColumn>> getInferenceScoreColumns_(const std::vector<InferenceTask>& tasks)
  {
    bool include_ipf_peptide_id = false;
    std::vector<InferenceScoreColumn> columns;

    const auto append_column = [&](const char* name, OptionalDoubleMember member) {
      const auto duplicate = std::find_if(columns.begin(), columns.end(), [&](const auto& column) { return std::string_view(column.name) == name; });
      if (duplicate == columns.end()) { columns.push_back({name, member}); }
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

  static bool parquetValuePresent_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  { return array != nullptr && ! array->IsNull(row); }

  static std::optional<double> parquetOptionalDouble_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    if (! parquetValuePresent_(array, row)) { return std::nullopt; }
    return ParquetFile::getDouble(array, row, 0.0, false);
  }

  static std::optional<Int64> parquetOptionalInt64_(const std::shared_ptr<arrow::Array>& array, const int64_t row)
  {
    if (! parquetValuePresent_(array, row)) { return std::nullopt; }
    return ParquetFile::getInt64(array, row, 0, false);
  }

  static std::shared_ptr<arrow::Array> getOptionalParquetColumn_(const std::unordered_map<std::string, std::shared_ptr<arrow::Array>>& columns,
                                                                 const std::string& name)
  {
    const auto it = columns.find(name);
    return it != columns.end() ? it->second : nullptr;
  }

  static std::shared_ptr<arrow::Table> getOSWPQRunsTable_(const OSWPQWorkspace& workspace)
  {
    if (workspace.runs_table_cache == nullptr) { workspace.runs_table_cache = ParquetFile::readTable(workspace.base_dir + "/runs/runs.parquet"); }
    return workspace.runs_table_cache;
  }

  // Deliberately not cached: callers process one run at a time, and keeping every
  // run's features table resident would make peak memory scale with the whole
  // experiment instead of the largest run.
  static std::shared_ptr<arrow::Table> getOSWPQFeatureTable_(const OSWPQWorkspace& workspace, const Int64 run_id)
  { return ParquetFile::readTable(workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet"); }

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
      if (columns_to_replace.contains(field->name())) { continue; }
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
      case InferenceLevel::Peptide:
        return "peptide_id";
      case InferenceLevel::Protein:
        return "protein_id";
      case InferenceLevel::Gene:
        return "gene_id";
      case InferenceLevel::Peptidoform:
        break;
    }
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Direct OSWPQ level-context helpers do not support peptidoform inference.");
  }

  static std::string inferenceParquetPath_(const OSWPQWorkspace& workspace, const InferenceLevel level)
  { return workspace.base_dir + "/inference/score_" + toString(level) + ".parquet"; }

  static bool betterLevelContextResult_(const LevelContextResultRow& candidate, const LevelContextResultRow& current)
  {
    if (candidate.pep != current.pep) { return candidate.pep < current.pep; }
    if (candidate.qvalue != current.qvalue) { return candidate.qvalue < current.qvalue; }
    if (candidate.score != current.score) { return candidate.score > current.score; }
    return candidate.entity_id < current.entity_id;
  }

  static void appendUniqueEntity_(std::vector<Int64>& entities, const Int64 entity_id)
  {
    if (std::find(entities.begin(), entities.end(), entity_id) == entities.end()) { entities.push_back(entity_id); }
  }

  static std::string joinStrings_(const std::vector<std::string>& values)
  {
    std::string joined;
    for (Size i = 0; i < values.size(); ++i)
    {
      if (i != 0) { joined += ";"; }
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
        if (protein_it != lookup.proteins.end()) { names.push_back(protein_it->second.accession); }
      }
      lookup.protein_names_by_peptide[peptide_id] = joinStrings_(names);
      if (protein_ids.size() == 1) { lookup.unique_protein_by_peptide[peptide_id] = protein_ids.front(); }
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
        if (gene_it != lookup.genes.end()) { names.push_back(gene_it->second.name); }
      }
      lookup.gene_names_by_peptide[peptide_id] = joinStrings_(names);
      if (gene_ids.size() == 1) { lookup.unique_gene_by_peptide[peptide_id] = gene_ids.front(); }
    }

    for (auto& [transition_id, transition] : lookup.transitions)
    {
      std::sort(transition.peptide_ids.begin(), transition.peptide_ids.end());
      transition.peptide_ids.erase(std::unique(transition.peptide_ids.begin(), transition.peptide_ids.end()), transition.peptide_ids.end());
    }
  }

  PreparedLibraryLookup_ buildPreparedLibraryLookupFromLightTargetedExperiment_(const OpenSwath::LightTargetedExperiment& targeted_exp,
                                                                                const bool load_transition_metadata) const
  {
    PreparedLibraryLookup_ lookup;
    const auto canonical_mapping = buildOpenSwathCanonicalLibraryMapping(targeted_exp);

    std::vector<std::string> peptide_sequences;
    peptide_sequences.reserve(targeted_exp.compounds.size() + targeted_exp.transitions.size());
    for (const auto& compound : targeted_exp.compounds)
    {
      if (compound.isPeptide()) { peptide_sequences.push_back(compound.sequence); }
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
      lookup.peptides.emplace(peptide_id, PreparedLibraryPeptide_ {peptide_id, unmodified_sequence, modified_sequence, false});
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
      lookup.proteins.emplace(protein_id, PreparedLibraryProtein_ {protein_id, accession, false});
    }

    std::unordered_map<std::string, Int64> gene_ids_by_name;
    gene_ids_by_name.reserve(targeted_exp.compounds.size());

    for (const auto& compound : targeted_exp.compounds)
    {
      const auto precursor_id_it = canonical_mapping.compound_to_precursor.find(compound.id);
      if (precursor_id_it == canonical_mapping.compound_to_precursor.end()) { continue; }

      const Int64 precursor_id = precursor_id_it->second;
      PreparedLibraryPrecursor_ precursor;
      precursor.precursor_id = precursor_id;
      precursor.traml_id = compound.id;
      precursor.group_label = compound.isPeptide() ? compound.peptide_group_label : "";
      const auto precursor_mz_it = canonical_mapping.precursor_mz_by_id.find(precursor_id);
      precursor.precursor_mz = precursor_mz_it != canonical_mapping.precursor_mz_by_id.end() ? precursor_mz_it->second : 0.0;
      precursor.charge = static_cast<Int32>(compound.charge);
      if (std::isfinite(compound.rt)) { precursor.library_rt = compound.rt; }
      if (compound.drift_time != -1 && std::isfinite(compound.drift_time)) { precursor.library_drift_time = compound.drift_time; }
      const auto precursor_decoy_it = canonical_mapping.precursor_decoy_by_id.find(precursor_id);
      precursor.decoy = precursor_decoy_it != canonical_mapping.precursor_decoy_by_id.end() ? precursor_decoy_it->second : false;
      lookup.precursors[precursor_id] = std::move(precursor);

      if (! compound.isPeptide()) { continue; }

      const auto peptide_id_it = peptide_ids_by_sequence.find(compound.sequence);
      if (peptide_id_it == peptide_ids_by_sequence.end()) { continue; }

      const Int64 peptide_id = peptide_id_it->second;
      appendUniqueEntity_(lookup.precursor_to_peptides[precursor_id], peptide_id);
      lookup.peptides[peptide_id].decoy = lookup.peptides[peptide_id].decoy || lookup.precursors.at(precursor_id).decoy;

      auto& protein_ids = lookup.peptide_to_proteins[peptide_id];
      for (const auto& protein_ref : compound.protein_refs)
      {
        const auto protein_id_it = protein_ids_by_accession.find(protein_ref);
        if (protein_id_it != protein_ids_by_accession.end()) { appendUniqueEntity_(protein_ids, protein_id_it->second); }
      }

      const std::string gene_name = compound.gene_name.empty() ? "NA" : compound.gene_name;
      auto [gene_it, inserted] = gene_ids_by_name.try_emplace(gene_name, static_cast<Int64>(gene_ids_by_name.size()));
      if (inserted) { lookup.genes.emplace(gene_it->second, PreparedLibraryGene_ {gene_it->second, gene_name, false}); }
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
      for (const auto& transition : targeted_exp.transitions)
      {
        const auto precursor_id_it = canonical_mapping.compound_to_precursor.find(transition.peptide_ref);
        if (precursor_id_it == canonical_mapping.compound_to_precursor.end()) { continue; }

        PreparedLibraryTransition_ transition_entry;
        transition_entry.transition_id = StringUtils::toInt64(transition.transition_name);
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
          if (peptide_id_it != peptide_ids_by_sequence.end()) { appendUniqueEntity_(transition_entry.peptide_ids, peptide_id_it->second); }
        }

        lookup.transitions.emplace(transition_entry.transition_id, std::move(transition_entry));
      }
    }

    finalizePreparedLibraryLookup_(lookup);
    return lookup;
  }

  static std::vector<Int64> mappedEntitiesForFeature_(const PreparedLibraryLookup_& lookup, const InferenceLevel level, const Int64 precursor_id)
  {
    std::vector<Int64> entities;
    const auto peptide_it = lookup.precursor_to_peptides.find(precursor_id);
    if (peptide_it == lookup.precursor_to_peptides.end()) { return entities; }

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
        if (protein_it != lookup.unique_protein_by_peptide.end()) { appendUniqueEntity_(entities, protein_it->second); }
        continue;
      }
      if (level == InferenceLevel::Gene)
      {
        const auto gene_it = lookup.unique_gene_by_peptide.find(peptide_id);
        if (gene_it != lookup.unique_gene_by_peptide.end()) { appendUniqueEntity_(entities, gene_it->second); }
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

      if (! candidate.has_value()) { continue; }
      if (! best.has_value() || betterLevelContextResult_(*candidate, *best)) { best = candidate; }
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
          if (row.run_id.has_value()) { maps.experiment_wide[{*row.run_id, row.entity_id}] = row; }
          break;
        case InferenceContext::RunSpecific:
          if (row.run_id.has_value()) { maps.run_specific[{*row.run_id, row.entity_id}] = row; }
          break;
      }
    }
    return maps;
  }

  static void
  writeLevelContextResultsParquet_(const OSWPQWorkspace& workspace, const InferenceLevel level, const std::vector<LevelContextResultRow>& results)
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
      if (row.run_id.has_value()) { ParquetFile::appendOrThrow(run_id_builder.Append(*row.run_id), "run_id"); }
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

    std::vector<std::shared_ptr<arrow::Field>> fields = {arrow::field("context", arrow::utf8(), false),
                                                         arrow::field("run_id", arrow::int64(), true),
                                                         arrow::field(entityIdColumnName_(level), arrow::int64(), false),
                                                         arrow::field("score", arrow::float64(), false),
                                                         arrow::field("pvalue", arrow::float64(), false),
                                                         arrow::field("qvalue", arrow::float64(), false),
                                                         arrow::field("pep", arrow::float64(), false)};
    std::vector<std::shared_ptr<arrow::Array>> arrays = {ParquetFile::finishArray(context_builder, "context"),
                                                         ParquetFile::finishArray(run_id_builder, "run_id"),
                                                         ParquetFile::finishArray(entity_id_builder, entityIdColumnName_(level)),
                                                         ParquetFile::finishArray(score_builder, "score"),
                                                         ParquetFile::finishArray(pvalue_builder, "pvalue"),
                                                         ParquetFile::finishArray(qvalue_builder, "qvalue"),
                                                         ParquetFile::finishArray(pep_builder, "pep")};

    ParquetFile::writeTable(arrow::Table::Make(arrow::schema(fields), arrays), inferenceParquetPath_(workspace, level));
  }

  static std::vector<LevelContextResultRow> readLevelContextResultsParquet_(const OSWPQWorkspace& workspace, const InferenceLevel level)
  {
    const std::string file_path = inferenceParquetPath_(workspace, level);
    if (! File::exists(file_path)) { return {}; }

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
      if (context == "global") { result.context = InferenceContext::Global; }
      else if (context == "experiment-wide") { result.context = InferenceContext::ExperimentWide; }
      else
      {
        result.context = InferenceContext::RunSpecific;
      }
      if (run_id_col != nullptr && ! run_id_col->IsNull(row)) { result.run_id = ParquetFile::getInt64(run_id_col, row, 0, false); }
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
    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    std::map<Int64, std::string> basenames;
    for (int64_t row = 0; row < runs_table->num_rows(); ++row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, row, 0, false);
      std::string filename = (filename_col != nullptr && ! filename_col->IsNull(row)) ? ParquetFile::getString(filename_col, row) : "";
      std::string basename = File::stemName(filename);
      if (basename.empty()) { basename = File::basename(filename); }
      if (basename.empty()) { basename = "RUN_ID " + StringUtils::toStr(run_id); }
      basenames[run_id] = basename;
    }
    return basenames;
  }

  static const std::array<const char*, 17>& featureMS1ParquetFields_()
  {
    static const std::array<const char*, 17> fields
      = {{"ms1_area_intensity", "ms1_apex_intensity", "ms1_exp_im", "ms1_delta_im", "var_ms1_massdev_score", "var_ms1_im_ms1_delta_score",
          "var_ms1_mi_score", "var_ms1_mi_contrast_score", "var_ms1_mi_combined_score", "var_ms1_isotope_correlation_score",
          "var_ms1_isotope_overlap_score", "var_ms1_xcorr_coelution", "var_ms1_xcorr_coelution_contrast", "var_ms1_xcorr_coelution_combined",
          "var_ms1_xcorr_shape", "var_ms1_xcorr_shape_contrast", "var_ms1_xcorr_shape_combined"}};
    return fields;
  }

  static const std::array<const char*, 37>& featureMS2ParquetFields_()
  {
    static const std::array<const char*, 37> fields = {{"ms2_area_intensity",
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
                                                        "var_ms2_im_log_intensity"}};
    return fields;
  }

  static const std::array<const char*, 41>& featureTransitionParquetFields_()
  {
    static const std::array<const char*, 41> fields = {{"area_intensity",
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
                                                        "points_across_half_height"}};
    return fields;
  }

  OSWPQWorkspace prepareOSWPQWorkspace_(const std::string& workflow_oswpq) const
  {
    OSWPQWorkspace workspace;
    workspace.output_path = workflow_oswpq;
    workspace.archive_input = ! File::isDirectory(workflow_oswpq);
    if (workspace.archive_input) { workspace.base_dir = ZipArchiveFile::unzipDirectory(workflow_oswpq, workspace.temp_dir); }
    else
    {
      workspace.base_dir = workflow_oswpq;
    }
    return workspace;
  }

  void commitOSWPQWorkspace_(OSWPQWorkspace& workspace) const
  {
    if (! workspace.dirty) { return; }
    if (workspace.archive_input)
    {
      OPENMS_LOG_INFO << "Repacking workflow.oswpq archive." << std::endl;
      if (File::isDirectory(workspace.output_path))
      {
        throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, workspace.output_path);
      }
      // Build the complete replacement beside the destination. Keep its basename
      // so the embedded sidecar index retains the archive's actual name.
      TempDir staging_dir(File::path(File::absolutePath(workspace.output_path)));
      const std::string staged = staging_dir.getPath() + "/" + File::basename(workspace.output_path);
      ZipArchiveFile::zipDirectory(workspace.base_dir, staged);
      ZipArchiveFile::writeSidecarIndex(staged);

      std::error_code error;
      std::filesystem::rename(to_path(staged), to_path(workspace.output_path), error);
      if (error)
      {
        // Windows may reject rename-over-existing. Keep a recoverable backup
        // outside staging_dir until the replacement has succeeded.
        const std::string backup = workspace.output_path + ".bak." + File::getUniqueName(false);
        if (! File::exists(workspace.output_path) || ! File::rename(workspace.output_path, backup, false))
        {
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, workspace.output_path);
        }
        if (! File::rename(staged, workspace.output_path, false))
        {
          if (! File::rename(backup, workspace.output_path, false))
          {
            OPENMS_LOG_ERROR << "Could not restore OSWPQ archive; original remains at: " << backup << '\n';
          }
          throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, workspace.output_path);
        }
        if (! File::remove(backup)) { OPENMS_LOG_WARN << "Could not remove OSWPQ archive backup: " << backup << '\n'; }
      }
      OPENMS_LOG_INFO << "Finished repacking workflow.oswpq archive." << std::endl;
    }
    workspace.dirty = false;
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

    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    std::map<std::pair<Int64, Int64>, LevelContextInputRow> best_rows;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = getOSWPQFeatureTable_(workspace, run_id);
      const auto precursor_id_array = ParquetFile::getColumn(features_table, "precursor_id");
      const auto score_ms2_array = ParquetFile::getOptionalColumn(features_table, "score_ms2_score");
      if (score_ms2_array == nullptr)
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Level-context inference on OSWPQ requires score_ms2_score in features.parquet.");
      }

      for (int64_t row = 0; row < features_table->num_rows(); ++row)
      {
        if (score_ms2_array->IsNull(row)) { continue; }

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
          if (context != InferenceContext::Global) { candidate.run_id = run_id; }
          candidate.group_id
            = context == InferenceContext::Global ? StringUtils::toStr(entity_id) : StringUtils::toStr(run_id) + "_" + StringUtils::toStr(entity_id);
          candidate.entity_id = entity_id;
          candidate.decoy = precursor_it->second.decoy;
          candidate.score = score;
          candidate.context = context;

          const auto existing = best_rows.find(key);
          if (existing == best_rows.end() || candidate.score > existing->second.score) { best_rows[key] = std::move(candidate); }
        }
      }
    }

    std::vector<LevelContextInputRow> rows;
    rows.reserve(best_rows.size());
    for (auto& [key, row] : best_rows)
    {
      rows.push_back(std::move(row));
    }

    OPENMS_LOG_INFO << "Read " << rows.size() << " best-score rows for " << toString(level) << " inference in '" << toString(context) << "' context."
                    << std::endl;
    return rows;
  }

  void applyLevelContextResultsToOSWPQ_(OSWPQWorkspace& workspace,
                                        const PreparedLibraryLookup_& lookup,
                                        const std::map<InferenceLevel, std::vector<LevelContextResultRow>>& results_by_level,
                                        const std::vector<InferenceTask>& tasks,
                                        ProgressLogger::LogType log_type) const
  {
    if (tasks.empty()) { return; }

    const auto [include_ipf_peptide_id, score_members] = getInferenceScoreColumns_(tasks);
    if (include_ipf_peptide_id)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ score sync currently supports peptide, protein, and gene inference only.");
    }
    if (score_members.empty()) { return; }

    std::map<InferenceLevel, LevelContextResultMaps_> result_maps;
    for (const auto& [level, results] : results_by_level)
    {
      result_maps[level] = buildLevelContextResultMaps_(results);
    }

    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_array = ParquetFile::getColumn(runs_table, "run_id");

    std::unordered_set<std::string> replace_columns;
    replace_columns.reserve(score_members.size());
    for (const auto& score_member : score_members)
    {
      replace_columns.insert(score_member.name);
    }

    ProgressLogger progress_logger;
    progress_logger.setLogType(log_type);
    progress_logger.startProgress(0, runs_table->num_rows(), "syncing level-context scores back to workflow.oswpq");

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_array, run_row, 0, false);
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = getOSWPQFeatureTable_(workspace, run_id);
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

        const auto assign_level = [&](const InferenceLevel level, const InferenceContext context, OptionalDoubleMember score_member,
                                      OptionalDoubleMember pvalue_member, OptionalDoubleMember qvalue_member, OptionalDoubleMember pep_member) {
          const auto level_it = result_maps.find(level);
          if (level_it == result_maps.end()) { return; }
          const auto entity_ids = mappedEntitiesForFeature_(lookup, level, precursor_id);
          const auto best_result = selectBestLevelContextResult_(entity_ids, run_id, level_it->second, context);
          if (! best_result.has_value()) { return; }
          score_row.*score_member = best_result->score;
          score_row.*pvalue_member = best_result->pvalue;
          score_row.*qvalue_member = best_result->qvalue;
          score_row.*pep_member = best_result->pep;
        };

        assign_level(InferenceLevel::Peptide, InferenceContext::Global, &OpenSwathFeatureScoreRow::score_peptide_global_score,
                     &OpenSwathFeatureScoreRow::score_peptide_global_pvalue, &OpenSwathFeatureScoreRow::score_peptide_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_global_pep);
        assign_level(InferenceLevel::Peptide, InferenceContext::ExperimentWide, &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pvalue, &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_experiment_wide_pep);
        assign_level(InferenceLevel::Peptide, InferenceContext::RunSpecific, &OpenSwathFeatureScoreRow::score_peptide_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_pvalue, &OpenSwathFeatureScoreRow::score_peptide_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_peptide_run_specific_pep);

        assign_level(InferenceLevel::Protein, InferenceContext::Global, &OpenSwathFeatureScoreRow::score_protein_global_score,
                     &OpenSwathFeatureScoreRow::score_protein_global_pvalue, &OpenSwathFeatureScoreRow::score_protein_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_global_pep);
        assign_level(InferenceLevel::Protein, InferenceContext::ExperimentWide, &OpenSwathFeatureScoreRow::score_protein_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pvalue, &OpenSwathFeatureScoreRow::score_protein_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_experiment_wide_pep);
        assign_level(InferenceLevel::Protein, InferenceContext::RunSpecific, &OpenSwathFeatureScoreRow::score_protein_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_pvalue, &OpenSwathFeatureScoreRow::score_protein_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_protein_run_specific_pep);

        assign_level(InferenceLevel::Gene, InferenceContext::Global, &OpenSwathFeatureScoreRow::score_gene_global_score,
                     &OpenSwathFeatureScoreRow::score_gene_global_pvalue, &OpenSwathFeatureScoreRow::score_gene_global_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_global_pep);
        assign_level(InferenceLevel::Gene, InferenceContext::ExperimentWide, &OpenSwathFeatureScoreRow::score_gene_experiment_wide_score,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pvalue, &OpenSwathFeatureScoreRow::score_gene_experiment_wide_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_experiment_wide_pep);
        assign_level(InferenceLevel::Gene, InferenceContext::RunSpecific, &OpenSwathFeatureScoreRow::score_gene_run_specific_score,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_pvalue, &OpenSwathFeatureScoreRow::score_gene_run_specific_qvalue,
                     &OpenSwathFeatureScoreRow::score_gene_run_specific_pep);

        for (Size column = 0; column < score_members.size(); ++column)
        {
          const auto value = score_row.*(score_members[column].member);
          if (value.has_value()) { ParquetFile::appendOrThrow(double_builders[column]->Append(*value), score_members[column].name); }
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
    auto update_min = [](auto& target, const auto& key, const double qvalue) {
      const auto existing = target.find(key);
      if (existing == target.end() || qvalue < existing->second) { target[key] = qvalue; }
    };

    for (const auto& row : results)
    {
      const auto peptide_it = entity_to_peptides.find(row.entity_id);
      if (peptide_it == entity_to_peptides.end()) { continue; }
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

  std::unordered_map<Int64, TransitionAggregation_>
  buildTransitionAggregations_(const OSWPQWorkspace& workspace, const PreparedLibraryLookup_& lookup, const double max_transition_pep) const
  {
    std::unordered_map<Int64, TransitionAggregation_> aggregations;
    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string feature_transition_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet";
      if (! File::exists(feature_transition_path)) { continue; }

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
          if (pep_col->IsNull(row)) { continue; }
          if (ParquetFile::getDouble(pep_col, row, 1.0, false) >= max_transition_pep) { continue; }
        }

        const Int64 transition_id = ParquetFile::getInt64(transition_id_col, row, 0, false);
        const auto transition_it = lookup.transitions.find(transition_id);
        if (transition_it == lookup.transitions.end()) { continue; }

        // Match OSWFile::readOpenSwathExportRows(): without transition-level
        // rescoring, aggregate every feature transition. Once transition PEP
        // scores are present, exclude decoy transitions and apply the PEP
        // threshold just like the SQLite SCORE_TRANSITION query.
        if (filter_by_transition_pep && transition_it->second.decoy) { continue; }

        const Int64 feature_id = ParquetFile::getInt64(feature_id_col, row, 0, false);
        auto& aggregation = aggregations[feature_id];
        aggregation.areas.push_back(StringUtils::toStr(ParquetFile::getDouble(area_col, row, 0.0, true)));
        aggregation.apices.push_back(StringUtils::toStr(ParquetFile::getDouble(apex_col, row, 0.0, true)));
        aggregation.annotations.push_back(StringUtils::toStr(transition_id) + "_" + transition_it->second.type
                                          + StringUtils::toStr(transition_it->second.ordinal) + "_"
                                          + StringUtils::toStr(transition_it->second.charge));
      }
    }

    return aggregations;
  }

  OpenSwathFeatureScoreTable
  readOSWPQFeatureScoreTable_(const OSWPQWorkspace& workspace, const PreparedLibraryLookup_& lookup, const OpenSwathParquetExportConfig& config) const
  {
    OpenSwathFeatureScoreTable table;
    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    bool discovered_dynamic_columns = false;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string filename = (filename_col != nullptr && ! filename_col->IsNull(run_row)) ? ParquetFile::getString(filename_col, run_row) : "";
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = getOSWPQFeatureTable_(workspace, run_id);
      const auto feature_id_col = ParquetFile::getColumn(features_table, "feature_id");
      const auto precursor_id_col = ParquetFile::getColumn(features_table, "precursor_id");

      std::unordered_map<std::string, std::shared_ptr<arrow::Array>> feature_columns;
      const auto add_feature_column
        = [&](const std::string& name) { feature_columns.emplace(name, ParquetFile::getOptionalColumn(features_table, name)); };
      for (const auto& name : std::array<std::string, 27> {"exp_rt",
                                                           "exp_im",
                                                           "norm_rt",
                                                           "delta_rt",
                                                           "left_width",
                                                           "right_width",
                                                           "exp_im_leftwidth",
                                                           "exp_im_rightwidth",
                                                           "score_ms1_score",
                                                           "score_ms1_peak_group_rank",
                                                           "score_ms1_pvalue",
                                                           "score_ms1_qvalue",
                                                           "score_ms1_pep",
                                                           "score_ms2_score",
                                                           "score_ms2_peak_group_rank",
                                                           "score_ms2_pvalue",
                                                           "score_ms2_qvalue",
                                                           "score_ms2_pep",
                                                           "ipf_peptide_id",
                                                           "score_ipf_precursor_peakgroup_pep",
                                                           "score_ipf_pep",
                                                           "score_ipf_qvalue",
                                                           "score_peptide_global_score",
                                                           "score_peptide_global_pvalue",
                                                           "score_peptide_global_qvalue",
                                                           "score_peptide_global_pep",
                                                           "score_peptide_experiment_wide_score"})
      {
        add_feature_column(name);
      }
      for (const auto& name : std::array<std::string, 20> {"score_peptide_experiment_wide_pvalue",
                                                           "score_peptide_experiment_wide_qvalue",
                                                           "score_peptide_experiment_wide_pep",
                                                           "score_peptide_run_specific_score",
                                                           "score_peptide_run_specific_pvalue",
                                                           "score_peptide_run_specific_qvalue",
                                                           "score_peptide_run_specific_pep",
                                                           "score_protein_global_score",
                                                           "score_protein_global_pvalue",
                                                           "score_protein_global_qvalue",
                                                           "score_protein_global_pep",
                                                           "score_protein_experiment_wide_score",
                                                           "score_protein_experiment_wide_pvalue",
                                                           "score_protein_experiment_wide_qvalue",
                                                           "score_protein_experiment_wide_pep",
                                                           "score_protein_run_specific_score",
                                                           "score_protein_run_specific_pvalue",
                                                           "score_protein_run_specific_qvalue",
                                                           "score_protein_run_specific_pep",
                                                           "score_gene_global_score"})
      {
        add_feature_column(name);
      }
      for (const auto& name : std::array<std::string, 11> {
             "score_gene_global_pvalue", "score_gene_global_qvalue", "score_gene_global_pep", "score_gene_experiment_wide_score",
             "score_gene_experiment_wide_pvalue", "score_gene_experiment_wide_qvalue", "score_gene_experiment_wide_pep",
             "score_gene_run_specific_score", "score_gene_run_specific_pvalue", "score_gene_run_specific_qvalue", "score_gene_run_specific_pep"})
      {
        add_feature_column(name);
      }

      if (! discovered_dynamic_columns)
      {
        for (const auto& name : featureMS1ParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(features_table, name) != nullptr) { table.feature_ms1_column_names.emplace_back(name); }
        }
        for (const auto& name : featureMS2ParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(features_table, name) != nullptr) { table.feature_ms2_column_names.emplace_back(name); }
        }
        discovered_dynamic_columns = true;
      }

      const auto delta_rt_col = getOptionalParquetColumn_(feature_columns, "delta_rt");
      const auto exp_im_col = getOptionalParquetColumn_(feature_columns, "exp_im");
      const auto exp_im_leftwidth_col = getOptionalParquetColumn_(feature_columns, "exp_im_leftwidth");
      const auto exp_im_rightwidth_col = getOptionalParquetColumn_(feature_columns, "exp_im_rightwidth");
      const auto exp_rt_col = getOptionalParquetColumn_(feature_columns, "exp_rt");
      const auto ipf_peptide_id_col = getOptionalParquetColumn_(feature_columns, "ipf_peptide_id");
      const auto left_width_col = getOptionalParquetColumn_(feature_columns, "left_width");
      const auto norm_rt_col = getOptionalParquetColumn_(feature_columns, "norm_rt");
      const auto right_width_col = getOptionalParquetColumn_(feature_columns, "right_width");
      const auto score_gene_experiment_wide_pep_col = getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_pep");
      const auto score_gene_experiment_wide_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_pvalue");
      const auto score_gene_experiment_wide_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_qvalue");
      const auto score_gene_experiment_wide_score_col = getOptionalParquetColumn_(feature_columns, "score_gene_experiment_wide_score");
      const auto score_gene_global_pep_col = getOptionalParquetColumn_(feature_columns, "score_gene_global_pep");
      const auto score_gene_global_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_global_pvalue");
      const auto score_gene_global_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_global_qvalue");
      const auto score_gene_global_score_col = getOptionalParquetColumn_(feature_columns, "score_gene_global_score");
      const auto score_gene_run_specific_pep_col = getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_pep");
      const auto score_gene_run_specific_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_pvalue");
      const auto score_gene_run_specific_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_qvalue");
      const auto score_gene_run_specific_score_col = getOptionalParquetColumn_(feature_columns, "score_gene_run_specific_score");
      const auto score_ipf_pep_col = getOptionalParquetColumn_(feature_columns, "score_ipf_pep");
      const auto score_ipf_precursor_peakgroup_pep_col = getOptionalParquetColumn_(feature_columns, "score_ipf_precursor_peakgroup_pep");
      const auto score_ipf_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_ipf_qvalue");
      const auto score_ms1_peak_group_rank_col = getOptionalParquetColumn_(feature_columns, "score_ms1_peak_group_rank");
      const auto score_ms1_pep_col = getOptionalParquetColumn_(feature_columns, "score_ms1_pep");
      const auto score_ms1_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_ms1_pvalue");
      const auto score_ms1_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_ms1_qvalue");
      const auto score_ms1_score_col = getOptionalParquetColumn_(feature_columns, "score_ms1_score");
      const auto score_ms2_peak_group_rank_col = getOptionalParquetColumn_(feature_columns, "score_ms2_peak_group_rank");
      const auto score_ms2_pep_col = getOptionalParquetColumn_(feature_columns, "score_ms2_pep");
      const auto score_ms2_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_ms2_pvalue");
      const auto score_ms2_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_ms2_qvalue");
      const auto score_ms2_score_col = getOptionalParquetColumn_(feature_columns, "score_ms2_score");
      const auto score_peptide_experiment_wide_pep_col = getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_pep");
      const auto score_peptide_experiment_wide_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_pvalue");
      const auto score_peptide_experiment_wide_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_qvalue");
      const auto score_peptide_experiment_wide_score_col = getOptionalParquetColumn_(feature_columns, "score_peptide_experiment_wide_score");
      const auto score_peptide_global_pep_col = getOptionalParquetColumn_(feature_columns, "score_peptide_global_pep");
      const auto score_peptide_global_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_global_pvalue");
      const auto score_peptide_global_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_global_qvalue");
      const auto score_peptide_global_score_col = getOptionalParquetColumn_(feature_columns, "score_peptide_global_score");
      const auto score_peptide_run_specific_pep_col = getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_pep");
      const auto score_peptide_run_specific_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_pvalue");
      const auto score_peptide_run_specific_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_qvalue");
      const auto score_peptide_run_specific_score_col = getOptionalParquetColumn_(feature_columns, "score_peptide_run_specific_score");
      const auto score_protein_experiment_wide_pep_col = getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_pep");
      const auto score_protein_experiment_wide_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_pvalue");
      const auto score_protein_experiment_wide_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_qvalue");
      const auto score_protein_experiment_wide_score_col = getOptionalParquetColumn_(feature_columns, "score_protein_experiment_wide_score");
      const auto score_protein_global_pep_col = getOptionalParquetColumn_(feature_columns, "score_protein_global_pep");
      const auto score_protein_global_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_global_pvalue");
      const auto score_protein_global_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_global_qvalue");
      const auto score_protein_global_score_col = getOptionalParquetColumn_(feature_columns, "score_protein_global_score");
      const auto score_protein_run_specific_pep_col = getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_pep");
      const auto score_protein_run_specific_pvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_pvalue");
      const auto score_protein_run_specific_qvalue_col = getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_qvalue");
      const auto score_protein_run_specific_score_col = getOptionalParquetColumn_(feature_columns, "score_protein_run_specific_score");
      std::vector<std::shared_ptr<arrow::Array>> feature_ms1_columns;
      feature_ms1_columns.reserve(table.feature_ms1_column_names.size());
      for (const auto& name : table.feature_ms1_column_names)
      {
        feature_ms1_columns.push_back(ParquetFile::getOptionalColumn(features_table, name));
      }
      std::vector<std::shared_ptr<arrow::Array>> feature_ms2_columns;
      feature_ms2_columns.reserve(table.feature_ms2_column_names.size());
      for (const auto& name : table.feature_ms2_column_names)
      {
        feature_ms2_columns.push_back(ParquetFile::getOptionalColumn(features_table, name));
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
        if (config.filters.exclude_decoys && precursor_it->second.decoy) { continue; }

        const auto peptide_mapping_it = lookup.precursor_to_peptides.find(precursor_id);
        if (peptide_mapping_it == lookup.precursor_to_peptides.end()) { continue; }

        for (const Int64 peptide_id : peptide_mapping_it->second)
        {
          const auto peptide_it = lookup.peptides.find(peptide_id);
          if (peptide_it == lookup.peptides.end()) { continue; }

          std::vector<std::optional<Int64>> protein_ids {std::nullopt};
          const auto protein_mapping_it = lookup.peptide_to_proteins.find(peptide_id);
          if (protein_mapping_it != lookup.peptide_to_proteins.end() && ! protein_mapping_it->second.empty())
          {
            protein_ids.clear();
            protein_ids.reserve(protein_mapping_it->second.size());
            for (const Int64 protein_id : protein_mapping_it->second)
            {
              protein_ids.push_back(protein_id);
            }
          }

          std::vector<std::optional<Int64>> gene_ids {std::nullopt};
          const auto gene_mapping_it = lookup.peptide_to_genes.find(peptide_id);
          if (gene_mapping_it != lookup.peptide_to_genes.end() && ! gene_mapping_it->second.empty())
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
              score_row.ipf_peptide_id = parquetOptionalInt64_(ipf_peptide_id_col, row);
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
              score_row.exp_rt = ParquetFile::getDouble(exp_rt_col, row, 0.0, true);
              score_row.exp_im = parquetOptionalDouble_(exp_im_col, row);
              score_row.norm_rt = ParquetFile::getDouble(norm_rt_col, row, 0.0, true);
              score_row.delta_rt = ParquetFile::getDouble(delta_rt_col, row, 0.0, true);
              score_row.left_width = ParquetFile::getDouble(left_width_col, row, 0.0, true);
              score_row.right_width = ParquetFile::getDouble(right_width_col, row, 0.0, true);
              score_row.im_left_width = parquetOptionalDouble_(exp_im_leftwidth_col, row);
              score_row.im_right_width = parquetOptionalDouble_(exp_im_rightwidth_col, row);

              for (const auto& column : feature_ms1_columns)
              {
                score_row.feature_ms1_values.push_back(ParquetFile::getDouble(column, row, 0.0, true));
              }
              for (const auto& column : feature_ms2_columns)
              {
                score_row.feature_ms2_values.push_back(ParquetFile::getDouble(column, row, 0.0, true));
              }

              score_row.score_ms1_score = parquetOptionalDouble_(score_ms1_score_col, row);
              const auto score_ms1_rank = parquetOptionalInt64_(score_ms1_peak_group_rank_col, row);
              if (score_ms1_rank.has_value()) score_row.score_ms1_rank = static_cast<Int32>(*score_ms1_rank);
              score_row.score_ms1_pvalue = parquetOptionalDouble_(score_ms1_pvalue_col, row);
              score_row.score_ms1_qvalue = parquetOptionalDouble_(score_ms1_qvalue_col, row);
              score_row.score_ms1_pep = parquetOptionalDouble_(score_ms1_pep_col, row);
              score_row.score_ms2_score = parquetOptionalDouble_(score_ms2_score_col, row);
              const auto score_ms2_rank = parquetOptionalInt64_(score_ms2_peak_group_rank_col, row);
              if (score_ms2_rank.has_value()) score_row.score_ms2_peak_group_rank = static_cast<Int32>(*score_ms2_rank);
              score_row.score_ms2_pvalue = parquetOptionalDouble_(score_ms2_pvalue_col, row);
              score_row.score_ms2_qvalue = parquetOptionalDouble_(score_ms2_qvalue_col, row);
              score_row.score_ms2_pep = parquetOptionalDouble_(score_ms2_pep_col, row);
              score_row.score_ipf_precursor_peakgroup_pep = parquetOptionalDouble_(score_ipf_precursor_peakgroup_pep_col, row);
              score_row.score_ipf_pep = parquetOptionalDouble_(score_ipf_pep_col, row);
              score_row.score_ipf_qvalue = parquetOptionalDouble_(score_ipf_qvalue_col, row);
              score_row.score_peptide_global_score = parquetOptionalDouble_(score_peptide_global_score_col, row);
              score_row.score_peptide_global_pvalue = parquetOptionalDouble_(score_peptide_global_pvalue_col, row);
              score_row.score_peptide_global_qvalue = parquetOptionalDouble_(score_peptide_global_qvalue_col, row);
              score_row.score_peptide_global_pep = parquetOptionalDouble_(score_peptide_global_pep_col, row);
              score_row.score_peptide_experiment_wide_score = parquetOptionalDouble_(score_peptide_experiment_wide_score_col, row);
              score_row.score_peptide_experiment_wide_pvalue = parquetOptionalDouble_(score_peptide_experiment_wide_pvalue_col, row);
              score_row.score_peptide_experiment_wide_qvalue = parquetOptionalDouble_(score_peptide_experiment_wide_qvalue_col, row);
              score_row.score_peptide_experiment_wide_pep = parquetOptionalDouble_(score_peptide_experiment_wide_pep_col, row);
              score_row.score_peptide_run_specific_score = parquetOptionalDouble_(score_peptide_run_specific_score_col, row);
              score_row.score_peptide_run_specific_pvalue = parquetOptionalDouble_(score_peptide_run_specific_pvalue_col, row);
              score_row.score_peptide_run_specific_qvalue = parquetOptionalDouble_(score_peptide_run_specific_qvalue_col, row);
              score_row.score_peptide_run_specific_pep = parquetOptionalDouble_(score_peptide_run_specific_pep_col, row);
              score_row.score_protein_global_score = parquetOptionalDouble_(score_protein_global_score_col, row);
              score_row.score_protein_global_pvalue = parquetOptionalDouble_(score_protein_global_pvalue_col, row);
              score_row.score_protein_global_qvalue = parquetOptionalDouble_(score_protein_global_qvalue_col, row);
              score_row.score_protein_global_pep = parquetOptionalDouble_(score_protein_global_pep_col, row);
              score_row.score_protein_experiment_wide_score = parquetOptionalDouble_(score_protein_experiment_wide_score_col, row);
              score_row.score_protein_experiment_wide_pvalue = parquetOptionalDouble_(score_protein_experiment_wide_pvalue_col, row);
              score_row.score_protein_experiment_wide_qvalue = parquetOptionalDouble_(score_protein_experiment_wide_qvalue_col, row);
              score_row.score_protein_experiment_wide_pep = parquetOptionalDouble_(score_protein_experiment_wide_pep_col, row);
              score_row.score_protein_run_specific_score = parquetOptionalDouble_(score_protein_run_specific_score_col, row);
              score_row.score_protein_run_specific_pvalue = parquetOptionalDouble_(score_protein_run_specific_pvalue_col, row);
              score_row.score_protein_run_specific_qvalue = parquetOptionalDouble_(score_protein_run_specific_qvalue_col, row);
              score_row.score_protein_run_specific_pep = parquetOptionalDouble_(score_protein_run_specific_pep_col, row);
              score_row.score_gene_global_score = parquetOptionalDouble_(score_gene_global_score_col, row);
              score_row.score_gene_global_pvalue = parquetOptionalDouble_(score_gene_global_pvalue_col, row);
              score_row.score_gene_global_qvalue = parquetOptionalDouble_(score_gene_global_qvalue_col, row);
              score_row.score_gene_global_pep = parquetOptionalDouble_(score_gene_global_pep_col, row);
              score_row.score_gene_experiment_wide_score = parquetOptionalDouble_(score_gene_experiment_wide_score_col, row);
              score_row.score_gene_experiment_wide_pvalue = parquetOptionalDouble_(score_gene_experiment_wide_pvalue_col, row);
              score_row.score_gene_experiment_wide_qvalue = parquetOptionalDouble_(score_gene_experiment_wide_qvalue_col, row);
              score_row.score_gene_experiment_wide_pep = parquetOptionalDouble_(score_gene_experiment_wide_pep_col, row);
              score_row.score_gene_run_specific_score = parquetOptionalDouble_(score_gene_run_specific_score_col, row);
              score_row.score_gene_run_specific_pvalue = parquetOptionalDouble_(score_gene_run_specific_pvalue_col, row);
              score_row.score_gene_run_specific_qvalue = parquetOptionalDouble_(score_gene_run_specific_qvalue_col, row);
              score_row.score_gene_run_specific_pep = parquetOptionalDouble_(score_gene_run_specific_pep_col, row);

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

    std::stable_sort(table.rows.begin(), table.rows.end(), [](const OpenSwathFeatureScoreRow& lhs, const OpenSwathFeatureScoreRow& rhs) {
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
    if (! config.include_transition_data) { return table; }

    std::unordered_map<Int64, std::vector<FeatureTransitionObservation_>> observations_by_transition;
    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    bool discovered_dynamic_columns = false;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string feature_transition_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/feature_transition.parquet";
      if (! File::exists(feature_transition_path)) { continue; }

      auto feature_transition_table = ParquetFile::readTable(feature_transition_path);
      if (! discovered_dynamic_columns)
      {
        for (const auto& name : featureTransitionParquetFields_())
        {
          if (ParquetFile::getOptionalColumn(feature_transition_table, name) != nullptr) { table.feature_transition_column_names.emplace_back(name); }
        }
        discovered_dynamic_columns = true;
      }

      const auto feature_id_col = ParquetFile::getColumn(feature_transition_table, "feature_id");
      const auto transition_id_col = ParquetFile::getColumn(feature_transition_table, "transition_id");
      std::unordered_map<std::string, std::shared_ptr<arrow::Array>> transition_columns;
      const auto add_transition_column
        = [&](const std::string& name) { transition_columns.emplace(name, ParquetFile::getOptionalColumn(feature_transition_table, name)); };
      for (const auto& name : table.feature_transition_column_names)
      {
        add_transition_column(name);
      }
      for (const auto& name : std::array<std::string, 5> {"score_transition_score", "score_transition_rank", "score_transition_pvalue",
                                                          "score_transition_qvalue", "score_transition_pep"})
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
          if (column != nullptr && ! column->IsNull(row)) { observation.values.push_back(ParquetFile::getDouble(column, row, 0.0, false)); }
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
      if (config.filters.exclude_decoys && transition.decoy) { continue; }

      std::vector<std::optional<Int64>> peptide_ids {std::nullopt};
      if (! transition.peptide_ids.empty())
      {
        peptide_ids.clear();
        peptide_ids.reserve(transition.peptide_ids.size());
        for (const Int64 peptide_id : transition.peptide_ids)
        {
          peptide_ids.push_back(peptide_id);
        }
      }

      std::vector<Int64> precursor_ids = transition.precursor_ids;
      if (precursor_ids.empty()) { precursor_ids.push_back(-1); }

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

    std::stable_sort(table.rows.begin(), table.rows.end(), [](const OpenSwathTransitionScoreRow& lhs, const OpenSwathTransitionScoreRow& rhs) {
      if (lhs.precursor_id != rhs.precursor_id) return lhs.precursor_id < rhs.precursor_id;
      if (lhs.transition_id != rhs.transition_id) return lhs.transition_id < rhs.transition_id;
      return lhs.feature_id.value_or(std::numeric_limits<Int64>::max()) < rhs.feature_id.value_or(std::numeric_limits<Int64>::max());
    });

    OPENMS_LOG_INFO << "Read " << table.rows.size() << " transition score rows." << std::endl;
    return table;
  }

  static bool oswpqHasIPFColumns_(const OSWPQWorkspace& workspace)
  {
    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      auto feature_table = getOSWPQFeatureTable_(workspace, run_id);
      if (ParquetFile::getOptionalColumn(feature_table, "score_ipf_qvalue") != nullptr
          || ParquetFile::getOptionalColumn(feature_table, "score_ipf_pep") != nullptr)
      {
        return true;
      }
    }
    return false;
  }

  std::vector<OpenSwathExportRow>
  readOSWPQExportRows_(const OSWPQWorkspace& workspace, const PreparedLibraryLookup_& lookup, const OpenSwathExportFilterConfig& config) const
  {
    if (config.use_alignment)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ export does not support alignment recovery because no alignment parquet tables are written yet.");
    }
    if (config.ipf_mode != OpenSwathIPFExportMode::Disable && oswpqHasIPFColumns_(workspace))
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    "Direct OSWPQ export currently supports standard OpenSWATH-style exports only. Use Export:*:ipf disable or "
                                    "workflow:working_format sqlite for IPF-aware export.");
    }

    const bool has_peptide_scores = File::exists(inferenceParquetPath_(workspace, InferenceLevel::Peptide));
    const bool has_protein_scores = File::exists(inferenceParquetPath_(workspace, InferenceLevel::Protein));
    const bool has_gene_scores = File::exists(inferenceParquetPath_(workspace, InferenceLevel::Gene));
    if (config.peptide && ! has_peptide_scores)
    {
      OPENMS_LOG_INFO
        << "Peptide-score export filtering requested, but no peptide inference table is present; leaving peptide filtering disabled for this export."
        << std::endl;
    }
    if (config.protein && ! has_protein_scores)
    {
      OPENMS_LOG_INFO
        << "Protein-score export filtering requested, but no protein inference table is present; leaving protein filtering disabled for this export."
        << std::endl;
    }
    if (config.gene && ! has_gene_scores)
    {
      OPENMS_LOG_INFO
        << "Gene-score export filtering requested, but no gene inference table is present; leaving gene filtering disabled for this export."
        << std::endl;
    }
    const auto peptide_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Peptide);
    const auto protein_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Protein);
    const auto gene_results = readLevelContextResultsParquet_(workspace, InferenceLevel::Gene);
    const ExportQValueMaps_ peptide_qvalues = buildPeptideQValueMaps_(peptide_results);
    const ExportQValueMaps_ protein_qvalues = buildAggregatedEntityQValueMaps_(protein_results, lookup.peptide_to_proteins);
    const ExportQValueMaps_ gene_qvalues = buildAggregatedEntityQValueMaps_(gene_results, lookup.peptide_to_genes);
    const auto transition_aggregations = config.transition_quantification ? buildTransitionAggregations_(workspace, lookup, config.max_transition_pep)
                                                                          : std::unordered_map<Int64, TransitionAggregation_> {};

    auto runs_table = getOSWPQRunsTable_(workspace);
    const auto run_id_col = ParquetFile::getColumn(runs_table, "run_id");
    const auto filename_col = ParquetFile::getOptionalColumn(runs_table, "filename");
    std::vector<OpenSwathExportRow> rows;

    for (int64_t run_row = 0; run_row < runs_table->num_rows(); ++run_row)
    {
      const Int64 run_id = ParquetFile::getInt64(run_id_col, run_row, 0, false);
      const std::string filename = (filename_col != nullptr && ! filename_col->IsNull(run_row)) ? ParquetFile::getString(filename_col, run_row) : "";
      const std::string features_path = workspace.base_dir + "/runs/run_id=" + StringUtils::toStr(run_id) + "/features.parquet";
      auto features_table = getOSWPQFeatureTable_(workspace, run_id);
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
        if (score_ms2_qvalue_col->IsNull(row)) { continue; }
        const double ms2_qvalue = ParquetFile::getDouble(score_ms2_qvalue_col, row, 1.0, false);
        if (! (ms2_qvalue < config.max_rs_peakgroup_qvalue)) { continue; }

        const Int64 precursor_id = ParquetFile::getInt64(precursor_id_col, row, 0, false);
        const auto precursor_it = lookup.precursors.find(precursor_id);
        if (precursor_it == lookup.precursors.end())
        {
          throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Missing prepared-library precursor metadata for precursor_id=" + StringUtils::toStr(precursor_id));
        }
        const auto peptide_mapping_it = lookup.precursor_to_peptides.find(precursor_id);
        if (peptide_mapping_it == lookup.precursor_to_peptides.end()) { continue; }

        for (const Int64 peptide_id : peptide_mapping_it->second)
        {
          const auto peptide_it = lookup.peptides.find(peptide_id);
          if (peptide_it == lookup.peptides.end()) { continue; }

          OpenSwathExportRow export_row;
          export_row.run_id = run_id;
          export_row.filename = filename;
          export_row.run_name = File::stemName(filename);
          if (export_row.run_name.empty()) { export_row.run_name = File::basename(filename); }
          if (export_row.run_name.empty()) { export_row.run_name = "RUN_ID " + StringUtils::toStr(run_id); }

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
          if (peptide_experiment_it != peptide_qvalues.experiment_wide.end())
            export_row.peptide_experiment_wide_qvalue = peptide_experiment_it->second;
          const auto peptide_run_it = peptide_qvalues.run_specific.find({run_id, peptide_id});
          if (peptide_run_it != peptide_qvalues.run_specific.end()) export_row.peptide_run_specific_qvalue = peptide_run_it->second;

          const auto protein_global_it = protein_qvalues.global.find(peptide_id);
          if (protein_global_it != protein_qvalues.global.end()) export_row.protein_global_qvalue = protein_global_it->second;
          const auto protein_experiment_it = protein_qvalues.experiment_wide.find({run_id, peptide_id});
          if (protein_experiment_it != protein_qvalues.experiment_wide.end())
            export_row.protein_experiment_wide_qvalue = protein_experiment_it->second;
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
    if (config.peptide && has_peptide_scores)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(),
                                [&](const auto& row) {
                                  return ! row.peptide_global_qvalue.has_value() || *row.peptide_global_qvalue >= config.max_global_peptide_qvalue;
                                }),
                 rows.end());
    }
    if (config.protein && has_protein_scores)
    {
      rows.erase(std::remove_if(rows.begin(), rows.end(),
                                [&](const auto& row) {
                                  return ! row.protein_global_qvalue.has_value() || *row.protein_global_qvalue >= config.max_global_protein_qvalue;
                                }),
                 rows.end());
    }
    if (config.gene && has_gene_scores)
    {
      rows.erase(std::remove_if(
                   rows.begin(), rows.end(),
                   [&](const auto& row) { return ! row.gene_global_qvalue.has_value() || *row.gene_global_qvalue >= config.max_global_gene_qvalue; }),
                 rows.end());
    }

    std::stable_sort(rows.begin(), rows.end(), [](const OpenSwathExportRow& lhs, const OpenSwathExportRow& rhs) {
      if (lhs.precursor_id != rhs.precursor_id) return lhs.precursor_id < rhs.precursor_id;
      if (lhs.feature_id != rhs.feature_id) return lhs.feature_id < rhs.feature_id;
      return lhs.peptide_id < rhs.peptide_id;
    });

    OPENMS_LOG_INFO << "Read " << rows.size() << " filtered export rows." << std::endl;
    return rows;
  }


  OSWPQWorkspace workspace_;
  PreparedLibraryLookup_ lookup_;
};

OSWParquetFile::OSWParquetFile(const std::string& filename, const OpenSwath::LightTargetedExperiment& library, bool load_transition_metadata):
    impl_(std::make_unique<Impl>())
{
  impl_->workspace_ = impl_->prepareOSWPQWorkspace_(filename);
  impl_->lookup_ = impl_->buildPreparedLibraryLookupFromLightTargetedExperiment_(library, load_transition_metadata);
}

OSWParquetFile::~OSWParquetFile() = default;
OSWParquetFile::OSWParquetFile(OSWParquetFile&&) = default;
OSWParquetFile& OSWParquetFile::operator=(OSWParquetFile&&) = default;

std::vector<LevelContextInputRow> OSWParquetFile::readLevelContextData(InferenceLevel level, InferenceContext context) const
{ return impl_->buildOSWPQLevelContextInputRows_(impl_->workspace_, impl_->lookup_, level, context); }

void OSWParquetFile::writeLevelContextResults(std::vector<InferenceResults> results)
{
  if (results.empty()) { return; }

  std::vector<Impl::InferenceTask> tasks;
  std::map<InferenceLevel, std::vector<LevelContextResultRow>> results_by_level;
  std::set<std::pair<InferenceLevel, InferenceContext>> seen;
  // Validate the entire request before modifying the workspace.
  for (auto& result : results)
  {
    Impl::entityIdColumnName_(result.level); // rejects unsupported levels
    if (! seen.emplace(result.level, result.context).second)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Duplicate OSWPQ inference level/context.", toString(result.level));
    }
    for (const auto& row : result.rows)
    {
      if (row.context != result.context)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "OSWPQ inference row has a different context from its batch.",
                                      toString(row.context));
      }
    }
    tasks.push_back({result.level, result.context});
    auto& rows = results_by_level[result.level];
    if (rows.empty()) { rows = std::move(result.rows); }
    else
    {
      rows.insert(rows.end(), std::make_move_iterator(result.rows.begin()), std::make_move_iterator(result.rows.end()));
      std::vector<LevelContextResultRow>().swap(result.rows);
    }
  }

  // A caller can explicitly preserve completed writes if a later update fails.
  impl_->workspace_.dirty = true;
  for (const auto& [level, rows] : results_by_level)
  {
    Impl::writeLevelContextResultsParquet_(impl_->workspace_, level, rows);
  }
  impl_->applyLevelContextResultsToOSWPQ_(impl_->workspace_, impl_->lookup_, results_by_level, tasks, getLogType());
}

std::vector<OpenSwathExportRow> OSWParquetFile::readOpenSwathExportRows(const OpenSwathExportFilterConfig& config) const
{ return impl_->readOSWPQExportRows_(impl_->workspace_, impl_->lookup_, config); }

OpenSwathFeatureScoreTable OSWParquetFile::readOpenSwathFeatureScoreTable(const OpenSwathParquetExportConfig& config) const
{ return impl_->readOSWPQFeatureScoreTable_(impl_->workspace_, impl_->lookup_, config); }

OpenSwathTransitionScoreTable OSWParquetFile::readOpenSwathTransitionScoreTable(const OpenSwathParquetExportConfig& config) const
{ return impl_->readOSWPQTransitionScoreTable_(impl_->workspace_, impl_->lookup_, config); }

std::map<Int64, std::string> OSWParquetFile::readRunBasenames() const
{ return Impl::readOSWPQRunBasenames_(impl_->workspace_); }

void OSWParquetFile::commit()
{ impl_->commitOSWPQWorkspace_(impl_->workspace_); }
} // namespace OpenMS
