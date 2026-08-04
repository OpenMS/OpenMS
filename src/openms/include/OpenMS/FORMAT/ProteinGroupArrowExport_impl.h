// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// INTERNAL HEADER -- deliberately not installed (not listed in the install rules for
// OpenMS/FORMAT). It exists so ProteinGroupArrowExport and QPXCollectionExport share ONE
// definition of what a QPX quantification unit is: the exporter decides what to write and the
// preflight decides whether the collection may be written at all, and the two must not be able
// to disagree about which designs are representable.

#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ArrowIOHelpers.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <algorithm>
#include <map>
#include <set>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace OpenMS::Internal
{

/// One label of a quantification unit and the sample it resolves to in the design.
struct QuantChannel
{
  unsigned label = 1;      ///< design label / channel number, 1-based
  Size sample = 0;         ///< design sample, retained to validate consistent run/sample resolution
  std::string qpx_label;   ///< QPX intensity label ("LFQ", "TMT126", ...)
};

/// One QPX quantification unit: exactly one OpenMS experimental-design fraction group.
struct QuantUnit
{
  UInt fraction_group = 0;             ///< OpenMS fraction_group defining this unit
  std::vector<std::string> runs;       ///< QPX run names, in fraction order
  std::vector<QuantChannel> channels;  ///< one per label the unit's runs carry
};

/// The quantification units of an experimental design, indexed by run name.
///
/// Active QPX 1.1 defines OpenMS @c fraction_group as @c grouped_runs: all files in one fraction
/// group form one quantity, and @c label is the orthogonal channel dimension. A sample may occur
/// in several fraction groups (technical replicates or a shared reference) without joining them.
/// PeptideAndProteinQuant preserves quantities at exactly this grain before protein aggregation.
class QuantificationUnits
{
public:
  QuantificationUnits() = default;

  /// @param design  the design that drove quantification
  /// @param known_cells the (QPX run name, label) cells the ConsensusMap has a column header
  ///        for. A design row for a file that is not in the map contributed nothing, so naming
  ///        it in grouped_runs would claim an aggregation that did not happen -- and would put a
  ///        run in the pg view that the psm and feature views cannot join to.
  ///
  ///        The set holds CELLS, not runs, because the two absences mean different things:
  ///         * a run with NO cell at all is simply not part of this map -- a design listing files
  ///           beyond the current input, which ProteomicsLFQ produces routinely by filtering the
  ///           design by basename. Dropped, with a diagnostic.
  ///         * a run that IS present but whose declared label has no header means the design and
  ///           the data disagree about the channel layout -- a TMT10-shaped design against a
  ///           6-plex map, say. Reading that cell would invent an intensity label ("LFQ" for
  ///           channel 1, otherwise the bare channel number) that joins to nothing, so it is
  ///           recorded and refused instead.
  QuantificationUnits(const ExperimentalDesign& design,
                      const std::set<std::pair<std::string, unsigned>>& known_cells)
  {
    std::set<std::string> known_runs;
    for (const auto& [run, label] : known_cells) { (void)label; known_runs.insert(run); }

    struct RunInfo
    {
      unsigned fraction_group = 0;
      unsigned fraction = 0;
      bool initialized = false;
      std::map<unsigned, Size> label_sample;   ///< this run's label -> sample
    };
    std::map<std::string, RunInfo> run_info;
    const Size n_samples = design.getNumberOfSamples();

    for (const auto& entry : design.getMSFileSection())
    {
      const std::string run = ArrowIOHelpers::qpxRunFileName(entry.path);
      if (run.empty()) { continue; }
      if (!known_runs.count(run)) { dropped_runs_.insert(run); continue; }
      if (entry.sample >= n_samples && invalid_sample_.empty())
      {
        invalid_sample_ = "run '" + run + "', label " + std::to_string(entry.label)
                        + " refers to sample " + std::to_string(entry.sample)
                        + ", but the sample section contains " + std::to_string(n_samples)
                        + " sample(s)";
      }
      if (!known_cells.count({run, entry.label}))
      {
        if (missing_cell_.empty())
        {
          missing_cell_ = "run '" + run + "' is in this map, but no column header describes its "
                          "label " + std::to_string(entry.label)
                        + ", which the experimental design declares";
        }
        continue;
      }
      auto& info = run_info[run];
      if (!info.initialized)
      {
        info.fraction_group = entry.fraction_group;
        info.fraction = entry.fraction;
        info.initialized = true;
      }
      else if (info.fraction_group != entry.fraction_group && run_in_several_groups_.empty())
      {
        run_in_several_groups_ = "run '" + run + "' belongs to fraction groups "
                               + std::to_string(info.fraction_group) + " and "
                               + std::to_string(entry.fraction_group);
      }

      // (path, label) is unique in a well-formed design, so a repeat with a different sample
      // means the design cannot say which sample this run's label belongs to.
      auto [it, inserted] = info.label_sample.emplace(entry.label, entry.sample);
      if (!inserted && it->second != entry.sample) { inconsistent_ = true; }
    }

    std::map<UInt, std::vector<std::string>> runs_by_fraction_group;
    for (const auto& [run, info] : run_info)
    {
      runs_by_fraction_group[info.fraction_group].push_back(run);
    }

    for (auto& [fraction_group, runs] : runs_by_fraction_group)
    {
      QuantUnit unit;
      unit.fraction_group = fraction_group;
      unit.runs = std::move(runs);
      std::sort(unit.runs.begin(), unit.runs.end(),
                [&run_info](const std::string& a, const std::string& b)
                {
                  const auto& ia = run_info.at(a);
                  const auto& ib = run_info.at(b);
                  return std::tie(ia.fraction, a) < std::tie(ib.fraction, b);
                });

      const std::string unit_name = "fraction group " + std::to_string(fraction_group);
      std::map<unsigned, Size> label_sample;
      for (const auto& run : unit.runs)
      {
        for (const auto& [label, sample] : run_info.at(run).label_sample)
        {
          auto [ls, ls_fresh] = label_sample.emplace(label, sample);
          if (!ls_fresh && ls->second != sample && ambiguous_label_.empty())
          {
            ambiguous_label_ = "label " + std::to_string(label) + " of " + unit_name
                             + " maps to sample " + std::to_string(ls->second) + " and sample "
                             + std::to_string(sample);
          }
        }
      }
      for (const auto& run : unit.runs)
      {
        for (const auto& [label, sample] : label_sample)
        {
          if (run_info.at(run).label_sample.count(label)) { continue; }
          if (!ragged_unit_.empty()) { continue; }
          ragged_unit_ = "run '" + run + "' of " + unit_name
                       + " has no design row for label " + std::to_string(label)
                       + ", which resolves to sample " + std::to_string(sample)
                       + " in the other run(s)";
        }
      }
      for (const auto& [label, sample] : label_sample) { unit.channels.push_back(QuantChannel{label, sample, ""}); }
      for (const auto& run : unit.runs)
      {
        auto [it, inserted] = unit_of_run_.emplace(run, units_.size());
        if (!inserted && run_in_several_groups_.empty())
        {
          run_in_several_groups_ = "run '" + run + "' belongs to more than one fraction group";
        }
      }
      units_.push_back(std::move(unit));
    }
  }

  bool empty() const { return units_.empty(); }
  /// A design assigning one (run, label) to several samples; the rows would be wrong.
  bool inconsistent() const { return inconsistent_; }
  /// Non-empty when one label of a unit means several samples; the row cannot be written.
  const std::string& ambiguousLabel() const { return ambiguous_label_; }
  /// Non-empty when a unit's runs do not all carry every label it publishes.
  const std::string& raggedUnit() const { return ragged_unit_; }
  /// Design files the ConsensusMap has no column header for, left out of every unit.
  const std::set<std::string>& droppedRuns() const { return dropped_runs_; }
  /// Non-empty when a run that IS in the map lacks a header for a label the design declares.
  const std::string& missingCell() const { return missing_cell_; }
  /// Non-empty when a design row refers outside its sample section.
  const std::string& invalidSample() const { return invalid_sample_; }
  /// Non-empty when one run is assigned to multiple fraction groups.
  const std::string& runInSeveralGroups() const { return run_in_several_groups_; }
  const std::vector<QuantUnit>& units() const { return units_; }

  /// Index of the unit containing @p run, or units().size() when the design does not list it.
  size_t indexOf(const std::string& run) const
  {
    auto it = unit_of_run_.find(run);
    return (it == unit_of_run_.end()) ? units_.size() : it->second;
  }

  /// Fill in every channel's QPX intensity label from @p label_for(run, channel_number).
  ///
  /// Every run in a fraction group must agree on the canonical QPX label for a numeric channel.
  /// @return false if a channel cannot be named, two runs disagree, or two channels resolve to
  /// the same QPX label; the caller must refuse.
  template <typename LabelFor>
  bool resolveLabels(LabelFor label_for, const std::string& context)
  {
    for (auto& unit : units_)
    {
      std::map<std::string, unsigned> label_to_channel;
      for (auto& channel : unit.channels)
      {
        channel.qpx_label = label_for(unit.runs.front(), channel.label);
        if (channel.qpx_label.empty()) { return false; }   // already logged
        for (size_t i = 1; i < unit.runs.size(); ++i)
        {
          const std::string other = label_for(unit.runs[i], channel.label);
          if (other == channel.qpx_label) { continue; }
          if (other.empty()) { return false; }             // already logged
          OPENMS_LOG_ERROR << context << ": runs '" << unit.runs.front() << "' and '"
                           << unit.runs[i] << "' belong to fraction group " << unit.fraction_group
                           << " but name label " << channel.label << " differently ('"
                           << channel.qpx_label << "' vs '" << other
                           << "'). A grouped row cannot choose one arbitrarily." << std::endl;
          return false;
        }
        auto [label_it, label_inserted] = label_to_channel.emplace(channel.qpx_label, channel.label);
        if (!label_inserted && label_it->second != channel.label)
        {
          OPENMS_LOG_ERROR << context << ": channels " << label_it->second << " and "
                           << channel.label << " of fraction group " << unit.fraction_group
                           << " both resolve to QPX label '" << channel.qpx_label
                           << "'. They would duplicate the protein-group primary key." << std::endl;
          return false;
        }
      }
    }
    return true;
  }

private:
  std::vector<QuantUnit> units_;
  std::unordered_map<std::string, size_t> unit_of_run_;
  std::set<std::string> dropped_runs_;
  std::string missing_cell_;
  std::string ambiguous_label_;
  std::string ragged_unit_;
  std::string invalid_sample_;
  std::string run_in_several_groups_;
  bool inconsistent_ = false;
};

/// The group's sample-level abundance array, or nullptr when it carries no quantification.
///
/// Looked up by NAME rather than by index: PeptideAndProteinQuant names it "abundances"
/// (PeptideAndProteinQuant.cpp) and the positional convention -- four float arrays in a fixed
/// order -- is an implementation detail of that one writer.
inline const ProteinIdentification::ProteinGroup::FloatDataArray* sampleAbundances(
  const ProteinIdentification::ProteinGroup& group)
{
  for (const auto& array : group.getFloatDataArrays())
  {
    if (array.getName() == "abundances") { return &array; }
  }
  return nullptr;
}

/// Parsed parallel arrays carrying protein quantities at (fraction_group, label) grain.
struct FractionGroupAbundanceData
{
  bool present = false;
  bool valid = true;
  std::string error;
  std::map<std::pair<UInt, UInt>, float> values;
};

/// Read PeptideAndProteinQuant's named unit-level arrays and validate their parallel-key contract.
inline FractionGroupAbundanceData fractionGroupAbundances(
  const ProteinIdentification::ProteinGroup& group)
{
  constexpr const char* abundance_name = "fraction_group_level_abundance";
  constexpr const char* fraction_group_name = "fraction_group_level_fraction_group";
  constexpr const char* label_name = "fraction_group_level_label";

  const ProteinIdentification::ProteinGroup::FloatDataArray* abundances = nullptr;
  const ProteinIdentification::ProteinGroup::IntegerDataArray* fraction_groups = nullptr;
  const ProteinIdentification::ProteinGroup::IntegerDataArray* labels = nullptr;
  for (const auto& array : group.getFloatDataArrays())
  {
    if (array.getName() == abundance_name) { abundances = &array; }
  }
  for (const auto& array : group.getIntegerDataArrays())
  {
    if (array.getName() == fraction_group_name) { fraction_groups = &array; }
    else if (array.getName() == label_name) { labels = &array; }
  }

  FractionGroupAbundanceData result;
  result.present = abundances != nullptr || fraction_groups != nullptr || labels != nullptr;
  if (!result.present) { return result; }
  if (abundances == nullptr || fraction_groups == nullptr || labels == nullptr)
  {
    result.valid = false;
    result.error = "the fraction-group abundance, fraction-group key and label key arrays are not all present";
    return result;
  }
  if (abundances->size() != fraction_groups->size() || abundances->size() != labels->size())
  {
    result.valid = false;
    result.error = "the parallel fraction-group abundance arrays have different lengths";
    return result;
  }

  for (Size i = 0; i < abundances->size(); ++i)
  {
    if ((*fraction_groups)[i] < 0 || (*labels)[i] <= 0)
    {
      result.valid = false;
      result.error = "a fraction-group key is negative or a label key is not 1-based";
      return result;
    }
    const auto key = std::make_pair(static_cast<UInt>((*fraction_groups)[i]),
                                    static_cast<UInt>((*labels)[i]));
    if (!result.values.emplace(key, (*abundances)[i]).second)
    {
      result.valid = false;
      result.error = "the fraction-group abundance arrays contain a duplicate (fraction_group, label) key";
      return result;
    }
  }
  return result;
}

/**
  @brief Validate one protein group's complete QPX 1.1 quantification annotation contract.

  @param[in] sample_abundances Optional legacy sample-level abundance array
  @param[in] quantities Parsed fraction-group/label abundance arrays
  @param[in] design_samples Number of samples in the experimental design
  @param[in] units Quantification units derived from the design and ConsensusMap headers
  @param[in] expected_quantity_keys Fraction-group/label keys required by @p units
  @return An empty string when exportable, otherwise a diagnostic describing the first violation
*/
inline std::string validateQuantificationAnnotations(
  const ProteinIdentification::ProteinGroup::FloatDataArray* sample_abundances,
  const FractionGroupAbundanceData& quantities,
  Size design_samples,
  const QuantificationUnits& units,
  const std::set<std::pair<UInt, UInt>>& expected_quantity_keys)
{
  if (!quantities.valid)
  {
    return "has invalid fraction-group abundance annotations: " + quantities.error;
  }
  if (sample_abundances != nullptr && sample_abundances->size() != design_samples)
  {
    return "carries " + std::to_string(sample_abundances->size())
         + " sample abundances but the experimental design describes "
         + std::to_string(design_samples)
         + " samples, so the design does not belong to this quantification";
  }

  const bool quantified = sample_abundances != nullptr || quantities.present;
  if (!quantified) { return {}; }
  if (!quantities.present)
  {
    return "carries sample abundances but no fraction-group/label abundances from the QPX 1.1 "
           "quantification path";
  }
  if (units.empty())
  {
    return "carries fraction-group abundances, but the experimental design yields no matching "
           "quantification unit";
  }
  if (quantities.values.size() != expected_quantity_keys.size())
  {
    return "carries " + std::to_string(quantities.values.size())
         + " fraction-group/label quantities but the experimental design requires "
         + std::to_string(expected_quantity_keys.size())
         + ", so the annotations and design do not describe the same quantification";
  }
  for (const auto& key : expected_quantity_keys)
  {
    if (!quantities.values.count(key))
    {
      return "has no quantity for fraction group " + std::to_string(key.first)
           + ", label " + std::to_string(key.second);
    }
  }
  return {};
}

} // namespace OpenMS::Internal
