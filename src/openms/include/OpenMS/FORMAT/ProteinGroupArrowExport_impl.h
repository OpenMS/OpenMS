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
#include <functional>
#include <map>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::Internal
{

/// One label of a quantification unit, and the sample its quantity is stored under.
struct QuantChannel
{
  unsigned label = 1;      ///< design label / channel number, 1-based
  Size sample = 0;         ///< index into the protein group's sample-level abundance array
  std::string qpx_label;   ///< QPX intensity label ("LFQ", "TMT126", ...)
};

/// One QPX quantification unit: the raw files aggregated into a single protein quantity,
/// together with the sample each of its labels resolves to.
struct QuantUnit
{
  std::vector<std::string> runs;       ///< QPX run names, in (fraction group, fraction) order
  std::vector<QuantChannel> channels;  ///< one per label the unit's runs carry
};

/// The quantification units of an experimental design, indexed by run name.
///
/// QPX 1.1 keys the pg view on the SET of files aggregated into one quantity (bigbio/qpx#220).
/// In OpenMS that set is defined by the SAMPLE, because the sample is what the quantity is
/// summed over: PeptideAndProteinQuant accumulates `total_abundances[sample] += abundance`
/// across every fraction, file, charge and channel of that sample (PeptideAndProteinQuant.cpp),
/// so the resulting per-file cells are not summands of the sample value and the files cannot be
/// separated after the fact. It also matches how the spec says a row is resolved:
/// "(any file in grouped_runs, label) -> run.samples[]".
///
/// A unit is a CONNECTED COMPONENT of the (run, sample) incidence: two runs belong together as
/// soon as they share one sample, transitively. That is the weakest grouping under which no
/// sample can straddle two units, which is what makes the spec's "any file in grouped_runs"
/// well defined and stops one aggregated quantity from being published twice.
///
/// Two rules it deliberately is NOT:
///  * the design's fraction_group -- share/OpenMS/examples/FRACTIONS/BSA_design_onetable_nonconsec.tsv
///    puts BSA1_F1 in fraction group 1 and BSA1_F2 in fraction group 4 while both are sample
///    BSA1, so keying on it would emit that one sample quantity twice;
///  * equality of a run's whole (label -> sample) map -- a design may simply be missing a row,
///    which OpenMS accepts on purpose (ExperimentalDesign_test.cpp: "missing fractions and wrong
///    orders should work now"). In ExperimentalDesign_input_2_wrong.tsv one fraction lacks its
///    label-1 row; requiring equal maps would split that plex into separate units that then
///    share samples, refusing a design the loader accepts.
///
/// Run names are File::stemName()'d design paths -- PeptideAndProteinQuant's own key
/// convention for the design -- so the two sides agree by construction.
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
      std::map<unsigned, Size> label_sample;   ///< this run's label -> sample
    };
    std::map<std::string, RunInfo> run_info;   // ordered, so the grouping below is deterministic

    for (const auto& entry : design.getMSFileSection())
    {
      const std::string run = ArrowIOHelpers::qpxRunFileName(entry.path);
      if (run.empty()) { continue; }
      if (!known_runs.count(run)) { dropped_runs_.insert(run); continue; }
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
      info.fraction_group = entry.fraction_group;
      info.fraction = entry.fraction;

      // (path, label) is unique in a well-formed design, so a repeat with a different sample
      // means the design cannot say which sample this run's label belongs to.
      auto [it, inserted] = info.label_sample.emplace(entry.label, entry.sample);
      if (!inserted && it->second != entry.sample) { inconsistent_ = true; }
    }

    // -- Connected components of the (run, sample) incidence, by union-find over the runs --
    std::vector<std::string> run_names;
    for (const auto& [run, info] : run_info) { (void)info; run_names.push_back(run); }
    std::vector<size_t> parent(run_names.size());
    for (size_t i = 0; i < parent.size(); ++i) { parent[i] = i; }
    std::function<size_t(size_t)> find = [&](size_t i) { while (parent[i] != i) { parent[i] = parent[parent[i]]; i = parent[i]; } return i; };
    std::map<Size, size_t> first_run_of_sample;
    for (size_t i = 0; i < run_names.size(); ++i)
    {
      for (const auto& [label, sample] : run_info.at(run_names[i]).label_sample)
      {
        (void)label;
        auto [it, fresh] = first_run_of_sample.emplace(sample, i);
        if (!fresh) { parent[find(i)] = find(it->second); }   // shares a sample -> same unit
      }
    }

    std::map<size_t, std::vector<std::string>> components;   // ordered by first run index
    for (size_t i = 0; i < run_names.size(); ++i) { components[find(i)].push_back(run_names[i]); }

    for (const auto& [root, runs] : components)
    {
      (void)root;
      QuantUnit unit;
      unit.runs = runs;
      // Fraction order, so the list reads as the acquisition did rather than alphabetically.
      std::sort(unit.runs.begin(), unit.runs.end(),
                [&run_info](const std::string& a, const std::string& b)
                {
                  const auto& ia = run_info.at(a);
                  const auto& ib = run_info.at(b);
                  return std::tie(ia.fraction_group, ia.fraction, a) <
                         std::tie(ib.fraction_group, ib.fraction, b);
                });

      // A unit must be RECTANGULAR: every (label -> sample) cell it carries has to exist, with
      // the same sample, in EVERY one of its runs. That is precisely the spec's resolution rule
      // "(any file in grouped_runs, label) -> run.samples[]" -- "any file" only means anything
      // when the mapping holds for all of them. Two ways it can fail, both refused because the
      // row would name files that did not contribute to the sample it publishes:
      //  * RAGGED: a run has no cell for a label the unit carries. Connectivity puts runs in one
      //    component as soon as they share ONE sample, so a design missing a row bridges runs
      //    whose remaining samples have smaller support -- the row would claim that run
      //    aggregated a sample it never saw.
      //  * AMBIGUOUS: two runs disagree about a label's sample. A common reference channel
      //    shared across plexes joins the plexes, and the other channels then disagree.
      // Plus the reverse direction: a sample reached by two labels of one unit would have its
      // single aggregated number published twice.
      // Naming a unit by its first run and a count keeps the diagnostic readable: a
      // fractionated plex has dozens of runs.
      const std::string unit_name = "the unit starting at '" + unit.runs.front() + "'"
                                  + (unit.runs.size() > 1 ? " (" + std::to_string(unit.runs.size()) + " runs)" : "");
      std::map<unsigned, Size> label_sample;
      std::map<Size, unsigned> sample_label;
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
          auto [sl, sl_fresh] = sample_label.emplace(sample, label);
          if (!sl_fresh && sl->second != label && duplicated_sample_.empty())
          {
            duplicated_sample_ = "sample " + std::to_string(sample) + " fills both label "
                               + std::to_string(sl->second) + " and label " + std::to_string(label)
                               + " of " + unit_name;
          }
        }
      }
      for (const auto& run : unit.runs)
      {
        for (const auto& [label, sample] : label_sample)
        {
          if (run_info.at(run).label_sample.count(label)) { continue; }
          if (!ragged_unit_.empty()) { continue; }
          ragged_unit_ = "run '" + run + "' of " + unit_name + " has no design row for label "
                       + std::to_string(label) + ", but the unit carries sample "
                       + std::to_string(sample) + " at that label";
        }
      }
      for (const auto& [label, sample] : label_sample) { unit.channels.push_back(QuantChannel{label, sample, ""}); }
      for (const auto& run : unit.runs) { unit_of_run_[run] = units_.size(); }
      units_.push_back(std::move(unit));
    }
  }

  bool empty() const { return units_.empty(); }
  /// A design assigning one (run, label) to several samples; the rows would be wrong.
  bool inconsistent() const { return inconsistent_; }
  /// Non-empty when one sample fills several cells of a unit; the row would double-count it.
  const std::string& duplicatedSample() const { return duplicated_sample_; }
  /// Non-empty when one label of a unit means several samples; the row cannot be written.
  const std::string& ambiguousLabel() const { return ambiguous_label_; }
  /// Non-empty when a unit's runs do not all carry every label it publishes.
  const std::string& raggedUnit() const { return ragged_unit_; }
  /// Design files the ConsensusMap has no column header for, left out of every unit.
  const std::set<std::string>& droppedRuns() const { return dropped_runs_; }
  /// Non-empty when a run that IS in the map lacks a header for a label the design declares.
  const std::string& missingCell() const { return missing_cell_; }
  const std::vector<QuantUnit>& units() const { return units_; }

  /// Index of the unit containing @p run, or units().size() when the design does not list it.
  size_t indexOf(const std::string& run) const
  {
    auto it = unit_of_run_.find(run);
    return (it == unit_of_run_.end()) ? units_.size() : it->second;
  }

  /// Fill in every channel's QPX intensity label from @p label_for(run, channel_number).
  ///
  /// A unit's runs are one sample measured several times, so they must agree on what a
  /// channel is called. They are all consulted rather than only the first: a disagreement means
  /// the fractions were labelled differently, which makes the single label the row can carry
  /// arbitrary, and that is worth saying out loud even though the export can continue.
  /// @return false if @p label_for could not name a channel; the caller must refuse the export.
  template <typename LabelFor>
  bool resolveLabels(LabelFor label_for, const std::string& context)
  {
    for (auto& unit : units_)
    {
      for (auto& channel : unit.channels)
      {
        channel.qpx_label = label_for(unit.runs.front(), channel.label);
        if (channel.qpx_label.empty()) { return false; }   // already logged
        for (size_t i = 1; i < unit.runs.size(); ++i)
        {
          const std::string other = label_for(unit.runs[i], channel.label);
          if (other == channel.qpx_label) { continue; }
          if (other.empty()) { return false; }             // already logged
          OPENMS_LOG_WARN << context << ": runs '" << unit.runs.front() << "' and '" << unit.runs[i]
                          << "' belong to one quantification unit but label channel "
                          << channel.label << " differently ('" << channel.qpx_label << "' vs '"
                          << other << "'). Keeping '" << channel.qpx_label << "'." << std::endl;
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
  std::string duplicated_sample_;
  std::string ambiguous_label_;
  std::string ragged_unit_;
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

} // namespace OpenMS::Internal
