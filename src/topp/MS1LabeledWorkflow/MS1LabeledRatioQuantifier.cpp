// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledRatioQuantifier.h"

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <iterator>
#include <cmath>
#include <limits>
#include <set>
#include <tuple>

using namespace std;

namespace OpenMS
{
  MS1LabeledRatioQuantifier::MS1LabeledRatioQuantifier() :
    DefaultParamHandler("MS1LabeledRatioQuantifier")
  {
    defaults_.setValue("reference_channel", 1, "Channel the ratios are formed against, as the 'Label' of the experimental design (1 = the light channel of '-labels').");
    defaults_.setMinInt("reference_channel", 1);
    defaults_.setValue("min_ratio_count", 2, "Minimum number of peptide ratios a protein group needs before a ratio is reported for it (MaxQuant's 'min. ratio count'). Groups below it are reported without a ratio, not with a less certain one.");
    defaults_.setMinInt("min_ratio_count", 1);
    defaults_.setValue("normalize", "true", "Additionally report every ratio divided by the median peptide ratio of its (fraction group, channel), i.e. assuming that most peptides do not change. The unnormalized ratios are reported either way.");
    defaults_.setValidStrings("normalize", {"true", "false"});
    defaultsToParam_();
  }

  double MS1LabeledRatioQuantifier::median_(std::vector<double>& values)
  {
    if (values.empty()) { return std::numeric_limits<double>::quiet_NaN(); }
    const Size middle = values.size() / 2;
    std::nth_element(values.begin(), values.begin() + middle, values.end());
    const double upper = values[middle];
    if (values.size() % 2 == 1) { return upper; }
    // even count: the mean of the two central values, so the median of two ratios is their mean
    const double lower = *std::max_element(values.begin(), values.begin() + middle);
    return (lower + upper) / 2.0;
  }

  const std::map<AASequence, std::vector<MS1LabeledRatioQuantifier::ChannelRatio>>& MS1LabeledRatioQuantifier::getPeptideRatios() const
  {
    return peptide_ratios_;
  }

  const std::map<std::string, std::vector<MS1LabeledRatioQuantifier::ChannelRatio>>& MS1LabeledRatioQuantifier::getProteinGroupRatios() const
  {
    return protein_group_ratios_;
  }

  void MS1LabeledRatioQuantifier::run(ConsensusMap& consensus, const ExperimentalDesign& design, ProteinIdentification& proteins)
  {
    peptide_ratios_.clear();
    protein_group_ratios_.clear();

    const unsigned reference_channel = static_cast<unsigned>(static_cast<int>(param_.getValue("reference_channel")));
    const Size min_ratio_count = static_cast<Size>(static_cast<int>(param_.getValue("min_ratio_count")));
    const bool normalize = param_.getValue("normalize").toBool();

    // -- (run, channel) of every column, and the fraction group of every run --
    const std::string& experiment_type = consensus.getExperimentType();
    std::map<UInt64, std::pair<std::string, unsigned>> column_to_run_channel;
    for (const auto& [map_index, header] : consensus.getColumnHeaders())
    {
      column_to_run_channel[map_index] = {File::basename(header.filename), header.getLabelAsUInt(experiment_type)};
    }
    std::map<std::string, unsigned> run_to_fraction_group;
    for (const auto& entry : design.getMSFileSection())
    {
      run_to_fraction_group[File::basename(entry.path)] = entry.fraction_group;
    }

    // -- evidence ratios: one per (feature, run, channel) --
    // Collected per peptide identity as well, so that the peptide medians below need no second pass.
    std::map<std::tuple<AASequence, unsigned, unsigned>, std::vector<double>> peptide_evidence_ratios;
    Size features_with_ratio = 0;
    Size features_without_reference = 0;

    for (ConsensusFeature& feature : consensus)
    {
      // intensity of every (run, channel) of this feature
      std::map<std::string, std::map<unsigned, double>> intensities;
      for (const FeatureHandle& handle : feature.getFeatures())
      {
        const auto column = column_to_run_channel.find(handle.getMapIndex());
        if (column == column_to_run_channel.end()) { continue; }
        intensities[column->second.first][column->second.second] += handle.getIntensity();
      }

      StringList ratio_runs;
      IntList ratio_channels;
      DoubleList ratio_values;
      bool reference_missing = false;

      for (const auto& [run, channels] : intensities)
      {
        const auto reference = channels.find(reference_channel);
        // A channel is a measurement only when it is positive and finite: a zero-intensity dummy
        // feature says the peptide is absent and a NaN one says it is not quantifiable; neither
        // yields a ratio (and an absent reference would yield infinity).
        if (reference == channels.end() || !(reference->second > 0.0) || !std::isfinite(reference->second))
        {
          reference_missing = true;
          continue;
        }
        IntList run_channels;
        DoubleList run_ratios;
        for (const auto& [channel, intensity] : channels)
        {
          if (channel == reference_channel) { continue; }
          if (!(intensity > 0.0) || !std::isfinite(intensity)) { continue; }
          run_channels.push_back(static_cast<Int>(channel));
          run_ratios.push_back(intensity / reference->second);
        }
        // The reference channel is reported as the 1.0 it is by construction, so that every ratio
        // annotation covers the complete set of channels -- but only where a channel was actually
        // measured against it: on its own, a multiplet seen in the reference channel alone says
        // nothing about a ratio, and a lone 1.0 would claim otherwise.
        if (run_ratios.empty()) { continue; }
        ratio_runs.push_back(run);
        ratio_channels.push_back(static_cast<Int>(reference_channel));
        ratio_values.push_back(1.0);
        for (Size i = 0; i < run_ratios.size(); ++i)
        {
          ratio_runs.push_back(run);
          ratio_channels.push_back(run_channels[i]);
          ratio_values.push_back(run_ratios[i]);
        }
      }

      if (ratio_values.empty())
      {
        if (reference_missing) { ++features_without_reference; }
        continue;
      }
      ++features_with_ratio;
      feature.setMetaValue("evidence_ratio_run", ratio_runs);
      feature.setMetaValue("evidence_ratio_channel", ratio_channels);
      feature.setMetaValue("evidence_ratio", ratio_values);

      // the peptide identity this evidence belongs to (unlabeled; see MS1LabelState)
      const auto& ids = feature.getPeptideIdentifications();
      if (ids.empty() || ids[0].getHits().empty()) { continue; }
      const AASequence& sequence = ids[0].getHits()[0].getSequence();
      for (Size i = 0; i < ratio_values.size(); ++i)
      {
        const auto fraction_group = run_to_fraction_group.find(ratio_runs[i]);
        if (fraction_group == run_to_fraction_group.end()) { continue; } // run not in the design
        peptide_evidence_ratios[{sequence, fraction_group->second, static_cast<unsigned>(ratio_channels[i])}].push_back(ratio_values[i]);
      }
    }

    // -- peptide ratios: median over the evidences of one peptide in one fraction group --
    std::map<std::pair<unsigned, unsigned>, std::vector<double>> ratios_per_cell; // for the normalization factors
    for (auto& [key, ratios] : peptide_evidence_ratios)
    {
      const auto& [sequence, fraction_group, channel] = key;
      ChannelRatio peptide_ratio;
      peptide_ratio.fraction_group = fraction_group;
      peptide_ratio.channel = channel;
      peptide_ratio.count = ratios.size();
      peptide_ratio.ratio = median_(ratios);
      peptide_ratios_[sequence].push_back(peptide_ratio);
      ratios_per_cell[{fraction_group, channel}].push_back(peptide_ratio.ratio);
    }

    // -- normalization factor per (fraction group, channel): the median peptide ratio --
    std::map<std::pair<unsigned, unsigned>, double> normalization;
    for (auto& [cell, ratios] : ratios_per_cell)
    {
      const double factor = median_(ratios);
      normalization[cell] = (normalize && std::isfinite(factor) && factor > 0.0) ? factor : 1.0;
    }
    const auto normalized = [&normalization](const ChannelRatio& r)
    {
      const auto factor = normalization.find({r.fraction_group, r.channel});
      return factor == normalization.end() ? r.ratio : r.ratio / factor->second;
    };
    for (auto& [sequence, ratios] : peptide_ratios_)
    {
      for (ChannelRatio& r : ratios) { r.normalized_ratio = normalized(r); }
    }

    // -- the peptide ratios back onto the features of that peptide --
    // A feature is one evidence, but the peptide ratio is the quantity MaxQuant reports at peptide
    // level, so it travels with every feature the peptide was measured in (and reaches the mzTab
    // peptide rows, which are per feature).
    for (ConsensusFeature& feature : consensus)
    {
      const auto& ids = feature.getPeptideIdentifications();
      if (ids.empty() || ids[0].getHits().empty()) { continue; }
      const auto ratios = peptide_ratios_.find(ids[0].getHits()[0].getSequence());
      if (ratios == peptide_ratios_.end()) { continue; }

      IntList fraction_groups, channels, counts;
      DoubleList values, normalized_values;
      for (const ChannelRatio& r : ratios->second)
      {
        fraction_groups.push_back(static_cast<Int>(r.fraction_group));
        channels.push_back(static_cast<Int>(r.channel));
        counts.push_back(static_cast<Int>(r.count));
        values.push_back(r.ratio);
        normalized_values.push_back(r.normalized_ratio);
      }
      feature.setMetaValue("peptide_ratio_fraction_group", fraction_groups);
      feature.setMetaValue("peptide_ratio_channel", channels);
      feature.setMetaValue("peptide_ratio", values);
      if (normalize) { feature.setMetaValue("peptide_ratio_normalized", normalized_values); }
      feature.setMetaValue("peptide_ratio_count", counts);
    }

    // -- protein group ratios: median over the peptide ratios of the group's own peptides --
    // A peptide belongs to a group when every protein it references is in that group: peptides
    // shared between groups are not evidence for either of them (a razor assignment, if wanted, has
    // already rewritten the references during inference).
    std::map<AASequence, std::set<std::string>> peptide_accessions;
    const auto collect = [&peptide_accessions](const PeptideIdentificationList& ids)
    {
      for (const PeptideIdentification& id : ids)
      {
        if (id.getHits().empty()) { continue; }
        const PeptideHit& hit = id.getHits()[0];
        const auto accessions = hit.extractProteinAccessionsSet();
        peptide_accessions[hit.getSequence()].insert(accessions.begin(), accessions.end());
      }
    };
    for (const ConsensusFeature& feature : consensus) { collect(feature.getPeptideIdentifications()); }
    collect(consensus.getUnassignedPeptideIdentifications());

    // Each peptide is assigned once, through the group of its first accession, and only if that group
    // covers every accession it references -- rather than testing every (group, peptide) pair.
    auto& groups = proteins.getIndistinguishableProteins();
    std::map<std::string, Size> accession_to_group;
    for (Size g = 0; g < groups.size(); ++g)
    {
      for (const std::string& accession : groups[g].accessions) { accession_to_group.emplace(accession, g); }
    }
    std::vector<std::map<std::pair<unsigned, unsigned>, std::vector<double>>> group_ratios(groups.size());
    for (const auto& [sequence, ratios] : peptide_ratios_)
    {
      const auto accessions = peptide_accessions.find(sequence);
      if (accessions == peptide_accessions.end() || accessions->second.empty()) { continue; }
      const auto candidate = accession_to_group.find(*accessions->second.begin());
      if (candidate == accession_to_group.end()) { continue; }
      const auto& group_accessions = groups[candidate->second].accessions;
      const bool all_in_group = std::all_of(accessions->second.begin(), accessions->second.end(),
                                            [&group_accessions](const std::string& accession)
                                            { return std::find(group_accessions.begin(), group_accessions.end(), accession) != group_accessions.end(); });
      if (!all_in_group) { continue; }
      for (const ChannelRatio& r : ratios) { group_ratios[candidate->second][{r.fraction_group, r.channel}].push_back(r.ratio); }
    }

    static const char* const ratio_array_names[] = {"fraction_group_level_ratio_fraction_group", "fraction_group_level_ratio_label",
                                                    "fraction_group_level_ratio_count", "fraction_group_level_ratio",
                                                    "fraction_group_level_ratio_normalized"};
    for (Size g = 0; g < groups.size(); ++g)
    {
      ProteinIdentification::ProteinGroup& group = groups[g];
      if (group.accessions.empty()) { continue; }

      // a second run() replaces the annotation instead of adding to it
      const auto is_ratio_array = [](const auto& array)
      { return std::find(std::begin(ratio_array_names), std::end(ratio_array_names), array.getName()) != std::end(ratio_array_names); };
      std::erase_if(group.getIntegerDataArrays(), is_ratio_array);
      std::erase_if(group.getFloatDataArrays(), is_ratio_array);

      ProteinIdentification::ProteinGroup::IntegerDataArray fraction_groups, labels, counts;
      fraction_groups.setName("fraction_group_level_ratio_fraction_group");
      labels.setName("fraction_group_level_ratio_label");
      counts.setName("fraction_group_level_ratio_count");
      ProteinIdentification::ProteinGroup::FloatDataArray values, normalized_values;
      values.setName("fraction_group_level_ratio");
      normalized_values.setName("fraction_group_level_ratio_normalized");

      std::vector<ChannelRatio> reported;
      for (auto& [cell, ratios] : group_ratios[g])
      {
        // Below the minimum the ratio is not reported at all: a protein quantified from one peptide
        // is exactly what the threshold exists to keep out of the result.
        if (ratios.size() < min_ratio_count) { continue; }
        ChannelRatio group_ratio;
        group_ratio.fraction_group = cell.first;
        group_ratio.channel = cell.second;
        group_ratio.count = ratios.size();
        group_ratio.ratio = median_(ratios);
        group_ratio.normalized_ratio = normalized(group_ratio);
        reported.push_back(group_ratio);

        fraction_groups.push_back(static_cast<Int>(group_ratio.fraction_group));
        labels.push_back(static_cast<Int>(group_ratio.channel));
        counts.push_back(static_cast<Int>(group_ratio.count));
        values.push_back(static_cast<float>(group_ratio.ratio));
        normalized_values.push_back(static_cast<float>(group_ratio.normalized_ratio));
      }
      if (reported.empty()) { continue; }
      protein_group_ratios_[group.accessions.front()] = reported;

      // next to the abundance arrays of PeptideAndProteinQuant, which use their own key arrays
      auto& integer_arrays = group.getIntegerDataArrays();
      integer_arrays.push_back(std::move(fraction_groups));
      integer_arrays.push_back(std::move(labels));
      integer_arrays.push_back(std::move(counts));
      auto& float_arrays = group.getFloatDataArrays();
      float_arrays.push_back(std::move(values));
      if (normalize) { float_arrays.push_back(std::move(normalized_values)); }
    }

    OPENMS_LOG_INFO << "Channel ratios: " << features_with_ratio << " multiplet(s) with a ratio, "
                    << peptide_ratios_.size() << " peptide(s) and " << protein_group_ratios_.size()
                    << " protein group(s) quantified by ratio";
    if (features_without_reference > 0)
    {
      OPENMS_LOG_INFO << "; " << features_without_reference << " multiplet(s) have no positive reference channel and no ratio";
    }
    OPENMS_LOG_INFO << "." << endl;
  }
} // namespace OpenMS
