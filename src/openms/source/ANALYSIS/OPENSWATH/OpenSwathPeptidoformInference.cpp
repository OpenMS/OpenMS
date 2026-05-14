// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <tuple>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  namespace
  {
    struct TransitionInferenceRow
    {
      Int64 feature_id = -1;
      Int64 transition_id = -1;
      Int64 hypothesis = -1;
      double pep = 1.0;
      double precursor_peakgroup_pep = 1.0;
      Int32 bmask = 0;
      Int32 num_peptidoforms = 0;
      Int64 alignment_group_id = -1;
    };

    std::vector<TransitionInferenceRow> propagateSignalAcrossRuns_(const std::vector<TransitionInferenceRow>& rows,
                                                                   const double confidence_threshold)
    {
      std::vector<TransitionInferenceRow> result;
      result.reserve(rows.size());

      std::unordered_map<Int64, std::vector<const TransitionInferenceRow*>> grouped_rows;
      grouped_rows.reserve(rows.size());
      for (const auto& row : rows)
      {
        // Negative group ids are reserved for rows that do not belong to any
        // across-run alignment group. Using a dedicated sentinel avoids clashes
        // with real DENSE_RANK()-based group ids on small synthetic OSWs where
        // group ids and feature ids can numerically overlap.
        if (row.alignment_group_id < 0)
        {
          result.push_back(row);
        }
        else
        {
          grouped_rows[row.alignment_group_id].push_back(&row);
        }
      }

      for (const auto& group : grouped_rows)
      {
        std::vector<Int64> feature_ids;
        feature_ids.reserve(group.second.size());
        for (const TransitionInferenceRow* row_ptr : group.second)
        {
          feature_ids.push_back(row_ptr->feature_id);
        }
        std::sort(feature_ids.begin(), feature_ids.end());
        feature_ids.erase(std::unique(feature_ids.begin(), feature_ids.end()), feature_ids.end());

        std::vector<TransitionInferenceRow> propagated_rows;
        for (Int64 target_feature_id : feature_ids)
        {
          for (const TransitionInferenceRow* row_ptr : group.second)
          {
            // PyProphet-style propagation keeps a feature's native evidence and adds
            // evidence imported from aligned features only when that source evidence
            // is already confident enough.
            if (row_ptr->feature_id == target_feature_id || row_ptr->pep <= confidence_threshold)
            {
              TransitionInferenceRow propagated = *row_ptr;
              propagated.feature_id = target_feature_id;
              propagated_rows.push_back(std::move(propagated));
            }
          }
        }

        std::sort(propagated_rows.begin(), propagated_rows.end(),
                  [](const TransitionInferenceRow& lhs, const TransitionInferenceRow& rhs)
                  {
                    return std::tie(lhs.feature_id, lhs.transition_id, lhs.hypothesis, lhs.bmask, lhs.num_peptidoforms, lhs.alignment_group_id) <
                           std::tie(rhs.feature_id, rhs.transition_id, rhs.hypothesis, rhs.bmask, rhs.num_peptidoforms, rhs.alignment_group_id);
                  });

        // Collapse propagated duplicates immediately so the later Bayesian model only sees the
        // compact minimum-PEP evidence per propagated key.
        for (const auto& row : propagated_rows)
        {
          if (!result.empty())
          {
            TransitionInferenceRow& back = result.back();
            if (std::tie(back.feature_id, back.transition_id, back.hypothesis, back.bmask, back.num_peptidoforms, back.alignment_group_id) ==
                std::tie(row.feature_id, row.transition_id, row.hypothesis, row.bmask, row.num_peptidoforms, row.alignment_group_id))
            {
              back.pep = std::min(back.pep, row.pep);
              back.precursor_peakgroup_pep = std::min(back.precursor_peakgroup_pep, row.precursor_peakgroup_pep);
              continue;
            }
          }
          result.push_back(row);
        }
      }

      return result;
    }

    std::vector<OpenSwathPeptidoformInference::BayesianModelRow> prepareTransitionBM_(const std::vector<TransitionInferenceRow>& rows)
    {
      std::vector<OpenSwathPeptidoformInference::BayesianModelRow> bm_rows;
      bm_rows.reserve(rows.size());
      for (const auto& row : rows)
      {
        if (row.hypothesis != -1 && row.num_peptidoforms <= 0)
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Transition inference requires num_peptidoforms > 0 for non-null hypotheses.");
        }

        OpenSwathPeptidoformInference::BayesianModelRow bm_row;
        bm_row.feature_id = row.feature_id;
        bm_row.hypothesis = row.hypothesis;
        if (row.hypothesis == -1)
        {
          // H0/null gets the precursor-level error mass as its prior.
          bm_row.prior = row.precursor_peakgroup_pep;
        }
        else
        {
          // Non-null peptidoform hypotheses share the remaining precursor-level
          // confidence equally, matching the IPF formulation in PyProphet.
          bm_row.prior = (1.0 - row.precursor_peakgroup_pep) / static_cast<double>(row.num_peptidoforms);
        }
        bm_row.evidence = (row.bmask == 1) ? (1.0 - row.pep) : row.pep;
        bm_rows.push_back(bm_row);
      }
      return bm_rows;
    }

    std::unordered_map<Int64, double> makePrecursorPepMap_(const std::vector<IPFPrecursorProbabilityRow>& rows)
    {
      std::unordered_map<Int64, double> pep_map;
      pep_map.reserve(rows.size());
      for (const auto& row : rows)
      {
        pep_map[row.feature_id] = row.precursor_peakgroup_pep;
      }
      return pep_map;
    }
  } // namespace

  std::vector<double> OpenSwathPeptidoformInference::computeModelFDR(const std::vector<double>& pep_values)
  {
    return Math::MultipleTesting::computeModelFDR(pep_values);
  }

  std::vector<OpenSwathPeptidoformInference::BayesianModelRow> OpenSwathPeptidoformInference::preparePrecursorBM(const std::vector<IPFPrecursorRow>& rows)
  {
    std::vector<BayesianModelRow> bm_rows;
    bm_rows.reserve(rows.size() * 4);
    std::unordered_set<Int64> features_with_precursor_evidence;
    features_with_precursor_evidence.reserve(rows.size());

    for (const auto& row : rows)
    {
      if (row.ms1_precursor_pep.has_value())
      {
        features_with_precursor_evidence.insert(row.feature_id);
        // hypothesis 1 = true precursor state, hypothesis 0 = false precursor state
        bm_rows.push_back({row.feature_id, 1, 1.0 - row.ms2_peakgroup_pep, 1.0 - *row.ms1_precursor_pep});
        bm_rows.push_back({row.feature_id, 0, row.ms2_peakgroup_pep, *row.ms1_precursor_pep});
      }
      if (row.ms2_precursor_pep.has_value())
      {
        features_with_precursor_evidence.insert(row.feature_id);
        // hypothesis 1 = true precursor state, hypothesis 0 = false precursor state
        bm_rows.push_back({row.feature_id, 1, 1.0 - row.ms2_peakgroup_pep, 1.0 - *row.ms2_precursor_pep});
        bm_rows.push_back({row.feature_id, 0, row.ms2_peakgroup_pep, *row.ms2_precursor_pep});
      }
    }

    std::unordered_set<Int64> missing_rows_added;
    missing_rows_added.reserve(rows.size());
    for (const auto& row : rows)
    {
      if (features_with_precursor_evidence.find(row.feature_id) == features_with_precursor_evidence.end() &&
          missing_rows_added.insert(row.feature_id).second)
      {
        // No precursor-specific evidence is available for this feature, so fall back
        // to the same synthetic true/false rows used by PyProphet.
        bm_rows.push_back({row.feature_id, 1, 1.0 - row.ms2_peakgroup_pep, 0.0});
        bm_rows.push_back({row.feature_id, 0, row.ms2_peakgroup_pep, 1.0});
      }
    }

    return bm_rows;
  }

  std::vector<OpenSwathPeptidoformInference::PosteriorRow> OpenSwathPeptidoformInference::applyBM(const std::vector<BayesianModelRow>& rows)
  {
    if (rows.empty())
    {
      return {};
    }

    std::vector<BayesianModelRow> sorted_rows = rows;
    std::sort(sorted_rows.begin(), sorted_rows.end(),
              [](const BayesianModelRow& lhs, const BayesianModelRow& rhs)
              {
                return std::tie(lhs.feature_id, lhs.hypothesis) < std::tie(rhs.feature_id, rhs.hypothesis);
              });

    struct LikelihoodPriorRow
    {
      Int64 feature_id = -1;
      Int64 hypothesis = -1;
      double likelihood_prior = 0.0;
    };

    std::vector<LikelihoodPriorRow> likelihood_rows;
    likelihood_rows.reserve(sorted_rows.size());
    Size i = 0;
    while (i < sorted_rows.size())
    {
      const Int64 feature_id = sorted_rows[i].feature_id;
      const Int64 hypothesis = sorted_rows[i].hypothesis;
      double evidence_product = 1.0;
      double prior = sorted_rows[i].prior;
      Size j = i;
      // For each feature/hypothesis pair, multiply all evidence terms together and
      // keep the smallest prior seen in the group, which matches the compact BM
      // reduction done in PyProphet.
      while (j < sorted_rows.size() &&
             sorted_rows[j].feature_id == feature_id &&
             sorted_rows[j].hypothesis == hypothesis)
      {
        evidence_product *= sorted_rows[j].evidence;
        prior = std::min(prior, sorted_rows[j].prior);
        ++j;
      }
      likelihood_rows.push_back({feature_id, hypothesis, evidence_product * prior});
      i = j;
    }

    std::vector<PosteriorRow> posterior_rows;
    posterior_rows.reserve(likelihood_rows.size());
    i = 0;
    while (i < likelihood_rows.size())
    {
      const Int64 feature_id = likelihood_rows[i].feature_id;
      double likelihood_sum = 0.0;
      Size j = i;
      while (j < likelihood_rows.size() && likelihood_rows[j].feature_id == feature_id)
      {
        likelihood_sum += likelihood_rows[j].likelihood_prior;
        ++j;
      }

      // Degenerate groups can end up with a zero or non-finite denominator. In that
      // case we keep a deterministic posterior of zero for every hypothesis instead
      // of crashing.
      const bool valid_sum = std::isfinite(likelihood_sum) && likelihood_sum > 0.0;
      for (Size k = i; k < j; ++k)
      {
        const double posterior = valid_sum ? (likelihood_rows[k].likelihood_prior / likelihood_sum) : 0.0;
        posterior_rows.push_back({feature_id, likelihood_rows[k].hypothesis, posterior});
      }
      i = j;
    }

    return posterior_rows;
  }

  std::vector<IPFPrecursorProbabilityRow> OpenSwathPeptidoformInference::precursorInference(const std::vector<IPFPrecursorRow>& precursor_rows,
                                                                                             const PeptidoformInferenceConfig& config) const
  {
    std::vector<IPFPrecursorProbabilityRow> inferred_precursors;
    inferred_precursors.reserve(precursor_rows.size());

    if (!(config.ipf_ms1_scoring || config.ipf_ms2_scoring))
    {
      for (const auto& row : precursor_rows)
      {
        if (row.ms2_peakgroup_pep < config.ipf_max_precursor_peakgroup_pep)
        {
          inferred_precursors.push_back({row.feature_id, row.ms2_peakgroup_pep});
        }
      }
      return inferred_precursors;
    }

    std::vector<IPFPrecursorRow> filtered_rows;
    filtered_rows.reserve(precursor_rows.size());
    for (const auto& row : precursor_rows)
    {
      IPFPrecursorRow filtered = row;
      // Precursor evidence is thresholded independently per evidence source before
      // entering the Bayesian model. Sources that fail the threshold are treated as
      // absent, exactly like in the PyProphet implementation.
      if (!config.ipf_ms1_scoring || !filtered.ms1_precursor_pep.has_value() || *filtered.ms1_precursor_pep >= config.ipf_max_precursor_pep)
      {
        filtered.ms1_precursor_pep.reset();
      }
      if (!config.ipf_ms2_scoring || !filtered.ms2_precursor_pep.has_value() || *filtered.ms2_precursor_pep >= config.ipf_max_precursor_pep)
      {
        filtered.ms2_precursor_pep.reset();
      }
      filtered_rows.push_back(std::move(filtered));
    }

    const std::vector<BayesianModelRow> bm_rows = preparePrecursorBM(filtered_rows);
    const std::vector<PosteriorRow> posterior_rows = applyBM(bm_rows);
    inferred_precursors.reserve(posterior_rows.size());
    for (const auto& row : posterior_rows)
    {
      if (row.hypothesis == 1)
      {
        const double precursor_peakgroup_pep = 1.0 - row.posterior;
        if (precursor_peakgroup_pep < config.ipf_max_precursor_peakgroup_pep)
        {
          inferred_precursors.push_back({row.feature_id, precursor_peakgroup_pep});
        }
      }
    }
    return inferred_precursors;
  }

  std::vector<IPFResultRow> OpenSwathPeptidoformInference::infer(const std::vector<IPFPrecursorRow>& precursor_rows,
                                                                 const std::vector<IPFTransitionRow>& transition_rows,
                                                                 const std::vector<IPFAlignmentRow>& alignment_rows,
                                                                 const PeptidoformInferenceConfig& config) const
  {
    const std::vector<IPFPrecursorProbabilityRow> precursor_probabilities = precursorInference(precursor_rows, config);
    if (precursor_probabilities.empty() || transition_rows.empty())
    {
      return {};
    }

    const std::unordered_map<Int64, double> precursor_pep_map = makePrecursorPepMap_(precursor_probabilities);
    std::unordered_map<Int64, Int64> alignment_group_map;
    alignment_group_map.reserve(alignment_rows.size());
    for (const auto& row : alignment_rows)
    {
      alignment_group_map[row.feature_id] = row.alignment_group_id;
    }

    // Transition-heavy jobs can dominate memory. Convert immediately to a narrow internal row
    // that keeps only the fields required for the Bayesian model and release intermediate maps
    // as soon as possible.
    std::vector<TransitionInferenceRow> merged_rows;
    merged_rows.reserve(transition_rows.size());
    std::unordered_map<Int64, Int32> num_peptidoforms_by_feature;
    num_peptidoforms_by_feature.reserve(transition_rows.size());
    for (const auto& row : transition_rows)
    {
      const auto precursor_it = precursor_pep_map.find(row.feature_id);
      if (precursor_it == precursor_pep_map.end())
      {
        continue;
      }

      TransitionInferenceRow merged_row;
      merged_row.feature_id = row.feature_id;
      merged_row.transition_id = row.transition_id;
      merged_row.hypothesis = row.peptide_id;
      merged_row.pep = row.pep;
      merged_row.precursor_peakgroup_pep = precursor_it->second;
      merged_row.bmask = row.bmask;
      merged_row.num_peptidoforms = row.num_peptidoforms;
      merged_row.alignment_group_id = row.alignment_group_id.has_value() ?
        *row.alignment_group_id :
        (alignment_group_map.count(row.feature_id) ? alignment_group_map[row.feature_id] : -1);
      merged_rows.push_back(std::move(merged_row));

      if (row.num_peptidoforms > 0 && num_peptidoforms_by_feature.find(row.feature_id) == num_peptidoforms_by_feature.end())
      {
        num_peptidoforms_by_feature[row.feature_id] = row.num_peptidoforms;
      }
    }

    if (merged_rows.empty())
    {
      return {};
    }

    if (config.propagate_signal_across_runs)
    {
      merged_rows = propagateSignalAcrossRuns_(merged_rows, config.across_run_confidence_threshold);
    }

    // At this point the heavy transition evidence has been reduced to a compact BM
    // input. The later steps only carry posterior/q-value sized vectors.
    const std::vector<BayesianModelRow> bm_rows = prepareTransitionBM_(merged_rows);
    const std::vector<PosteriorRow> posterior_rows = applyBM(bm_rows);
    if (posterior_rows.empty())
    {
      return {};
    }

    std::vector<double> pep_values;
    pep_values.reserve(posterior_rows.size());
    for (const auto& row : posterior_rows)
    {
      pep_values.push_back(1.0 - row.posterior);
    }

    std::vector<double> qvalues(posterior_rows.size(), 1.0);
    if (config.ipf_grouped_fdr)
    {
      std::unordered_map<Int32, std::vector<Size>> grouped_indices;
      grouped_indices.reserve(num_peptidoforms_by_feature.size());
      for (Size idx = 0; idx < posterior_rows.size(); ++idx)
      {
        const auto num_it = num_peptidoforms_by_feature.find(posterior_rows[idx].feature_id);
        const Int32 num_peptidoforms = (num_it == num_peptidoforms_by_feature.end()) ? 0 : num_it->second;
        grouped_indices[num_peptidoforms].push_back(idx);
      }
      for (const auto& group : grouped_indices)
      {
        std::vector<double> grouped_peps;
        grouped_peps.reserve(group.second.size());
        for (Size idx : group.second)
        {
          grouped_peps.push_back(pep_values[idx]);
        }
        const std::vector<double> grouped_qvalues = computeModelFDR(grouped_peps);
        for (Size i = 0; i < group.second.size(); ++i)
        {
          qvalues[group.second[i]] = grouped_qvalues[i];
        }
      }
    }
    else
    {
      qvalues = computeModelFDR(pep_values);
    }

    std::vector<IPFResultRow> results;
    results.reserve(posterior_rows.size());
    for (Size i_row = 0; i_row < posterior_rows.size(); ++i_row)
    {
      const PosteriorRow& row = posterior_rows[i_row];
      if (row.hypothesis == -1)
      {
        continue;
      }
      const auto precursor_it = precursor_pep_map.find(row.feature_id);
      if (precursor_it == precursor_pep_map.end())
      {
        continue;
      }
      results.push_back({row.feature_id, row.hypothesis, precursor_it->second, qvalues[i_row], pep_values[i_row]});
    }
    std::sort(results.begin(), results.end(),
              [](const IPFResultRow& lhs, const IPFResultRow& rhs)
              {
                return std::tie(lhs.feature_id, lhs.peptide_id) < std::tie(rhs.feature_id, rhs.peptide_id);
              });
    return results;
  }
} // namespace OpenMS
