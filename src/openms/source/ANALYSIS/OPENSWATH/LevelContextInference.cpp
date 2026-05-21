// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/MATH/STATISTICS/MultipleTesting.h>

#include <algorithm>
#include <cmath>
#include <limits>
#include <map>
#include <numeric>
#include <tuple>
#include <vector>

namespace OpenMS
{
  namespace
  {
    using ErrorTableRow = std::tuple<double, double, double, double>;

    struct Metrics
    {
      std::vector<double> svalues;
    };

    std::vector<Int64> countNumPositives_(const std::vector<double>& values)
    {
      std::vector<Int64> counts(values.size(), 0);
      if (values.empty())
      {
        return counts;
      }

      // The input is a sorted list of p-values/cutoffs. For each cutoff we store how
      // many observations remain positive at or above that cutoff. Tied cutoffs share
      // the same "maximum-rank" count, matching the tie handling used in PyProphet.
      Size tie_group_begin = 0;
      Size tie_group_end = 0;
      Int64 remaining_positives = static_cast<Int64>(values.size());
      while (tie_group_begin < values.size())
      {
        while (tie_group_end < values.size() && values[tie_group_begin] == values[tie_group_end])
        {
          counts[tie_group_end] = remaining_positives;
          ++tie_group_end;
        }
        --remaining_positives;
        ++tie_group_begin;
      }
      return counts;
    }

    Size findNearestMatchAscending_(const std::vector<double>& basis, const double sample)
    {
      if (basis.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Nearest-match lookup requires a non-empty basis vector.");
      }

      if (basis.size() == 1)
      {
        return 0;
      }

      Size low = 0;
      Size high = basis.size() - 1;
      Size best = std::numeric_limits<Size>::max();

      if (basis[low] == sample)
      {
        best = low;
      }
      else if (basis[high] == sample)
      {
        best = high;
      }
      else
      {
        while (low < high - 1)
        {
          const Size mid = (low + high) / 2;
          if (basis[mid] == sample)
          {
            best = mid;
          }
          if (basis[mid] < sample)
          {
            low = mid;
          }
          else
          {
            high = mid;
          }
        }
        if (best == std::numeric_limits<Size>::max())
        {
          const double low_dist = std::fabs(basis[low] - sample);
          const double high_dist = std::fabs(basis[high] - sample);
          best = (low_dist < high_dist) ? low : high;
        }
      }

      while (best > 0 && basis[best - 1] == basis[best])
      {
        --best;
      }
      return best;
    }

    Metrics statMetrics_(const std::vector<double>& p_values, const double pi0, const bool pfdr)
    {
      Metrics metrics;
      const Size num_total = p_values.size();
      metrics.svalues.assign(num_total, 0.0);
      if (num_total == 0)
      {
        return metrics;
      }

      const std::vector<Int64> num_positives = countNumPositives_(p_values);
      const double num_total_d = static_cast<double>(num_total);
      const double num_null = pi0 * num_total_d;
      const double sens_denom = std::max(0.0, num_total_d - num_null);
      std::vector<double> sens(num_total, 0.0);
      for (Size i = 0; i < num_total; ++i)
      {
        // num_positives[i] is the number of observations called positive at the
        // current cutoff. Under the two-group model, pi0 controls the expected null
        // fraction, which lets us decompose positives/negatives into expected
        // TP/FP/TN/FN counts mirroring PyProphet.
        double tp = static_cast<double>(num_positives[i]) - num_null * p_values[i];
        double fp = num_null * p_values[i];
        double tn = num_null * (1.0 - p_values[i]);
        double fn = (num_total_d - static_cast<double>(num_positives[i])) - num_null * (1.0 - p_values[i]);
        double fdr = 0.0;
        double fnr = 0.0;
        if (num_positives[i] != 0)
        {
          fdr = fp / static_cast<double>(num_positives[i]);
        }
        const double num_negatives = num_total_d - static_cast<double>(num_positives[i]);
        if (num_negatives != 0.0)
        {
          fnr = fn / num_negatives;
        }
        if (pfdr)
        {
          // Optional pFDR scaling mirrors PyProphet's positive-FDR correction for
          // both FDR and FNR style quantities.
          const double p = p_values[i];
          const double fdr_scale = 1.0 - std::pow(1.0 - p, num_total_d);
          if (p == 0.0)
          {
            fdr = 1.0 / num_total_d;
          }
          else if (fdr_scale != 0.0)
          {
            fdr /= fdr_scale;
          }

          const double fnr_scale = 1.0 - std::pow(p, num_total_d);
          if (p == 0.0)
          {
            fnr = 1.0 / num_total_d;
          }
          else if (fnr_scale != 0.0)
          {
            fnr /= fnr_scale;
          }
        }

        if (sens_denom != 0.0)
        {
          sens[i] = tp / sens_denom;
        }
        sens[i] = std::clamp(sens[i], 0.0, 1.0);
        fdr = std::clamp(fdr, 0.0, 1.0);
        fnr = std::clamp(fnr, 0.0, 1.0);
        if (num_positives[i] == 0)
        {
          fdr = 0.0;
          fnr = 0.0;
        }
        (void)tn;
      }

      double running_max = 0.0;
      for (SignedSize i = static_cast<SignedSize>(sens.size()) - 1; i >= 0; --i)
      {
        // s-values are the monotone running maximum of sensitivity when moving from
        // stringent to permissive cutoffs.
        running_max = std::max(running_max, sens[static_cast<Size>(i)]);
        metrics.svalues[static_cast<Size>(i)] = running_max;
      }
      return metrics;
    }

    struct ErrorStatistics
    {
      std::vector<double> cutoffs;
      std::vector<double> pvalues;
      std::vector<double> qvalues;
      std::vector<double> peps;
    };

    ErrorStatistics buildErrorStatistics_(const std::vector<double>& target_scores,
                                          const std::vector<double>& decoy_scores,
                                          const ErrorEstimationConfig& config)
    {
      if (target_scores.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Target score distribution is empty.");
      }
      if (decoy_scores.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decoy score distribution is empty.");
      }

      std::vector<double> sorted_targets;
      sorted_targets.reserve(target_scores.size());
      for (double score : target_scores)
      {
        if (std::isfinite(score))
        {
          sorted_targets.push_back(score);
        }
      }
      std::sort(sorted_targets.begin(), sorted_targets.end());

      std::vector<double> sorted_decoys;
      sorted_decoys.reserve(decoy_scores.size());
      for (double score : decoy_scores)
      {
        if (std::isfinite(score))
        {
          sorted_decoys.push_back(score);
        }
      }
      std::sort(sorted_decoys.begin(), sorted_decoys.end());

      if (sorted_targets.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Target score distribution only contained NaN values.");
      }
      if (sorted_decoys.empty())
      {
        throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Decoy score distribution only contained NaN values.");
      }

      // OpenMS MultipleTesting expects compact target and decoy score vectors. The
      // target scores become the cutoff basis we later map result rows back onto.
      std::vector<double> pvalues = config.parametric ?
        Math::MultipleTesting::pNorm(sorted_targets, sorted_decoys) :
        Math::MultipleTesting::pEmp(sorted_targets, sorted_decoys);
      const Math::Pi0Result pi0 = Math::MultipleTesting::pi0Est(
        pvalues,
        config.pi0_lambda,
        Math::MultipleTesting::toPi0Method(config.pi0_method),
        config.pi0_smooth_df,
        config.pi0_smooth_log_pi0
      );
      std::vector<double> qvalues = Math::MultipleTesting::qValue(pvalues, pi0.pi0, config.pfdr);
      std::vector<double> peps = Math::MultipleTesting::lfdr(
        pvalues,
        pi0.pi0,
        config.lfdr_truncate,
        config.lfdr_monotone,
        Math::MultipleTesting::toLfdrTransform(config.lfdr_transformation),
        config.lfdr_adj,
        config.lfdr_eps
      );

      // Keep the PyProphet metric construction order even though only p/q/pep are returned.
      (void)statMetrics_(pvalues, pi0.pi0, config.pfdr);

      return {std::move(sorted_targets), std::move(pvalues), std::move(qvalues), std::move(peps)};
    }

    template <typename IndexedRowRange>
    void assignErrorStatistics_(const IndexedRowRange& rows,
                                const ErrorEstimationConfig& config,
                                std::vector<LevelContextResultRow>& results)
    {
      std::vector<double> target_scores;
      std::vector<double> decoy_scores;
      target_scores.reserve(rows.size());
      decoy_scores.reserve(rows.size());
      for (const auto& indexed_row : rows)
      {
        if (indexed_row.second.decoy)
        {
          decoy_scores.push_back(indexed_row.second.score);
        }
        else
        {
          target_scores.push_back(indexed_row.second.score);
        }
      }

      const ErrorStatistics error_stats = buildErrorStatistics_(target_scores, decoy_scores, config);
      for (const auto& indexed_row : rows)
      {
        const Size result_index = indexed_row.first;
        const LevelContextInputRow& row = indexed_row.second;
        LevelContextResultRow& result = results[result_index];
        result.run_id = row.run_id;
        result.entity_id = row.entity_id;
        result.context = row.context;
        result.score = row.score;

        if (!std::isfinite(row.score))
        {
          result.pvalue = std::numeric_limits<double>::quiet_NaN();
          result.qvalue = std::numeric_limits<double>::quiet_NaN();
          result.pep = std::numeric_limits<double>::quiet_NaN();
          continue;
        }

        // Error statistics are computed on the sorted target-score cutoffs. Each row,
        // whether target or decoy, is projected onto the nearest available cutoff to
        // retrieve the corresponding p-value, q-value, and PEP estimate.
        const Size match_index = findNearestMatchAscending_(error_stats.cutoffs, row.score);
        result.pvalue = error_stats.pvalues[match_index];
        result.qvalue = error_stats.qvalues[match_index];
        result.pep = error_stats.peps[match_index];
      }
    }
  } // namespace

  std::vector<LevelContextResultRow> LevelContextInference::infer(const std::vector<LevelContextInputRow>& input,
                                                                  const LevelContextInferenceConfig& config)
  {
    if (config.level == InferenceLevel::Peptidoform)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Peptidoform inference is not supported by LevelContextInference.");
    }
    if (input.empty())
    {
      return {};
    }

    std::vector<LevelContextResultRow> results(input.size());
    std::vector<std::pair<Size, LevelContextInputRow>> indexed_rows;
    indexed_rows.reserve(input.size());
    for (Size i = 0; i < input.size(); ++i)
    {
      indexed_rows.emplace_back(i, input[i]);
    }

    if (config.context == InferenceContext::RunSpecific)
    {
      std::map<Int64, std::vector<std::pair<Size, LevelContextInputRow>>> grouped_rows;
      for (const auto& indexed_row : indexed_rows)
      {
        if (!indexed_row.second.run_id.has_value())
        {
          throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Run-specific inference requires RUN_ID values for all rows.");
        }
        grouped_rows[*indexed_row.second.run_id].push_back(indexed_row);
      }

      for (const auto& run_group : grouped_rows)
      {
        try
        {
          assignErrorStatistics_(run_group.second, config.error, results);
        }
        catch (const Exception::BaseException& e)
        {
          throw Exception::Precondition(
            __FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Run-specific inference failed for RUN_ID " + String(run_group.first) + ": " + e.getMessage()
          );
        }
      }
    }
    else
    {
      assignErrorStatistics_(indexed_rows, config.error, results);
    }
    return results;
  }
} // namespace OpenMS
