// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIonLikelihoodModel.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Macros.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace OpenMS
{
  FragmentIonLikelihoodModel::FragmentIonLikelihoodModel() :
    FragmentIonLikelihoodModel(20.0)
  {
  }

  FragmentIonLikelihoodModel::FragmentIonLikelihoodModel(double pseudo_count) :
    signal_counts_(CONTEXTS * OUTCOMES, 0.0),
    noise_counts_(CONTEXTS * OUTCOMES, 0.0),
    pseudo_count_(pseudo_count)
  {
    if (! std::isfinite(pseudo_count) || pseudo_count <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The pseudo-count must be finite and positive.");
    }
  }

  std::vector<Size> FragmentIonLikelihoodModel::intensityRanks(const MSSpectrum& spectrum)
  {
    std::vector<Size> order(spectrum.size());
    std::iota(order.begin(), order.end(), Size(0));
    // Ties are broken by peak order (ascending m/z for a sorted spectrum) to keep ranks deterministic.
    std::stable_sort(order.begin(), order.end(), [&spectrum](Size a, Size b) { return spectrum[a].getIntensity() > spectrum[b].getIntensity(); });
    std::vector<Size> ranks(spectrum.size());
    for (Size rank = 0; rank < order.size(); ++rank)
    {
      ranks[order[rank]] = rank + 1;
    }
    return ranks;
  }

  Size FragmentIonLikelihoodModel::rankOutcome(Size rank)
  {
    // Bins double roughly with the rank: the top peaks carry most of the evidence, the tail little.
    if (rank <= 2) return 0;
    if (rank <= 5) return 1;
    if (rank <= 10) return 2;
    if (rank <= 20) return 3;
    if (rank <= 40) return 4;
    if (rank <= 80) return 5;
    return 6;
  }

  FragmentIonLikelihoodModel::Context FragmentIonLikelihoodModel::contextOf(bool prefix, int precursor_charge, int fragment_charge,
                                                                            Size fragment_length, Size peptide_length)
  {
    Context context;
    context.series = prefix ? 0 : 1;
    context.precursor_bucket = precursor_charge <= 2 ? 0 : (precursor_charge == 3 ? 1 : 2);
    context.fragment_charge = fragment_charge <= 1 ? 0 : 1;
    const Size bin = peptide_length == 0 ? 0 : (POSITION_BINS * fragment_length) / peptide_length;
    context.position_bin = std::min(bin, POSITION_BINS - 1);
    return context;
  }

  bool FragmentIonLikelihoodModel::parseIonName(const std::string& name, bool& prefix, Size& ordinal)
  {
    if (name.empty()) return false;
    // Cross-link annotations put the ion type after a '$'; plain annotations start with it.
    Size pos = 0;
    if (const Size dollar = name.find('$'); dollar != std::string::npos && dollar + 1 < name.size())
    {
      pos = dollar + 1;
    }
    const char series = name[pos];
    if (series == 'a' || series == 'b' || series == 'c') { prefix = true; }
    else if (series == 'x' || series == 'y' || series == 'z') { prefix = false; }
    else { return false; }
    ++pos;
    // z. (z+1) and z' (z+2) radical ions carry a marker between letter and ordinal.
    while (pos < name.size() && (name[pos] == '.' || name[pos] == '\'')) ++pos;
    const Size begin = pos;
    while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9') ++pos;
    if (pos == begin) return false;
    ordinal = static_cast<Size>(std::stoul(name.substr(begin, pos - begin)));
    return ordinal > 0;
  }

  Size FragmentIonLikelihoodModel::index_(const Context& context, Size outcome)
  {
    const Size cell = ((context.series * PRECURSOR_BUCKETS + context.precursor_bucket) * FRAGMENT_CHARGES + context.fragment_charge) * POSITION_BINS
                      + context.position_bin;
    return cell * OUTCOMES + outcome;
  }

  void FragmentIonLikelihoodModel::matchIons_(const MSSpectrum& spectrum,
                                              const std::vector<Size>& ranks,
                                              const MSSpectrum& theoretical,
                                              Size peptide_length,
                                              int precursor_charge,
                                              double tolerance,
                                              bool ppm,
                                              std::vector<Ion_>& ions)
  {
    if (! std::isfinite(tolerance) || tolerance <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Matching tolerance must be finite and positive.");
    }
    if (theoretical.getStringDataArrays().empty() || theoretical.getStringDataArrays()[0].size() != theoretical.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Theoretical spectrum needs one ion name per peak", "IonNames");
    }
    if (ranks.size() != spectrum.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "One intensity rank per experimental peak required", "ranks");
    }
    const auto& names = theoretical.getStringDataArrays()[0];
    // TheoreticalSpectrumGenerator writes the fragment charges as the first integer array; without it every ion counts as singly charged.
    const auto* charges = (! theoretical.getIntegerDataArrays().empty() && theoretical.getIntegerDataArrays()[0].size() == theoretical.size())
                            ? &theoretical.getIntegerDataArrays()[0]
                            : nullptr;
    ions.clear();
    ions.reserve(theoretical.size());
    for (Size i = 0; i < theoretical.size(); ++i)
    {
      bool prefix = false;
      Size ordinal = 0;
      if (! parseIonName(names[i], prefix, ordinal)) continue;
      const int fragment_charge = charges != nullptr ? (*charges)[i] : 1;
      const double mz = theoretical[i].getMZ();
      const double width = ppm ? mz * tolerance * 1e-6 : tolerance;
      Size outcome = ABSENT;
      if (! spectrum.empty())
      {
        // Closest peak within the window, as the scorers match.
        const Int nearest = spectrum.findNearest(mz, width);
        if (nearest >= 0) outcome = rankOutcome(ranks[static_cast<Size>(nearest)]);
      }
      ions.push_back({contextOf(prefix, precursor_charge, fragment_charge, ordinal, peptide_length), outcome});
    }
  }

  void FragmentIonLikelihoodModel::addObservations(const MSSpectrum& spectrum,
                                                   const std::vector<Size>& ranks,
                                                   const MSSpectrum& theoretical,
                                                   Size peptide_length,
                                                   int precursor_charge,
                                                   double tolerance,
                                                   bool ppm,
                                                   bool noise)
  {
    std::vector<Ion_> ions;
    matchIons_(spectrum, ranks, theoretical, peptide_length, precursor_charge, tolerance, ppm, ions);
    std::vector<double>& counts = noise ? noise_counts_ : signal_counts_;
    for (const Ion_& ion : ions)
    {
      counts[index_(ion.context, ion.outcome)] += 1.0;
    }
    if (noise) { ++noise_psms_; }
    else { ++signal_psms_; }
    finalized_ = false;
  }

  std::vector<double> FragmentIonLikelihoodModel::smooth_(const std::vector<double>& counts, double pseudo_count)
  {
    // Three back-off levels above the full context: (series, precursor bucket, fragment charge),
    // (series) and the global outcome distribution with a flat Laplace prior. Each level's
    // probability blends its own counts with its parent's probability, weighted by the pseudo-count,
    // so sparse cells inherit from their parents and unseen cells stay finite.
    std::vector<double> global(OUTCOMES, 0.0);
    for (Size cell = 0; cell < CONTEXTS; ++cell)
    {
      for (Size o = 0; o < OUTCOMES; ++o) global[o] += counts[cell * OUTCOMES + o];
    }
    const double global_total = std::accumulate(global.begin(), global.end(), 0.0);
    std::vector<double> global_p(OUTCOMES);
    for (Size o = 0; o < OUTCOMES; ++o) global_p[o] = (global[o] + 1.0) / (global_total + static_cast<double>(OUTCOMES));

    auto blend = [pseudo_count](const std::vector<double>& own, const std::vector<double>& parent) {
      const double total = std::accumulate(own.begin(), own.end(), 0.0);
      std::vector<double> p(OUTCOMES);
      for (Size o = 0; o < OUTCOMES; ++o) p[o] = (own[o] + pseudo_count * parent[o]) / (total + pseudo_count);
      return p;
    };

    std::vector<double> logp(CONTEXTS * OUTCOMES, 0.0);
    for (Size series = 0; series < SERIES; ++series)
    {
      std::vector<double> series_counts(OUTCOMES, 0.0);
      for (Size bucket = 0; bucket < PRECURSOR_BUCKETS; ++bucket)
      {
        for (Size charge = 0; charge < FRAGMENT_CHARGES; ++charge)
        {
          for (Size position = 0; position < POSITION_BINS; ++position)
          {
            const Size cell = index_({series, bucket, charge, position}, 0);
            for (Size o = 0; o < OUTCOMES; ++o) series_counts[o] += counts[cell + o];
          }
        }
      }
      const std::vector<double> series_p = blend(series_counts, global_p);
      for (Size bucket = 0; bucket < PRECURSOR_BUCKETS; ++bucket)
      {
        for (Size charge = 0; charge < FRAGMENT_CHARGES; ++charge)
        {
          std::vector<double> charge_counts(OUTCOMES, 0.0);
          for (Size position = 0; position < POSITION_BINS; ++position)
          {
            const Size cell = index_({series, bucket, charge, position}, 0);
            for (Size o = 0; o < OUTCOMES; ++o) charge_counts[o] += counts[cell + o];
          }
          const std::vector<double> charge_p = blend(charge_counts, series_p);
          for (Size position = 0; position < POSITION_BINS; ++position)
          {
            const Size cell = index_({series, bucket, charge, position}, 0);
            const std::vector<double> own(counts.begin() + cell, counts.begin() + cell + OUTCOMES);
            const std::vector<double> p = blend(own, charge_p);
            for (Size o = 0; o < OUTCOMES; ++o) logp[cell + o] = std::log(p[o]);
          }
        }
      }
    }
    return logp;
  }

  void FragmentIonLikelihoodModel::finalize()
  {
    signal_logp_ = smooth_(signal_counts_, pseudo_count_);
    noise_logp_ = smooth_(noise_counts_, pseudo_count_);
    finalized_ = true;
  }

  void FragmentIonLikelihoodModel::requireTrained_() const
  {
    if (! finalized_)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "finalize() must be called after the last observation and before scoring.");
    }
  }

  void FragmentIonLikelihoodModel::requireValidContext_(const Context& context)
  {
    // Contexts come from callers (also through pyOpenMS); every dimension must lie inside the tables.
    if (context.series >= SERIES || context.precursor_bucket >= PRECURSOR_BUCKETS || context.fragment_charge >= FRAGMENT_CHARGES
        || context.position_bin >= POSITION_BINS)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Context dimension out of range (see contextOf())", "Context");
    }
  }

  double FragmentIonLikelihoodModel::logLikelihoodRatio(const Context& context, Size outcome) const
  {
    requireTrained_();
    requireValidContext_(context);
    const Size i = index_(context, std::min(outcome, ABSENT));
    return signal_logp_[i] - noise_logp_[i];
  }

  double FragmentIonLikelihoodModel::presenceProbability(const Context& context) const
  {
    requireTrained_();
    requireValidContext_(context);
    return 1.0 - std::exp(signal_logp_[index_(context, ABSENT)]);
  }

  FragmentIonLikelihoodModel::Features FragmentIonLikelihoodModel::score(const MSSpectrum& spectrum,
                                                                         const std::vector<Size>& ranks,
                                                                         const MSSpectrum& theoretical,
                                                                         Size peptide_length,
                                                                         int precursor_charge,
                                                                         double tolerance,
                                                                         bool ppm,
                                                                         Size top_k) const
  {
    requireTrained_();
    std::vector<Ion_> ions;
    matchIons_(spectrum, ranks, theoretical, peptide_length, precursor_charge, tolerance, ppm, ions);
    Features features;
    features.theoretical_ions = ions.size();
    if (ions.empty()) return features;

    double presence_total = 0.0, presence_matched = 0.0;
    std::vector<std::pair<double, bool>> predicted; // (P(present), matched) per ion
    predicted.reserve(ions.size());
    for (const Ion_& ion : ions)
    {
      const bool matched = ion.outcome != ABSENT;
      features.log_likelihood_ratio += logLikelihoodRatio(ion.context, ion.outcome);
      const double presence = presenceProbability(ion.context);
      presence_total += presence;
      if (matched)
      {
        presence_matched += presence;
        ++features.matched_ions;
      }
      predicted.emplace_back(presence, matched);
    }
    features.explained_presence = presence_total > 0.0 ? presence_matched / presence_total : 0.0;

    // Of the ions the run predicts most confidently, how many did this PSM deliver? Ties in the
    // presence probability are resolved by ion order to keep the feature deterministic.
    std::stable_sort(predicted.begin(), predicted.end(), [](const auto& a, const auto& b) { return a.first > b.first; });
    const Size k = std::min(top_k, predicted.size());
    if (k > 0)
    {
      const Size observed = static_cast<Size>(std::count_if(predicted.begin(), predicted.begin() + k, [](const auto& p) { return p.second; }));
      features.top_predicted_observed = static_cast<double>(observed) / static_cast<double>(k);
    }
    return features;
  }
} // namespace OpenMS
