// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Oliver Kohlbacher $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FragmentIonLikelihoodModel.h>

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/Exception.h>

#include <algorithm>
#include <cmath>
#include <numeric>

namespace OpenMS
{
  // ---------------------------------------------------------------------------
  // PeakLists
  // ---------------------------------------------------------------------------

  void FragmentIonLikelihoodModel::PeakLists::reset(Size spectra)
  {
    clear();
    mz_.resize(spectra);
    rank_bins_.resize(spectra);
  }

  void FragmentIonLikelihoodModel::PeakLists::clear()
  {
    std::vector<std::vector<float>>().swap(mz_);
    std::vector<std::vector<std::uint8_t>>().swap(rank_bins_);
  }

  void FragmentIonLikelihoodModel::PeakLists::assign(Size index, const MSSpectrum& spectrum)
  {
    if (index >= mz_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, index, mz_.size());
    }
    if (! spectrum.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The spectrum must be sorted by m/z.");
    }
    const Size n = spectrum.size();
    std::vector<float> mz(n);
    std::vector<std::uint8_t> bins(n, static_cast<std::uint8_t>(RANK_BINS - 1));
    for (Size i = 0; i < n; ++i) mz[i] = static_cast<float>(spectrum[i].getMZ());

    // Intensity ranks (1 = most intense; ties by ascending m/z). Only the ranks up to the last bin boundary (80)
    // distinguish bins, so a partial sort under the same total order replaces the full stable sort.
    constexpr Size ranked_bins_end = 80;
    std::vector<Size> order(n);
    std::iota(order.begin(), order.end(), Size(0));
    const Size ranked = std::min(n, ranked_bins_end);
    std::partial_sort(order.begin(), order.begin() + ranked, order.end(), [&spectrum](Size a, Size b) {
      const float ia = spectrum[a].getIntensity(), ib = spectrum[b].getIntensity();
      return ia > ib || (ia == ib && a < b);
    });
    for (Size rank = 0; rank < ranked; ++rank) bins[order[rank]] = static_cast<std::uint8_t>(FragmentIonLikelihoodModel::rankBin(rank + 1));

    mz_[index].swap(mz);
    rank_bins_[index].swap(bins);
  }

  Int FragmentIonLikelihoodModel::PeakLists::findNearest(Size index, double mz, double tolerance) const
  {
    const std::vector<float>& list = mz_[index];
    if (list.empty()) return -1;
    // as MSSpectrum::findNearest: the first peak at or above mz, or the one before it if that is nearer
    const auto it = std::lower_bound(list.begin(), list.end(), mz, [](float peak, double value) { return static_cast<double>(peak) < value; });
    Size nearest;
    if (it == list.begin()) { nearest = 0; }
    else if (it == list.end()) { nearest = list.size() - 1; }
    else
    {
      const Size right = static_cast<Size>(it - list.begin());
      nearest = std::fabs(static_cast<double>(list[right]) - mz) < std::fabs(static_cast<double>(list[right - 1]) - mz) ? right : right - 1;
    }
    const double found = list[nearest];
    return (found >= mz - tolerance && found <= mz + tolerance) ? static_cast<Int>(nearest) : -1;
  }

  Size FragmentIonLikelihoodModel::PeakLists::totalPeaks() const
  {
    Size total = 0;
    for (const auto& list : mz_) total += list.size();
    return total;
  }

  Size FragmentIonLikelihoodModel::PeakLists::memoryUsage() const
  {
    Size bytes = mz_.capacity() * sizeof(std::vector<float>) + rank_bins_.capacity() * sizeof(std::vector<std::uint8_t>);
    for (const auto& list : mz_) bytes += list.capacity() * sizeof(float);
    for (const auto& list : rank_bins_) bytes += list.capacity() * sizeof(std::uint8_t);
    return bytes;
  }

  // ---------------------------------------------------------------------------
  // FragmentIonLikelihoodModel
  // ---------------------------------------------------------------------------

  FragmentIonLikelihoodModel::FragmentIonLikelihoodModel(ContextSet context_set, double pseudo_count) :
    context_set_(context_set),
    sites_(context_set == ContextSet::RICH ? CLEAVAGE_SITES : 1),
    complement_(context_set == ContextSet::RICH ? COMPLEMENT_STATES : 1),
    errors_(context_set == ContextSet::RICH ? ERROR_BINS : 1),
    contexts_(SERIES * PRECURSOR_BUCKETS * FRAGMENT_CHARGES * POSITION_BINS * sites_ * complement_),
    outcomes_(RANK_BINS * errors_ + 1),
    signal_counts_(contexts_ * outcomes_, 0.0),
    noise_counts_(contexts_ * outcomes_, 0.0),
    pseudo_count_(pseudo_count)
  {
    if (! std::isfinite(pseudo_count) || pseudo_count <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The pseudo-count must be finite and positive.");
    }
  }

  Size FragmentIonLikelihoodModel::rankBin(Size rank)
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

  Size FragmentIonLikelihoodModel::errorBin(double relative_error)
  {
    return relative_error < 0.25 ? 0 : (relative_error < 0.5 ? 1 : 2);
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
    Size value = 0;
    const Size begin = pos;
    while (pos < name.size() && name[pos] >= '0' && name[pos] <= '9')
    {
      value = value * 10 + static_cast<Size>(name[pos] - '0');
      ++pos;
    }
    if (pos == begin || value == 0) return false;
    ordinal = value;
    return true;
  }

  void FragmentIonLikelihoodModel::matchIons(const PeakLists& peaks,
                                             Size index,
                                             const MSSpectrum& theoretical,
                                             const AASequence& peptide,
                                             int precursor_charge,
                                             double tolerance,
                                             bool ppm,
                                             std::vector<Ion>& ions) const
  {
    if (! std::isfinite(tolerance) || tolerance <= 0.0)
    {
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Matching tolerance must be finite and positive.");
    }
    if (theoretical.getStringDataArrays().empty() || theoretical.getStringDataArrays()[0].size() != theoretical.size())
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Theoretical spectrum needs one ion name per peak", "IonNames");
    }
    const auto& names = theoretical.getStringDataArrays()[0];
    // TheoreticalSpectrumGenerator writes the fragment charges as the first integer array.
    const auto* charges = (! theoretical.getIntegerDataArrays().empty() && theoretical.getIntegerDataArrays()[0].size() == theoretical.size())
                            ? &theoretical.getIntegerDataArrays()[0]
                            : nullptr;
    const Size n = peptide.size();
    const Size bucket = precursor_charge <= 2 ? 0 : (precursor_charge == 3 ? 1 : 2);
    const bool rich = context_set_ == ContextSet::RICH;

    // Pass 1: match every ion; its context waits for the complementary ions (RICH), so it holds the
    // provisional key (ordinal, series, charge bin) for now.
    ions.clear();
    ions.reserve(theoretical.size());
    std::vector<std::uint8_t> matched_singly(rich ? 2 * (n + 1) : 0, 0); // per (series, ordinal): a charge-1 ion matched
    for (Size i = 0; i < theoretical.size(); ++i)
    {
      bool prefix = false;
      Size ordinal = 0;
      if (! parseIonName(names[i], prefix, ordinal) || ordinal >= n) continue;
      const int fragment_charge = charges != nullptr ? (*charges)[i] : 1;
      const Size charge_bin = fragment_charge <= 1 ? 0 : 1;
      const double mz = theoretical[i].getMZ();
      const double width = ppm ? mz * tolerance * 1e-6 : tolerance;
      Size outcome = absentOutcome();
      const Int nearest = peaks.findNearest(index, mz, width);
      if (nearest >= 0)
      {
        const Size peak = static_cast<Size>(nearest);
        const Size error = rich ? errorBin(std::fabs(peaks.mz(index, peak) - mz) / width) : 0;
        outcome = peaks.rankBin(index, peak) * errors_ + error;
        if (rich && charge_bin == 0) matched_singly[(prefix ? 0 : 1) * (n + 1) + ordinal] = 1;
      }
      ions.push_back({static_cast<UInt32>((ordinal << 2) | (prefix ? 0u : 2u) | charge_bin), static_cast<UInt32>(outcome)});
    }

    // Pass 2: full contexts
    for (Ion& ion : ions)
    {
      const Size ordinal = ion.context >> 2;
      const bool prefix = (ion.context & 2u) == 0;
      const Size charge_bin = ion.context & 1u;
      const Size series = prefix ? 0 : 1;
      const Size position = std::min(POSITION_BINS - 1, (POSITION_BINS * ordinal) / n);
      Size site = 0, complement = 0;
      if (rich)
      {
        // the bond between residues nside and nside + 1
        const Size nside = prefix ? ordinal - 1 : n - ordinal - 1;
        const char c_residue = peptide[nside + 1].getOneLetterCode()[0];
        const char n_residue = peptide[nside].getOneLetterCode()[0];
        site = c_residue == 'P' ? 1 : ((n_residue == 'D' || n_residue == 'E') ? 2 : 0);
        complement = matched_singly[(prefix ? 1 : 0) * (n + 1) + (n - ordinal)];
      }
      const Size context = ((((series * PRECURSOR_BUCKETS + bucket) * FRAGMENT_CHARGES + charge_bin) * POSITION_BINS + position) * sites_ + site)
                           * complement_ + complement;
      ion.context = static_cast<UInt32>(context);
    }
  }

  void FragmentIonLikelihoodModel::addObservations(const std::vector<Ion>& ions, bool noise)
  {
    std::vector<double>& counts = noise ? noise_counts_ : signal_counts_;
    for (const Ion& ion : ions)
    {
      counts[ion.context * outcomes_ + ion.outcome] += 1.0;
    }
    if (noise) { ++noise_psms_; }
    else { ++signal_psms_; }
    finalized_ = false;
  }

  Size FragmentIonLikelihoodModel::level1Context_(Size context) const
  {
    // context = (((((series * B + bucket) * C + charge) * P + position) * sites + site) * complement + complement_state)
    const Size rest = context / (complement_ * sites_ * POSITION_BINS);
    const Size charge = rest % FRAGMENT_CHARGES;
    const Size series = rest / (FRAGMENT_CHARGES * PRECURSOR_BUCKETS);
    // BASIC: series; RICH: series and fragment charge
    return context_set_ == ContextSet::RICH ? series * FRAGMENT_CHARGES + charge : series;
  }

  Size FragmentIonLikelihoodModel::level2Context_(Size context) const
  {
    const Size site_complement = context % (complement_ * sites_);
    const Size series_bucket_charge = context / (complement_ * sites_ * POSITION_BINS);
    // BASIC: series, precursor and fragment charge; RICH: in addition cleavage site and complement
    return series_bucket_charge * (sites_ * complement_) + site_complement;
  }

  std::vector<double> FragmentIonLikelihoodModel::smooth_(const std::vector<double>& counts) const
  {
    // Back-off levels above the full context: level 2 (BASIC: series, precursor and fragment charge; RICH: in
    // addition cleavage site and complement), level 1 (BASIC: series; RICH: series and fragment charge) and the
    // global outcome distribution with a flat Laplace prior. Each level's probability blends its own counts with
    // its parent's probability, weighted by the pseudo-count, so sparse cells inherit from their parents and unseen
    // cells stay finite. Counts are integers, so every sum below is exact and independent of summation order.
    const Size level1_size = context_set_ == ContextSet::RICH ? SERIES * FRAGMENT_CHARGES : SERIES;
    const Size level2_size = SERIES * PRECURSOR_BUCKETS * FRAGMENT_CHARGES * sites_ * complement_;
    std::vector<double> global(outcomes_, 0.0), level1(level1_size * outcomes_, 0.0), level2(level2_size * outcomes_, 0.0);
    std::vector<Size> level2_parent(level2_size, 0);
    for (Size context = 0; context < contexts_; ++context)
    {
      const Size l1 = level1Context_(context), l2 = level2Context_(context);
      level2_parent[l2] = l1;
      for (Size o = 0; o < outcomes_; ++o)
      {
        const double c = counts[context * outcomes_ + o];
        global[o] += c;
        level1[l1 * outcomes_ + o] += c;
        level2[l2 * outcomes_ + o] += c;
      }
    }
    const double global_total = std::accumulate(global.begin(), global.end(), 0.0);
    std::vector<double> global_p(outcomes_);
    for (Size o = 0; o < outcomes_; ++o) global_p[o] = (global[o] + 1.0) / (global_total + static_cast<double>(outcomes_));

    // p = (own + pseudo_count * parent) / (total + pseudo_count), cell by cell
    const double pc = pseudo_count_;
    const Size n_out = outcomes_;
    auto blend = [pc, n_out](const double* own, const double* parent, double* p) {
      double total = 0.0;
      for (Size o = 0; o < n_out; ++o) total += own[o];
      for (Size o = 0; o < n_out; ++o) p[o] = (own[o] + pc * parent[o]) / (total + pc);
    };
    std::vector<double> level1_p(level1.size()), level2_p(level2.size());
    for (Size l1 = 0; l1 < level1_size; ++l1) blend(&level1[l1 * n_out], global_p.data(), &level1_p[l1 * n_out]);
    for (Size l2 = 0; l2 < level2_size; ++l2) blend(&level2[l2 * n_out], &level1_p[level2_parent[l2] * n_out], &level2_p[l2 * n_out]);
    std::vector<double> logp(contexts_ * n_out);
    for (Size context = 0; context < contexts_; ++context)
    {
      blend(&counts[context * n_out], &level2_p[level2Context_(context) * n_out], &logp[context * n_out]);
      for (Size o = 0; o < n_out; ++o) logp[context * n_out + o] = std::log(logp[context * n_out + o]);
    }
    return logp;
  }

  void FragmentIonLikelihoodModel::finalize()
  {
    signal_logp_ = smooth_(signal_counts_);
    noise_logp_ = smooth_(noise_counts_);
    finalized_ = true;
  }

  void FragmentIonLikelihoodModel::requireTrained_() const
  {
    if (! finalized_)
    {
      throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "finalize() must be called after the last observation and before scoring.");
    }
  }

  double FragmentIonLikelihoodModel::logLikelihoodRatio(Size context, Size outcome) const
  {
    requireTrained_();
    const Size i = std::min(context, contexts_ - 1) * outcomes_ + std::min(outcome, absentOutcome());
    return signal_logp_[i] - noise_logp_[i];
  }

  double FragmentIonLikelihoodModel::presenceProbability(Size context) const
  {
    requireTrained_();
    return 1.0 - std::exp(signal_logp_[std::min(context, contexts_ - 1) * outcomes_ + absentOutcome()]);
  }

  FragmentIonLikelihoodModel::Features FragmentIonLikelihoodModel::score(const std::vector<Ion>& ions, Size top_k) const
  {
    requireTrained_();
    Features features;
    features.theoretical_ions = ions.size();
    if (ions.empty()) return features;

    const Size absent = absentOutcome();
    double presence_total = 0.0, presence_matched = 0.0;
    std::vector<std::pair<double, bool>> predicted; // (P(present), matched) per ion
    predicted.reserve(ions.size());
    for (const Ion& ion : ions)
    {
      const Size cell = ion.context * outcomes_;
      const bool matched = ion.outcome != absent;
      features.log_likelihood_ratio += signal_logp_[cell + ion.outcome] - noise_logp_[cell + ion.outcome];
      const double presence = 1.0 - std::exp(signal_logp_[cell + absent]);
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
