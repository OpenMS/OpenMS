// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/ML/SVM/SimpleSVM.h>
#include "TransferFDRModel.h"
#include "Util.h"

#include <algorithm>
#include <functional>
#include <numeric>
#include <ranges>

namespace OpenMS::PipEcho
{

/******************************************************************************/
// Minimum number of acceptors (targets + decoys) needed to run PEP
// analysis.
const static size_t MIN_TOTAL_ACCEPTORS = 100;

/******************************************************************************/
// Minimum number of decoys needed to run PEP analysis.
const static size_t MIN_DECOYS = 20;

/******************************************************************************/
// The number of groups to use for cross validation.
const static size_t CROSS_VALIDATION_GROUPS = 3;

/******************************************************************************/
// The number of training and predicting iterations to use.
const static size_t TRAINING_ROUNDS = 10;

/******************************************************************************/
// The percentage of acceptors that will be used as true positives,
// as selected from a vector of acceptors sorted by score.
const static double TRUE_POSITIVE_CUTOFF = 0.25;

/******************************************************************************/
// A wrapper around the SVM wrapper.  File-local (internal linkage); not part of
// the module's header API.
namespace
{
struct predictors_t
{
public:
  /**
   * Construct a new SVM helper class.
   *
   * @param cutoff Don't label features below this value.
   * @param get A function to call to get the score.
   */
  predictors_t(double cutoff,
               std::function<double(const TransferFDRModel::acceptor_t&)> get,
               std::function<bool(double, double)> cmp,
               bool use_im,
               bool use_rt_score);

  /**
   * Train the model with the given data.
   *
   * Returns `true` if this method was able to train the SVM.
   */
  bool train(const TransferFDRModel::group_t& training);

  /**
   * Generate predictions for the given group and update the
   * acceptor PEP values.
   *
   * Returns `true` if this method was able to make predictions.
   */
  bool predict(TransferFDRModel::group_t& group);

public:
  std::size_t targets = 0;
  std::size_t decoys = 0;

private:
  constexpr static const double LABEL_TARGET = 1.0;
  constexpr static const double LABEL_DECOY = 0.0;

  // Encode the given acceptor into the given predictor map.
  void encode(const TransferFDRModel::acceptor_t&, SimpleSVM::PredictorMap&, bool);

private:
  double cutoff;
  std::function<double(const TransferFDRModel::acceptor_t&)> getter;
  std::function<bool(double, double)> cmp;
  bool use_im;
  bool use_rt_score;
  std::size_t index = 0;

  SimpleSVM svm;
  SimpleSVM::PredictorMap training_predictors;
  SimpleSVM::PredictorMap prediction_predictors;
  std::map<std::size_t, double> labels;
};
} // namespace

/******************************************************************************/
TransferFDRModel::TransferFDRModel(const std::vector<std::shared_ptr<Acceptor>>& acceptors,
         std::size_t min_decoys, std::size_t max_training_points):
    min_decoys_(min_decoys),
    max_training_points_(max_training_points)
{
  auto push_acceptor
    = [&](const Acceptor& acceptor, std::optional<Acceptor::scored_t>& scored,
          DonorType dtype) -> void {
    if (scored.has_value())
    {
      acceptors_.push_back(
        std::make_shared<acceptor_t>(acceptor, *scored, dtype));
    }
  };

  for (const auto& acceptor : acceptors)
  {
    push_acceptor(*acceptor, acceptor->target, DonorType::Target);
    push_acceptor(*acceptor, acceptor->decoy, DonorType::Decoy);
  }

  // Ion mobility is used as an SVM feature only when the whole experiment
  // carries it -- i.e. every acceptor has a valid (non-negative) IM score.
  // This keeps the SVM predictor columns rectangular (mixed inputs fall back
  // to the ion-mobility-free feature set).
  use_im_ = ! acceptors_.empty()
            && std::ranges::all_of(acceptors_, [](const acceptor_ptr_t& a) {
                 return a->score.im_diff_score >= 0.0;
               });

  // The calibrated RT-agreement score replaces the raw rt_diff_error predictor
  // only when EVERY acceptor carries one (>= 0); otherwise the column would mix
  // calibrated [0,1] scores with raw seconds. This is true exactly on the local
  // adaptive RT path with a calibrated RtScoreMode; the baseline leaves the -1
  // sentinel, so the raw rt_diff_error is used and the output stays unchanged.
  use_rt_score_ = ! acceptors_.empty()
                  && std::ranges::all_of(acceptors_, [](const acceptor_ptr_t& a) {
                       return a->score.rt_score >= 0.0;
                     });
}

/******************************************************************************/
const TransferFDRModel::group_t& TransferFDRModel::run(double fdr_cutoff)
{
  // Conservative FDR-estimability gate.  Computed up front -- before any
  // throw-prone work and independent of whether the ML model trains -- so an
  // SVM hiccup or a small (but FDR-resolvable) experiment still falls through
  // to the valid MBR-score fallback rather than being dropped here.
  //
  // We drop ALL transfers when the requested FDR cannot be resolved by the
  // available decoys: required_decoys = max(min_decoys_, ceil(1/fdr_cutoff)),
  // expressed below without ceil/division rounding traps.  A cutoff of 1.0
  // disables FDR control (a q-value of 1.0 accepts everything), so transfers
  // are always kept then.  We count MBR decoy transfers -- the kind of decoy
  // the transfer FDR is estimated from.
  const std::size_t decoy_count
    = std::ranges::count_if(acceptors_, std::not_fn(&acceptor_t::is_target));
  const std::size_t target_count = acceptors_.size() - decoy_count;

  const bool unresolvable
    = fdr_cutoff <= 0.0 || decoy_count < min_decoys_
      || static_cast<double>(decoy_count) * fdr_cutoff < 1.0;

  if (target_count > 0 && fdr_cutoff < 1.0 && unresolvable)
  {
    OPENMS_LOG_WARN
      << "WARNING: PIP-ECHO: only " << decoy_count
      << " MBR decoy transfer(s) were generated, which cannot estimate the "
      << "match-between-runs transfer FDR at the requested cutoff " << fdr_cutoff
      << " (a reliable estimate needs at least " << min_decoys_
      << " decoy(s) and enough decoys to resolve the cutoff, i.e. "
         "1/decoys <= fdr). "
      << "The transfer FDR could NOT be controlled. As a conservative "
      << "fallback, all " << target_count
      << " transferred (acceptor) feature(s) are being DROPPED; all direct "
      << "identifications (donors) are RETAINED. Generate more decoys (more "
      << "runs/identifications), raise 'fdr', or lower 'min_decoys' to enable "
      << "MBR transfers." << std::endl;

    // Drop every acceptor (target AND decoy) and return empty.  This is safe:
    // Impl::generate_consensus_map builds the ConsensusMap purely from the
    // donor ident_map (direct identifications), which never reads acceptor_t,
    // so the skipped compute_qvalues()/status-update loop below has no
    // observable effect.
    acceptors_.clear();
    return acceptors_;
  }

  bool pep_okay = false;

  // SimpleSVM::setup/predict can throw (e.g. too few label classes).  Catch it
  // so we fall back to MBR scores instead of aborting the whole tool.
  try
  {
    pep_okay = internal_run();
  }
  catch (const Exception::BaseException& e)
  {
    OPENMS_LOG_WARN << "WARNING: PEP computation threw (" << e.getName() << ": "
                    << e.what() << ")." << std::endl;
    pep_okay = false;
  }

  if (! pep_okay)
  {
    OPENMS_LOG_WARN << "WARNING: unable to calculate ML PEP values, "
                    << "falling back to MBR scores." << std::endl;
  }

  compute_qvalues(pep_okay);

  const auto [first, last]
    = std::ranges::remove_if(acceptors_, [&fdr_cutoff](const auto& a) {
        return a->q_value > fdr_cutoff;
      });

  std::size_t before = acceptors_.size();
  acceptors_.erase(first, last);

  OPENMS_LOG_INFO << "PIP-ECHO: Filtering acceptors at FDR " << fdr_cutoff
                  << " reduced acceptors from " << before << " to "
                  << acceptors_.size() << std::endl;

  // One last pass to update the `Feature` target/decoy status:
  for (auto& acceptor : acceptors_)
  {
    // Only update the status for known decoys.  We don't really
    // "know" if a target is real or not so just leave it alone.
    if (! acceptor->is_target())
    {
      auto hit = Util::feature_hit(acceptor->acceptor.get().feature);

      if (hit.has_value())
      {
        hit->setTargetDecoyType(PeptideHit::TargetDecoyType::DECOY);
      }
    }
  }

  // The FDR was estimable (we passed the gate above), but the fallback filter
  // may still have removed every transfer.  Report that loudly rather than
  // leaving only the soft "falling back to MBR scores" note, so silent
  // near-total transfer loss is honest.
  const std::size_t kept_targets
    = std::ranges::count_if(acceptors_, &acceptor_t::is_target);
  if (kept_targets == 0 && target_count > 0)
  {
    OPENMS_LOG_WARN
      << "WARNING: PIP-ECHO: the match-between-runs fallback retained 0 of "
      << target_count << " candidate transfer(s) at FDR " << fdr_cutoff
      << "; all transferred features were filtered out. Only direct "
      << "identifications (donors) are retained." << std::endl;
  }

  return acceptors_;
}

/******************************************************************************/
bool TransferFDRModel::internal_run()
{
  size_t decoy_count
    = std::ranges::count_if(acceptors_, std::not_fn(&acceptor_t::is_target));

  OPENMS_LOG_INFO << "PIP-ECHO: Computing PEP values via SVM with "
                  << acceptors_.size() << " total acceptors, " << decoy_count
                  << " of which are decoys." << std::endl;

  if (acceptors_.size() < MIN_TOTAL_ACCEPTORS || decoy_count < MIN_DECOYS)
  {
    return false;
  }

  std::vector<group_t> groups;
  groups.reserve(CROSS_VALIDATION_GROUPS);
  group_acceptors(groups);

  // First round uses MBR score in descending order.
  if (! round(1, groups, &acceptor_t::mbr_score, std::ranges::greater {}))
  {
    return false;
  }

  // All other rounds use the PEP in ascending order.
  for (std::size_t i {}; i < (TRAINING_ROUNDS - 1); ++i)
  {
    if (! round(i + 1, groups, &acceptor_t::pep_score, std::ranges::less {}))
    {
      return false;
    }
  }

  return true;
};

/******************************************************************************/
/**
 * Separate features into several groups in such a way as to ensure an
 * equal distribution of high scoring targets.
 */
void TransferFDRModel::group_acceptors(std::vector<group_t>& groups)
{
  // This code would be so much nicer if we could use C++23.
  namespace rg = std::ranges;
  namespace vw = std::views;

  auto split = [](const auto& view) -> std::size_t {
    std::size_t n = view.size();
    std::size_t r = n / CROSS_VALIDATION_GROUPS;
    if (n % CROSS_VALIDATION_GROUPS) ++r;
    return r;
  };

  auto chunk = [](const auto& view, std::size_t size, std::size_t i) {
    return view | vw::drop(size * i) | vw::take(size);
  };

  auto decoys = rg::partition(acceptors_, &acceptor_t::is_target);
  std::size_t decoys_per_group = split(decoys);

  // Targets are a little different, we need them to be evenly
  // distributed by score so that each group has access to labeled
  // targets (those in the upper 25%).
  auto target_pool = rg::subrange(acceptors_.begin(), decoys.begin());
  std::size_t targets_per_group = split(target_pool);
  rg::sort(target_pool, std::greater {}, Util::dref_fn(&acceptor_t::mbr_score));

  std::vector<group_t> target_groups;
  target_groups.reserve(CROSS_VALIDATION_GROUPS);

  // Prepare groups.
  for (std::size_t i = 0; i < CROSS_VALIDATION_GROUPS; ++i)
  {
    group_t group;
    group.reserve(targets_per_group);
    target_groups.push_back(group);
  }

  // Distribute target features into the temporary `target_groups`
  // vector, which will then later be used to fill in the final
  // `groups` vector.
  for (std::size_t i {}; auto& target : target_pool)
  {
    target_groups[i].push_back(target);
    i = (i + 1) % CROSS_VALIDATION_GROUPS;
  }

  // View to pull from each of the distributed target groups.
  auto targets = vw::join(target_groups);

  // Now we can fill in the `groups` vectors.
  for (std::size_t i : vw::iota(std::size_t {}, CROSS_VALIDATION_GROUPS))
  {
    group_t group;
    group.reserve(decoys_per_group + targets_per_group);

    // NB: Targets go first to preserve sorting order.
    auto g_targets = chunk(targets, targets_per_group, i);
    rg::copy(g_targets, std::back_inserter(group));

    auto g_decoys = chunk(decoys, decoys_per_group, i);
    rg::copy(g_decoys, std::back_inserter(group));

    groups.push_back(group);
  }
}

/******************************************************************************/
/// NOTE: It is undefined behavior to use a comparison function that
/// doesn't conform to the *Compare* requirement (i.e. uses `>=` or
/// `<=`).
///
/// https://cppreference.com/w/cpp/named_req/Compare.html
bool TransferFDRModel::round(size_t round_number,
                std::vector<group_t>& groups,
                std::function<double(const acceptor_t&)> getter,
                std::function<bool(double, double)> cmp)
{
  OPENMS_LOG_INFO << "PIP-ECHO: Target/decoy train-predict round "
                  << round_number << "/" << TRAINING_ROUNDS << std::endl;

  group_t training; // All groups except the prediction group.
  auto sizes = groups | std::views::transform(&group_t::size);
  training.reserve(std::accumulate(sizes.begin(), sizes.end(), 0));

  // Deterministic stratified subsample used to cap the per-fold SVM training
  // set at `max_training_points_` (0 = unlimited). Keeps all decoys (scarce,
  // they define the FDR null) and all labelled positives (targets above the
  // cutoff -- the positive training signal), then fills any remaining budget
  // with a score-spread sample of the below-cutoff targets so the predictor
  // scaling stays representative. No RNG -> reproducible. Prediction and the
  // FDR/q-values use EVERY acceptor regardless, so q-values stay unbiased.
  auto subsample_training =
    [](const group_t& full, double cutoff,
       const std::function<double(const acceptor_t&)>& getter,
       const std::function<bool(double, double)>& cmp, std::size_t cap) -> group_t {
    group_t decoys, positives, rest; // preserve the incoming score-sorted order
    for (const acceptor_ptr_t& a : full)
    {
      if (! a->is_target()) { decoys.push_back(a); }
      else if (std::invoke(cmp, std::invoke(getter, *a), cutoff)) { positives.push_back(a); }
      else { rest.push_back(a); }
    }
    // Append `k` score-spread elements of `src` (uniform stride across its full
    // range), or all of `src` when k >= src.size().
    auto append_stride = [](group_t& dst, const group_t& src, std::size_t k) {
      const std::size_t n = src.size();
      if (k == 0 || n == 0) { return; }
      if (k >= n) { dst.insert(dst.end(), src.begin(), src.end()); return; }
      for (std::size_t j = 0; j < k; ++j) { dst.push_back(src[(j * n) / k]); }
    };
    // Spend the budget on the labelled classes first: decoys + positives are
    // what actually train the SVM (the below-cutoff 'rest' is unlabelled and
    // only affects predictor scaling). Keep BOTH classes viable -- a fold with
    // too few labelled positives OR decoys cannot train (predictors_t::train
    // needs > 4 of each) and would silently fall back to MBR scores. When the
    // labelled set does not fit, keep the MINORITY class whole (it almost always
    // fits) and give the remaining budget to the majority; only when both
    // classes individually exceed half the cap do we split evenly. A
    // size-proportional split would starve a small minority (e.g. few decoys
    // against many positives), defeating the cap.
    std::size_t dec_keep = decoys.size();
    std::size_t pos_keep = positives.size();
    if (decoys.size() + positives.size() > cap)
    {
      const std::size_t minority = std::min(decoys.size(), positives.size());
      const std::size_t minority_keep = (minority <= cap / 2) ? minority : cap / 2;
      const std::size_t majority_keep = cap - minority_keep;
      if (decoys.size() <= positives.size()) { dec_keep = minority_keep; pos_keep = majority_keep; }
      else { pos_keep = minority_keep; dec_keep = majority_keep; }
    }
    group_t sampled;
    sampled.reserve(cap);
    append_stride(sampled, decoys, dec_keep);
    append_stride(sampled, positives, pos_keep);
    append_stride(sampled, rest, cap - sampled.size());
    return sampled;
  };

  for (std::size_t i : std::views::iota(std::size_t {}, groups.size()))
  {
    training.clear();

    for (std::size_t j : std::views::iota(std::size_t {}, groups.size()))
    {
      if (i != j)
      {
        training.insert(training.end(), groups[j].begin(), groups[j].end());
      }
    }

    std::ranges::sort(training, cmp, Util::dref_fn(getter));

    std::size_t cutoff_index
      = std::floor(training.size() * TRUE_POSITIVE_CUTOFF);
    double cutoff = std::invoke(getter, *training[cutoff_index]);

    // Cap the training set (cutoff is computed on the full set above, so the
    // labels stay consistent with the unsubsampled cutoff).
    const std::size_t full_training_size = training.size();
    if (max_training_points_ != 0 && training.size() > max_training_points_)
    {
      training = subsample_training(training, cutoff, getter, cmp, max_training_points_);
    }

    OPENMS_LOG_INFO << "PIP-ECHO: Acceptor group " << i + 1 << "/"
                    << groups.size() << " will be the prediction group with "
                    << groups[i].size() << " acceptors, training on "
                    << training.size() << " of " << full_training_size
                    << " acceptors (cutoff value " << cutoff << ")" << std::endl;

    if (! train_predict(training, groups[i], cutoff, getter, cmp))
    {
      return false;
    }
  }

  return true;
}

/******************************************************************************/
bool TransferFDRModel::train_predict(const group_t& training,
                        group_t& predict,
                        double cutoff,
                        std::function<double(const acceptor_t&)> get,
                        std::function<bool(double, double)> cmp)
{
  predictors_t predictors(cutoff, get, cmp, use_im_, use_rt_score_);
  if (! predictors.train(training)) return false;
  return predictors.predict(predict);
}

/******************************************************************************/
void TransferFDRModel::compute_qvalues(bool pep_values_are_valid)
{
  // Sort ascending by PEP score, and then descending by MBR score.
  std::ranges::sort(acceptors_, [&](auto a, auto b) {
    if (! pep_values_are_valid || a->pep_score == b->pep_score)
    {
      return a->score.mbr_score > b->score.mbr_score;
    }
    return a->pep_score < b->pep_score;
  });

  double decoy_mbr_count {}, decoy_peptide_count {}, decoy_double_count {},
    total_so_far {};

  for (const auto& acceptor : acceptors_)
  {
    ++total_so_far;

    bool is_peptide_decoy
      = Util::feature_is_decoy(acceptor->acceptor.get().feature);
    bool is_mbr_decoy = ! acceptor->is_target();

    if (is_peptide_decoy && ! is_mbr_decoy) { ++decoy_peptide_count; }
    else if (! is_peptide_decoy && is_mbr_decoy) { ++decoy_mbr_count; }
    else if (is_peptide_decoy && is_mbr_decoy) { ++decoy_double_count; }

    double errors = std::fmax(0, decoy_peptide_count - decoy_double_count);
    acceptor->q_value = (1 + decoy_mbr_count + errors) / total_so_far;
  }

  correct_qvalues();
}

/******************************************************************************/
// Standard Q value correct.  Ensures that as you iterate over the
// list of scored acceptors the Q value either increases or stays the
// same.
void TransferFDRModel::correct_qvalues()
{
  if (acceptors_.size() < 2) return;
  std::size_t i = acceptors_.size() - 2;

  // If a Q value is greater than the one that comes after it we force
  // it to be the same value as the one that comes after it.  This
  // ensures that Q values increase or stay the same.
  for (;;)
  {
    if (acceptors_[i]->q_value > acceptors_[i + 1]->q_value)
    {
      acceptors_[i]->q_value = acceptors_[i + 1]->q_value;
    }

    if (i == 0) { break; }
    else
    {
      --i;
    }
  }
}

/******************************************************************************/
namespace
{
predictors_t::predictors_t(double cutoff,
                           std::function<double(const TransferFDRModel::acceptor_t&)> get,
                           std::function<bool(double, double)> cmp,
                           bool use_im,
                           bool use_rt_score):
    cutoff(cutoff),
    getter(get),
    cmp(cmp),
    use_im(use_im),
    use_rt_score(use_rt_score)
{
  Param svm_param = svm.getParameters();

  // FIXME: Are these good values?
  svm_param.setValue("kernel", "linear");
  svm_param.setValue("log2_C", ListUtils::create<double>("-5,-1,1,5,7,11,15"));
  svm_param.setValue("log2_p",
                     ListUtils::create<double>(
                       "-15,-9,-6,-3.32192809489,0,3.32192809489,6,9,15"));

  svm.setParameters(svm_param);
}

/******************************************************************************/
bool predictors_t::train(const TransferFDRModel::group_t& training)
{
  for (const TransferFDRModel::acceptor_ptr_t& acceptor : training)
  {
    encode(*acceptor, training_predictors, true);
  }

  // Need a certain number of labeled feature or this won't work.
  if (targets <= 4 || decoys <= 4)
  {
    OPENMS_LOG_WARN << "WARNING: there are not enough targets (" << targets
                    << ") or decoys (" << decoys << ") for proper training."
                    << std::endl;

    return false;
  }

  svm.setup(training_predictors, labels);
  return true;
}

/******************************************************************************/
bool predictors_t::predict(TransferFDRModel::group_t& group)
{
  for (const TransferFDRModel::acceptor_ptr_t& acceptor : group)
  {
    encode(*acceptor, prediction_predictors, false);
  }

  std::vector<SimpleSVM::Prediction> predictions;
  svm.predict(prediction_predictors, predictions);

  for (std::size_t i {}; i < predictions.size(); ++i)
  {
    auto item = predictions[i].probabilities.find(LABEL_TARGET);

    if (item == predictions[i].probabilities.end())
    {
      OPENMS_LOG_WARN
        << "WARNING: Predictions map does not contain expected labels: ";

      // Missing C++23's `join_with` :(
      for (auto key : predictions[i].probabilities | std::views::keys)
        OPENMS_LOG_WARN << key << " ";

      OPENMS_LOG_WARN << std::endl;
      return false;
    }
    else
    {
      group[i]->pep_score = 1 - item->second;
    }
  }

  return true;
}

/******************************************************************************/
void predictors_t::encode(const TransferFDRModel::acceptor_t& acceptor,
                          SimpleSVM::PredictorMap& predictors,
                          bool create_labels)
{
  predictors["intensity"].push_back(acceptor.score.intensity);
  // RT predictor: the calibrated [0,1] agreement score on the local path (when
  // every acceptor has one), else the raw |Δrt| (baseline, byte-identical). The
  // column key is kept stable; only the value's meaning/scale changes with mode.
  predictors["rt_diff_error"].push_back(
    use_rt_score ? acceptor.score.rt_score : acceptor.score.rt_diff_error);
  predictors["mass_error"].push_back(acceptor.score.mass_error);
  // Isotope-distribution agreement is always [0, 1] (0.5 sentinel when absent),
  // so the column is rectangular and can be pushed unconditionally.
  predictors["isotope_score"].push_back(acceptor.score.isotope_score);
  // Ion mobility is added as a predictor only when used experiment-wide, so the
  // predictor columns stay equal-length (required by SimpleSVM).
  if (use_im)
  {
    predictors["im_diff_score"].push_back(acceptor.score.im_diff_score);
  }
  predictors["mbr_score"].push_back(acceptor.score.mbr_score);

  if (create_labels)
  {
    switch (acceptor.donor_type)
    {
      case DonorType::Target:
        if (std::invoke(cmp, std::invoke(getter, acceptor), cutoff))
        {
          labels[index] = LABEL_TARGET;
          ++targets;
        }
        break;

      case DonorType::Decoy:
        labels[index] = LABEL_DECOY;
        ++decoys;
        break;
    }

    ++index;
  }
}
} // namespace

} // namespace OpenMS::PipEcho
