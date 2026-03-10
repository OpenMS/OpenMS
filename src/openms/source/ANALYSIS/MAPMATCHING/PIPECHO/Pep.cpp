// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "OpenMS/CONCEPT/LogStream.h"
#include "OpenMS/ML/SVM/SimpleSVM.h"
#include "Pep.h"
#include "source/ANALYSIS/MAPMATCHING/PIPECHO/Util.h"

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
// A wrapper around the SVM wrapper.
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
               std::function<double(const Pep::acceptor_t&)> get,
               std::function<bool(double, double)> cmp);

  /**
   * Train the model with the given data.
   *
   * Returns `true` if this method was able to train the SVM.
   */
  bool train(const Pep::group_t& training);

  /**
   * Generate predictions for the given group and update the
   * acceptor PEP values.
   *
   * Returns `true` if this method was able to make predictions.
   */
  bool predict(Pep::group_t& group);

public:
  std::size_t targets = 0;
  std::size_t decoys = 0;

private:
  constexpr static const double LABEL_TARGET = 1.0;
  constexpr static const double LABEL_DECOY = 0.0;

  // Encode the given acceptor into the given predictor map.
  void encode(const Pep::acceptor_t&, SimpleSVM::PredictorMap&, bool);

private:
  double cutoff;
  std::function<double(const Pep::acceptor_t&)> getter;
  std::function<bool(double, double)> cmp;
  std::size_t index = 0;

  SimpleSVM svm;
  SimpleSVM::PredictorMap training_predictors;
  SimpleSVM::PredictorMap prediction_predictors;
  std::map<std::size_t, double> labels;
};

/******************************************************************************/
Pep::Pep(const std::vector<std::shared_ptr<Acceptor>>& acceptors)
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
}

/******************************************************************************/
const Pep::group_t& Pep::run(double fdr_cutoff)
{
  bool pep_okay = internal_run();

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

  // One last pass to update the `Feature` target/decoy status and q-value:
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

    // FIXME: Record the q-value.
  }

  return acceptors_;
}

/******************************************************************************/
bool Pep::internal_run()
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
void Pep::group_acceptors(std::vector<group_t>& groups)
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
bool Pep::round(size_t round_number,
                std::vector<group_t>& groups,
                std::function<double(const acceptor_t&)> getter,
                std::function<bool(double, double)> cmp)
{
  OPENMS_LOG_INFO << "PIP-ECHO: Target/decoy train-predict round "
                  << round_number << "/" << TRAINING_ROUNDS << std::endl;

  group_t training; // All groups except the prediction group.
  auto sizes = groups | std::views::transform(&group_t::size);
  training.reserve(std::accumulate(sizes.begin(), sizes.end(), 0));

  for (std::size_t i : std::views::iota(std::size_t {}, groups.size()))
  {
    training.clear();

    for (std::size_t j : std::views::iota(std::size_t {}, groups.size()))
    {
      if (i != j)
      {
        training.insert(training.end(), groups[i].begin(), groups[i].end());
      }
    }

    std::ranges::sort(training, cmp, Util::dref_fn(getter));

    std::size_t cutoff_index
      = std::floor(training.size() * TRUE_POSITIVE_CUTOFF);
    double cutoff = std::invoke(getter, *training[cutoff_index]);

    OPENMS_LOG_INFO << "PIP-ECHO: Acceptor group " << i + 1 << "/"
                    << groups.size() << " will be the prediction group with "
                    << groups[i].size() << " acceptors"
                    << " with the remaining " << training.size()
                    << " acceptors in the training set"
                    << " with selected cutoff value of " << cutoff << std::endl;

    if (! train_predict(training, groups[i], cutoff, getter, cmp))
    {
      return false;
    }
  }

  return true;
}

/******************************************************************************/
bool Pep::train_predict(const group_t& training,
                        group_t& predict,
                        double cutoff,
                        std::function<double(const acceptor_t&)> get,
                        std::function<bool(double, double)> cmp)
{
  predictors_t predictors(cutoff, get, cmp);
  if (! predictors.train(training)) return false;
  return predictors.predict(predict);
}

/******************************************************************************/
void Pep::compute_qvalues(bool pep_values_are_valid)
{
  // Sort ascending by PEP score, and then descending by MBR score.
  std::ranges::sort(acceptors_, [&](auto a, auto b) {
    if (! pep_values_are_valid || std::fabs(a->pep_score - b->pep_score) < 0.01)
    {
      return a->score.mbr_score > b->score.mbr_score;
    }
    else
    {
      return a->pep_score < b->pep_score;
    }
  });

  double decoy_mbr_count {}, decoy_peptide_count {}, decoy_double_count {},
    total_so_far {};

  for (auto acceptor : acceptors_)
  {
    ++total_so_far;

    bool is_peptide_decoy = false; // FIXME
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
void Pep::correct_qvalues()
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
predictors_t::predictors_t(double cutoff,
                           std::function<double(const Pep::acceptor_t&)> get,
                           std::function<bool(double, double)> cmp):
    cutoff(cutoff),
    getter(get),
    cmp(cmp)
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
bool predictors_t::train(const Pep::group_t& training)
{
  for (const Pep::acceptor_ptr_t& acceptor : training)
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
bool predictors_t::predict(Pep::group_t& group)
{
  for (const Pep::acceptor_ptr_t& acceptor : group)
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
void predictors_t::encode(const Pep::acceptor_t& acceptor,
                          SimpleSVM::PredictorMap& predictors,
                          bool create_labels)
{
  predictors["intensity"].push_back(acceptor.score.intensity);
  predictors["rt_diff_error"].push_back(acceptor.score.rt_diff_error);
  predictors["mass_error"].push_back(acceptor.score.mass_error);
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

} // namespace OpenMS::PipEcho
