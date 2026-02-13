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
               std::function<double(const Pep::acceptor_t&)> get);

  /**
   * Train the model with the given data.
   */
  void train(auto&& training);

  /**
   * Generate predictions for the given group and update the
   * acceptor PEP values.
   *
   * Returns `true` if this method was able to make predictions.
   */
  bool predict(Pep::group_t& group);

public:
  std::size_t targets;
  std::size_t decoys;

private:
  constexpr static const double LABEL_TARGET = 1.0;
  constexpr static const double LABEL_DECOY = 0.0;

  // Encode the given acceptor into the given predictor map.
  void encode(const Pep::acceptor_t&, SimpleSVM::PredictorMap&, bool);

private:
  double cutoff;
  std::function<double(const Pep::acceptor_t&)> getter;
  std::size_t index;

  SimpleSVM svm;
  SimpleSVM::PredictorMap training_predictors;
  SimpleSVM::PredictorMap prediction_predictors;
  std::map<std::size_t, double> labels;
};

/******************************************************************************/
Pep::Pep(const std::vector<Acceptor*> acceptors)
{
  auto push_acceptor_t
    = [](Acceptor* acceptor, std::vector<Acceptor::scored_t>& src,
         bool is_target, group_t& dest) -> void {
    for (auto& scored : src)
    {
      dest.push_back(std::make_shared<acceptor_t>(acceptor, scored, is_target));
    }
  };

  for (auto acceptor : acceptors)
  {
    push_acceptor_t(acceptor, acceptor->targets, true, this->acceptors);
    push_acceptor_t(acceptor, acceptor->decoys, false, this->acceptors);
  }
}

/******************************************************************************/
bool Pep::run()
{
  size_t decoy_count
    = std::ranges::count_if(acceptors, std::not_fn(&acceptor_t::is_target));

  if (acceptors.size() > MIN_TOTAL_ACCEPTORS && decoy_count > MIN_DECOYS)
  {
    std::vector<group_t> groups;
    groups.reserve(CROSS_VALIDATION_GROUPS);
    group_acceptors(groups);

    // First round uses MBR score.
    if (! round(groups, [](auto& a) { return a.score.mbr_score; }))
    {
      return false;
    }

    // All other rounds use the PEP.
    for (std::size_t i {}; i < (TRAINING_ROUNDS - 1); ++i)
    {
      if (! round(groups, &acceptor_t::pep_score)) { return false; }
    }

    // FIXME: Filter and set FDR.

    return true;
  }

  return false;
}

/******************************************************************************/
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

  // FIXME: Should we sort these by MBR score too?
  auto decoys = rg::partition(acceptors, &acceptor_t::is_target);
  std::size_t decoys_per_group = split(decoys);

  // Targets are a little different, we need them to be evenly
  // distributed by score so that each group has access to labeled
  // targets (those in the upper 25%).
  auto target_pool = rg::subrange(acceptors.begin(), decoys.begin());
  std::size_t targets_per_group = split(target_pool);

  rg::sort(target_pool, [](auto& a, auto& b) {
    return a->score.mbr_score > b->score.mbr_score;
  });

  std::vector<group_t> target_groups;
  target_groups.reserve(CROSS_VALIDATION_GROUPS);

  for (std::size_t i = 0; i < CROSS_VALIDATION_GROUPS; ++i)
  {
    group_t group;
    group.reserve(targets_per_group);
    target_groups.push_back(group);
  }

  for (std::size_t i {}; auto& target : target_pool)
  {
    target_groups[i].push_back(target);
    i = (i + 1) % CROSS_VALIDATION_GROUPS;
  }

  // View to pull from each of the distributed target groups.
  auto targets = vw::join(target_groups);

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
bool Pep::round(std::vector<group_t>& groups,
                std::function<double(const acceptor_t&)> getter)
{
  for (std::size_t i : std::views::iota(std::size_t {}, groups.size()))
  {
    group_t& group = groups[i];

    // FIXME: Does the sort direction change between MBR and PEP?
    // It seems like it would have to.
    std::ranges::sort(group, [&getter](const auto& a, const auto& b) {
      return getter(*a) > getter(*b);
    });

    std::size_t cutoff_index = std::floor(group.size() * TRUE_POSITIVE_CUTOFF);
    double cutoff = getter(*group[cutoff_index]);

    // This group becomes the prediction set, and all other groups
    // become the training data.
    auto others
      = std::views::iota(std::size_t {}, CROSS_VALIDATION_GROUPS)
        | std::views::filter([&i](auto j) { return i != j; })
        | std::views::transform([&groups](auto n) { return groups[n]; })
        | std::views::join;

    if (! train_predict(others, group, cutoff, getter)) { return false; }
  }

  return true;
}

/******************************************************************************/
bool Pep::train_predict(auto&& training,
                        group_t& predict,
                        double cutoff,
                        std::function<double(const acceptor_t&)> get)
{
  predictors_t predictors(cutoff, get);
  predictors.train(training);
  return predictors.predict(predict);
}

/******************************************************************************/
predictors_t::predictors_t(double cutoff,
                           std::function<double(const Pep::acceptor_t&)> get):
    cutoff(cutoff),
    getter(get)
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
void predictors_t::train(auto&& training)
{
  for (const Pep::acceptor_ptr_t& acceptor : training)
  {
    encode(*acceptor, training_predictors, true);
    svm.setup(training_predictors, labels);
  }
}

/******************************************************************************/
bool predictors_t::predict(Pep::group_t& group)
{
  // Need a certain number of labeled feature or this won't work.
  if (targets <= 4 || decoys <= 4) return false;

  for (const Pep::acceptor_ptr_t& acceptor : group)
  {
    encode(*acceptor, prediction_predictors, false);
  }

  std::vector<SimpleSVM::Prediction> predictions;
  svm.predict(prediction_predictors, predictions);

  // FIXME: Double check this!
  for (std::size_t i {}; i < predictions.size(); ++i)
  {
    auto item = predictions[i].probabilities.find(LABEL_TARGET);

    if (item == predictions[i].probabilities.end())
    {
      OPENMS_LOG_WARN
        << "WARNING: Predictions map does not contain expected labels: ";

      auto view
        = predictions[i].probabilities | std::views::transform([](auto pair) {
            return pair.first;
          }); // Missing C++23's `join_with` :(

      for (auto key : view)
        OPENMS_LOG_WARN << key << " ";

      OPENMS_LOG_WARN << std::endl;
      return false;
    }
    else { group[i]->pep_score = 1 - item->second; }
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
    if (acceptor.is_target && getter(acceptor) >= cutoff)
    {
      labels[index] = LABEL_TARGET;
      ++targets;
    }
    else if (! acceptor.is_target)
    {
      labels[index] = LABEL_DECOY;
      ++decoys;
    }

    ++index;
  }
}

} // namespace OpenMS::PipEcho
