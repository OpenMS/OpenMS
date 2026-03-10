// Copyright (c) 2025-present, OpenMS Inc. -- EKU Tuebingen
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "PeakTypes.h"

#include <functional>
#include <memory>

namespace OpenMS::PipEcho
{

class Pep
{
public:
  /// Everything we need to know about an acceptor feature.
  struct acceptor_t
  {
    acceptor_t(const Peak& acceptor,
               const Acceptor::scored_t& scored,
               DonorType donor_type):
        acceptor(std::cref(acceptor)),
        donor_ident(scored.second->ident()),
        donor_charge(scored.second->feature.getCharge()),
        score(scored.first),
        donor_type(donor_type) {};

    /// Helper function.
    bool is_target() const
    {
      return donor_type == DonorType::Target;
    };

    /// Helper function.
    double mbr_score() const
    {
      return score.mbr_score;
    };

    std::reference_wrapper<const Peak> acceptor;
    std::string donor_ident;
    BaseFeature::ChargeType donor_charge;
    Score score;
    DonorType donor_type;
    double pep_score = 1.0;
    double q_value = 1.0;
  };

  /// Just for cleaner looking code.
  using acceptor_ptr_t = std::shared_ptr<acceptor_t>;

  /// A group of `acceptor_t` objects.
  using group_t = std::vector<acceptor_ptr_t>;

public:
  Pep(const std::vector<std::shared_ptr<Acceptor>>&);

  /**
   * Run the PEP analysis.
   *
   * Returns a vector of scored acceptors that survived FDR control.
   *
   * Parameters:
   *
   *  - The FDR cutoff value.
   */
  const group_t& run(double);

private:
  /// Try to compute PEP values.
  bool internal_run();

  /// Group the acceptors into roughly equal sized buckets.
  void group_acceptors(std::vector<group_t>&);

  /// One round of training and predicting for the given groups.
  bool round(std::size_t,
             std::vector<group_t>&,
             std::function<double(const acceptor_t&)>,
             std::function<bool(double, double)>);

  /// Train a model using the data in the first view, and generate
  /// predictions for the data in the second parameter..
  bool train_predict(const group_t& training,
                     group_t& predict,
                     double cutoff,
                     std::function<double(const acceptor_t&)> get,
                     std::function<bool(double, double)> cmp);

  /// Compute and assign Q values to each acceptor.
  void compute_qvalues(bool);

  /// Standard Q value correction.
  void correct_qvalues();

private:
  group_t acceptors_;
};

} // namespace OpenMS::PipEcho
