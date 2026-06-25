// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#pragma once

#include "FeatureTypes.h"

#include <functional>
#include <memory>

namespace OpenMS::PipEcho
{

/**
 * Scores match-between-runs (MBR) transfers and applies FDR control to them.
 *
 * Given the candidate target and decoy acceptor transfers, an SVM is trained on
 * the target/decoy labels to assign each transfer a posterior error probability
 * (PEP), from which q-values are derived. run() returns only the transfers that
 * survive the requested FDR cutoff (with a conservative fallback when the FDR is
 * not estimable -- see run()).
 */
class TransferFDRModel
{
public:
  /// Everything we need to know about an acceptor feature.
  struct acceptor_t
  {
    acceptor_t(const FeatureRef& acceptor,
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

    std::reference_wrapper<const FeatureRef> acceptor;
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
  /**
   * Construct the PEP analysis over the given acceptors.
   *
   * @param acceptors   The target and decoy acceptor transfers to score.
   * @param min_decoys  FDR-estimability floor: the minimum number of decoy
   *                    transfers required before any transfer is kept (see
   *                    run()).
   */
  TransferFDRModel(const std::vector<std::shared_ptr<Acceptor>>&, std::size_t min_decoys);

  /**
   * Run the PEP analysis.
   *
   * Returns a vector of scored acceptors that survived FDR control.
   *
   * If the transfer FDR cannot be estimated reliably at the requested cutoff
   * -- too few decoys to resolve it (1/decoys > cutoff) or fewer than
   * `min_decoys` -- ALL transfers are dropped (conservative fallback) and an
   * empty group is returned; the caller keeps only direct identifications. A
   * cutoff of 1.0 disables FDR control and keeps all transfers regardless.
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

  /// FDR-estimability floor (see run()).
  std::size_t min_decoys_;

  /// Whether ion mobility is used as an SVM feature. Ion mobility is an
  /// all-or-nothing property of the acquisition, so this is true only when
  /// *every* acceptor carries a valid IM score; that keeps the SVM predictor
  /// columns rectangular (SimpleSVM requires equal-length predictor vectors).
  bool use_im_ = false;
};

} // namespace OpenMS::PipEcho
