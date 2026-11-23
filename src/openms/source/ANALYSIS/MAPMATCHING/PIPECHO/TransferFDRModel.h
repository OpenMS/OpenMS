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
 * Given the candidate target and decoy acceptor transfers, the in-tree
 * Percolator (semi-supervised, cross-validated linear SVM) is trained on the
 * target/decoy labels to assign each transfer an SVM score, from whose ranking
 * q-values are derived. run() returns only the transfers that survive the
 * requested FDR cutoff (with a conservative fallback when the FDR is not
 * estimable -- see run()).
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

    std::reference_wrapper<const FeatureRef> acceptor;
    std::string donor_ident;
    BaseFeature::ChargeType donor_charge;
    Score score;
    DonorType donor_type;
    double svm_score = 0.0; ///< Percolator SVM score (higher = more target-like)
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
   * @param max_training_points  Upper bound on the number of transfers the
   *                    SVM is trained on (0 = unlimited; Percolator's
   *                    'subset_max_train'). Only the model fit is subsampled;
   *                    every transfer is scored and enters the q-values.
   */
  TransferFDRModel(const std::vector<std::shared_ptr<Acceptor>>&, std::size_t min_decoys,
                   std::size_t max_training_points = 0);

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

  /// Compute and assign Q values to each acceptor, ranked by the SVM score when
  /// it is valid, else by the MBR score.
  void compute_qvalues(bool svm_scores_are_valid);

  /// Standard Q value correction.
  void correct_qvalues();

private:
  group_t acceptors_;

  /// FDR-estimability floor (see run()).
  std::size_t min_decoys_;

  /// Upper bound on the SVM training-set size (0 = unlimited). Only the model
  /// fit is subsampled; scoring and FDR use all transfers.
  std::size_t max_training_points_;

  /// Whether ion mobility is used as an SVM feature. Ion mobility is an
  /// all-or-nothing property of the acquisition, so this is true only when
  /// *every* acceptor carries a valid IM score; that keeps the SVM predictor
  /// feature rows rectangular (Percolator requires equal-length rows).
  bool use_im_ = false;

  /// Whether the RT predictor is the CALIBRATED rt_score ([0,1] agreement) rather
  /// than the raw rt_diff_error. Like ion mobility this is all-or-nothing: true
  /// only when every acceptor carries a calibrated rt_score (>= 0), which happens
  /// exactly on the local adaptive RT path with a calibrated RtScoreMode. The
  /// baseline/global path leaves rt_score at the -1 sentinel -> raw rt_diff_error.
  bool use_rt_score_ = false;
};

} // namespace OpenMS::PipEcho
