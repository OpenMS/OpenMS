// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/Percolator.h>
#include <OpenMS/ANALYSIS/ID/PercolatorTypes.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include "TransferFDRModel.h"
#include "Util.h"

#include <algorithm>
#include <exception>
#include <functional>

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

  // Percolator::rescore can throw (e.g. degenerate training data).  Catch it
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
  catch (const std::exception& e)
  {
    // The vendored Percolator code can also throw plain std::exceptions.
    OPENMS_LOG_WARN << "WARNING: PEP computation threw (" << e.what() << ")."
                    << std::endl;
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

  OPENMS_LOG_INFO << "PIP-ECHO: Computing PEP values via Percolator with "
                  << acceptors_.size() << " total acceptors, " << decoy_count
                  << " of which are decoys." << std::endl;

  if (acceptors_.size() < MIN_TOTAL_ACCEPTORS || decoy_count < MIN_DECOYS)
  {
    return false;
  }

  // Semi-supervised target/decoy training (iterative, 3-fold cross-validated
  // linear SVM) via the in-tree Percolator. Percolator starts from the MBR score
  // as the initial direction, trains on the target transfers that pass its
  // training FDR, and returns an SVM score for every transfer (its PEP is not
  // used: the q-values below rank by the SVM score).
  RescoreInput input;
  input.feature_names = {"intensity", "rt_diff_error", "mass_error", "isotope_score"};
  if (use_im_) { input.feature_names.push_back("im_diff_score"); }
  input.feature_names.push_back("mbr_score");

  input.features.reserve(acceptors_.size());
  input.is_decoy.reserve(acceptors_.size());
  // Each transfer is its own row: no cv_group_keys. Rows that share a key also
  // share a "scan", and Percolator's target-decoy competition then keeps only
  // the best row per scan, which would pit a target and a decoy transfer onto
  // the same acceptor feature against each other.

  for (const acceptor_ptr_t& acceptor : acceptors_)
  {
    const Score& s = acceptor->score;
    std::vector<double> row {s.intensity,
                             // The calibrated [0, 1] RT agreement replaces the raw
                             // |delta RT| when every transfer carries one.
                             use_rt_score_ ? s.rt_score : s.rt_diff_error,
                             s.mass_error, s.isotope_score};
    // Ion mobility is used only when every transfer has a score for it.
    if (use_im_) { row.push_back(s.im_diff_score); }
    row.push_back(s.mbr_score);

    input.features.push_back(std::move(row));
    input.is_decoy.push_back(! acceptor->is_target());
  }

  Percolator percolator;
  Param param = percolator.getParameters();
  param.setValue("initial_direction", "mbr_score");
  // 0 = unlimited; Percolator trains on a random subset when the cap is hit.
  param.setValue("subset_max_train", static_cast<int>(max_training_points_));
  percolator.setParameters(param);

  const RescoreOutput output = percolator.rescore(input);

  if (percolator.usedDefaultDirection())
  {
    OPENMS_LOG_INFO << "PIP-ECHO: Percolator could not improve on the MBR score;"
                       " transfers are ranked by the MBR score direction." << std::endl;
  }

  for (std::size_t i = 0; i < acceptors_.size(); ++i)
  {
    acceptors_[i]->svm_score = output.scores[i];
  }

  return true;
}

/******************************************************************************/
void TransferFDRModel::compute_qvalues(bool svm_scores_are_valid)
{
  // Sort descending by the Percolator SVM score, and then descending by MBR
  // score.
  std::ranges::sort(acceptors_, [&](auto a, auto b) {
    if (! svm_scores_are_valid || a->svm_score == b->svm_score)
    {
      return a->score.mbr_score > b->score.mbr_score;
    }
    return a->svm_score > b->svm_score;
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

} // namespace OpenMS::PipEcho
