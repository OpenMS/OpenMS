// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Peter J. Jones $
// $Authors: Peter J. Jones $
// --------------------------------------------------------------------------

#include "FeatureTypes.h"
#include "Impl.h"
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include "TransferFDRModel.h"
#include "Score.h"
#include "Util.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <memory>
#include <optional>
#include <ranges>
#include <set>
#include <unordered_map>

namespace OpenMS::PipEcho
{

// File-local helpers (internal linkage); not part of the module's header API.
namespace
{

/******************************************************************************/
/**
 * Return a key that can be used to link features with the same
 * amino acid sequence and charge state.
 */
std::string feature_sequence_key(const std::string& ident,
                                 BaseFeature::ChargeType charge)
{
  return ident + "_" + std::to_string(charge);
}

/******************************************************************************/
/**
 * Return a key that can be used to link features with the same
 * amino acid sequence and charge state.
 */
std::string feature_sequence_key(const Feature& feature)
{
  auto hit = PipEcho::Util::feature_hit(feature);

  if (hit)
  {
    return feature_sequence_key(hit->getSequence().toString(),
                                feature.getCharge());
  }

  std::string msg("donor feature missing peptide sequence");
  throw(Exception::MissingInformation(__FILE__, __LINE__,
                                      OPENMS_PRETTY_FUNCTION, msg));
}

/******************************************************************************/
/// Ensure that a feature is ready for analysis.
void prepare_feature(Feature& feature)
{
  // Ensure everything is sorted for later use.
  for (auto& peptide : feature.getPeptideIdentifications())
  {
    peptide.sort();
  }

  feature.sortPeptideIdentifications();
}

} // namespace

/******************************************************************************/
// [LOCAL-RT] Per-(donor_run, acceptor_run) local RT calibration built from
// peptides identified (=donors) in BOTH runs. Models the residual local RT error
// remaining AFTER the global b-spline alignment (the quantity the single global
// window -- sized by the worst-case anchor residual -- fails to localise).
// predict() re-centers on the locally predicted acceptor RT and sizes the
// half-window from the local residual scatter; sparse queries use a tight fixed
// fallback (never the global window). Enabled by the 'local_rt:enabled' parameter.
// Built per pair in the single-threaded matching loop (see group()).
// Inspired by FlashLfqEngine.PredictRetentionTime: gather up to
// `n_anchors` shared anchors per side WITHIN `locality_s` of the donor RT,
// re-center on their median RT shift, and size the half-window as c*stddev,
// HARD-CAPPED at cap_s (FlashLFQ MaxMbrRtWindow = 1.0 min => +/-30 s). Sparse
// (<2 local anchors) -> a tight fixed fallback window, NEVER the global window.
// This is the key fix: PIP-ECHO's window was sized from the global worst-case
// alignment residual (here +/-1206 s); FlashLFQ never searches wider than ~30 s.
struct LocalRtModel
{
  std::vector<double> anchor_rt;   ///< donor-run RT of shared anchors (ascending)
  std::vector<double> residual;    ///< acceptor_run_RT - donor_run_RT (parallel)
  // Defaults are FlashLFQ-INSPIRED but GENEROUS, not the literal constants:
  // FlashLFQ's +/-30 s cap assumes well-aligned bulk data + a calibrated RT
  // score in a geomean. PIP-ECHO feeds RAW rt to an SVM on poorly-aligned
  // single-cell data, where +/-30 s over-tightens (censors the RT dynamic range
  // and starves decoys). So cap/locality are loosened; the window stays data-
  // driven (c*sigma) and the cap is only a backstop. See adaptive widening below.
  double locality_s = 120.0;       ///< anchors within +/- this of donor RT (FlashLFQ 0.5 min, loosened)
  std::size_t n_anchors = 3;       ///< max anchors per side (FlashLFQ NumberOfAnchorPeptidesForMbr)
  double c = 3.0;                  ///< half-window = c*stddev (FlashLFQ width 6*SD => 3*SD half)
  double floor_s = 5.0;            ///< min half-window (guard degenerate near-zero scatter)
  double cap_s = 100.0;            ///< generous half-window backstop (NOT FlashLFQ's 30 s)
  double fallback_half_s = 15.0;   ///< half-window when <2 local anchors (FlashLFQ 0.25 min)
  double widen = 1.0;              ///< adaptive-widening multiplier (decoy estimability)
  double hard_ceiling = 1e9;       ///< absolute max half-window (= global rt window)

  /// supported = had >=2 local anchors (sigma-sized window); else a tight fallback.
  struct Pred { double center; double half_window; bool supported; };

  /// Apply the adaptive-widening factor, bounded by the global ceiling.
  double scale(double half) const { return std::min(half * widen, hard_ceiling); }

  /// Nearest <=n_anchors-per-side residuals within locality_s of r, leave-one-out
  /// (the query donor's own anchor, at ~r, is skipped). Used by predict() and by
  /// the auto-parameter estimation.
  std::vector<double> gather(double r) const
  {
    const std::size_t n = anchor_rt.size();
    std::vector<double> res;
    if (n == 0) { return res; }
    const std::size_t idx = static_cast<std::size_t>(
      std::lower_bound(anchor_rt.begin(), anchor_rt.end(), r) - anchor_rt.begin());
    constexpr double self_eps = 1e-6;
    res.reserve(2 * n_anchors);
    for (std::size_t i = idx, k = 0; i < n && k < n_anchors; ++i)  // forward (>= r)
    {
      const double d = std::fabs(anchor_rt[i] - r);
      if (d > locality_s) { break; }
      if (d <= self_eps) { continue; }
      res.push_back(residual[i]); ++k;
    }
    for (std::size_t i = idx, k = 0; i > 0 && k < n_anchors; )     // backward (< r)
    {
      --i;
      const double d = std::fabs(r - anchor_rt[i]);
      if (d > locality_s) { break; }
      if (d <= self_eps) { continue; }
      res.push_back(residual[i]); ++k;
    }
    return res;
  }

  /// Sample standard deviation of the local RT shifts at r (nullopt if <2 local
  /// anchors). [LOCAL-RT auto] used to estimate the global window scales.
  std::optional<double> local_sigma_at(double r) const
  {
    std::vector<double> res = gather(r);
    if (res.size() < 2) { return std::nullopt; }
    double mean = 0.0;
    for (double v : res) { mean += v; }
    mean /= static_cast<double>(res.size());
    double var = 0.0;
    for (double v : res) { var += (v - mean) * (v - mean); }
    var /= static_cast<double>(res.size() - 1);
    const double sd = std::sqrt(var);
    return std::isfinite(sd) ? std::optional<double>(sd) : std::nullopt;
  }

  /// Median gap between consecutive anchors (0 if <2). [LOCAL-RT auto]
  double median_anchor_gap() const
  {
    if (anchor_rt.size() < 2) { return 0.0; }
    std::vector<double> gaps;
    gaps.reserve(anchor_rt.size() - 1);
    for (std::size_t i = 1; i < anchor_rt.size(); ++i) { gaps.push_back(anchor_rt[i] - anchor_rt[i - 1]); }
    std::sort(gaps.begin(), gaps.end());
    return gaps[gaps.size() / 2];
  }

  Pred predict(double r) const
  {
    if (anchor_rt.empty()) { return {r, scale(fallback_half_s), false}; }
    std::vector<double> res = gather(r);
    if (res.empty()) { return {r, scale(fallback_half_s), false}; }
    std::sort(res.begin(), res.end());
    // Guard against absurd shifts from unreliable anchors: never center the search
    // outside the donor's global RT window.
    const double shift = std::clamp(res[res.size() / 2], -hard_ceiling, hard_ceiling);
    if (res.size() < 2) { return {r + shift, scale(fallback_half_s), false}; }

    double mean = 0.0;
    for (double v : res) { mean += v; }
    mean /= static_cast<double>(res.size());
    double var = 0.0;
    for (double v : res) { var += (v - mean) * (v - mean); }
    var /= static_cast<double>(res.size() - 1);      // sample stddev (FlashLFQ)
    double sd = std::sqrt(var);
    if (! std::isfinite(sd)) { sd = 0.0; }            // NaN/inf guard -> floor window
    const double half = std::clamp(c * sd, floor_s, cap_s);
    return {r + shift, scale(half), true};
  }
};

namespace
{
// Build a LocalRtModel from the two runs' donors (peptides identified in both),
// keyed by sequence+charge; residual = acceptor-run RT - donor-run RT.
LocalRtModel build_local_rt_model(const DonorMap& donor_donors,
                                  const DonorMap& acceptor_donors)
{
  // One robust anchor per peptide (seq+charge) per run: median RT over duplicate
  // features, so multi-feature/charge-state peptides don't inject arbitrary
  // residuals that inflate the local-sigma tail driving the auto cap (Codex review).
  auto median_rt_by_key = [](const auto& storage) {
    std::unordered_map<std::string, std::vector<double>> tmp;
    for (const auto& d : storage) { tmp[feature_sequence_key(d->feature)].push_back(d->feature.getRT()); }
    std::unordered_map<std::string, double> med;
    med.reserve(tmp.size());
    for (auto& [k, v] : tmp) { std::sort(v.begin(), v.end()); med[k] = v[v.size() / 2]; }
    return med;
  };
  const auto donor_rt = median_rt_by_key(donor_donors.storage);
  const auto acc_rt = median_rt_by_key(acceptor_donors.storage);
  std::vector<std::pair<double, double>> pairs;  // (donor_rt, residual) one per shared key
  for (const auto& [key, drt] : donor_rt)
  {
    auto it = acc_rt.find(key);
    if (it == acc_rt.end()) { continue; }
    pairs.emplace_back(drt, it->second - drt);
  }
  std::sort(pairs.begin(), pairs.end());
  LocalRtModel m;  // params filled by set_local_rt_model
  m.anchor_rt.reserve(pairs.size());
  m.residual.reserve(pairs.size());
  for (const auto& p : pairs)
  {
    m.anchor_rt.push_back(p.first);
    m.residual.push_back(p.second);
  }
  return m;
}
} // namespace

/******************************************************************************/
Impl::Impl(const Param& params, const std::pair<double, double>& mz_range):
    mz_max_diff(params),
    mz_grid_center(0.5),
    rt_sec_max_window(params.getValue("distance_RT:max_difference")),
    mbr_fdr(params.getValue("fdr")),
    // setMinInt("min_decoys", 1) guarantees a positive value, so the cast
    // through int to size_t cannot underflow.
    min_decoys(static_cast<std::size_t>(
      static_cast<int>(params.getValue("min_decoys")))),
    // setMinInt("max_training_points", 0) guarantees a non-negative value.
    max_training_points(static_cast<std::size_t>(
      static_cast<int>(params.getValue("max_training_points")))),
    rng_seed_(static_cast<std::mt19937::result_type>(
      static_cast<int>(params.getValue("random_seed")))),
    rng_(rng_seed_)
{
  OPENMS_LOG_INFO << "MBR via PIP-ECHO(" << mz_range.first << ", "
                  << mz_range.second << ", " << mz_grid_center << ", "
                  << mz_max_diff.mz_diff(mz_range.second - mz_range.first)
                  << ", " << rt_sec_max_window << ", " << mbr_fdr << ")"
                  << std::endl;

  // Local adaptive RT window (experimental; default off -> the global window).
  use_local_rt_ = params.getValue("local_rt:enabled").toString() == "true";
  if (use_local_rt_)
  {
    lrt_cap_ = params.getValue("local_rt:max_window");
    lrt_floor_ = params.getValue("local_rt:min_window");
    lrt_locality_ = params.getValue("local_rt:anchor_window");
    lrt_n_anchors_ = static_cast<std::size_t>(static_cast<int>(params.getValue("local_rt:anchors")));
    lrt_c_ = params.getValue("local_rt:sigma_scale");
    lrt_fallback_half_ = params.getValue("local_rt:fallback_window");
    use_auto_ = params.getValue("local_rt:auto").toString() == "true";
    host_fwhm_ = params.getValue("local_rt:median_fwhm");
    OPENMS_LOG_INFO << "PIP-ECHO: local adaptive RT window enabled"
                    << (use_auto_ ? " (auto window estimation)" : "") << "." << std::endl;
  }
}

/******************************************************************************/
void Impl::set_local_rt_model(const DonorMap& donor_donors,
                              const DonorMap& acceptor_donors)
{
  if (! use_local_rt_) { return; }
  // FlashLFQ-inspired but generous (single-cell has poorer alignment + raw-RT-to-SVM
  // vs FlashLFQ's bulk + calibrated-RT-in-geomean): the window is data-driven c*sigma,
  // the cap is a generous backstop, and adaptive widening (group()) handles decoy
  // estimability. All values come from the 'local_rt:' parameters.
  LocalRtModel m = build_local_rt_model(donor_donors, acceptor_donors);
  m.locality_s = lrt_locality_;
  m.n_anchors = lrt_n_anchors_;
  m.c = lrt_c_;
  m.floor_s = lrt_floor_;
  m.cap_s = lrt_cap_;
  m.fallback_half_s = lrt_fallback_half_;
  m.widen = lrt_widen_;                                     // adaptive-widening multiplier
  m.hard_ceiling = rt_sec_max_window;                      // never exceed the global window
  if (m.floor_s > m.cap_s) { m.floor_s = m.cap_s; }         // guard std::clamp(lo>hi) UB
  local_rt_ = std::make_shared<const LocalRtModel>(std::move(m));
}

/******************************************************************************/
// [LOCAL-RT auto] Estimate the global (experiment-wide) window scales from the
// data so one config works across short and long gradients. Self-contained from
// feature RTs (local-sigma + anchor density); median_fwhm is an OPTIONAL host
// refinement (features carry no usable width). Sets lrt_cap_/floor_/locality_/
// fallback_; no-op unless 'local_rt:auto' is enabled.
void Impl::estimate_auto_params(const RunMap& runs)
{
  if (! use_local_rt_ || ! use_auto_) { return; }

  constexpr double MARGIN = 2.0;   // locality = gap_med * anchors * margin
  constexpr double K_FWHM = 4.0;   // cap FWHM lower-guard multiple
  const double c = lrt_c_;

  // Build each (donor_run, acceptor_run) pair model once; gather inter-anchor gaps.
  std::vector<LocalRtModel> models;
  std::vector<double> gaps;
  for (const auto& a_run : runs)
  {
    for (const auto& d_run : runs)
    {
      if (d_run.first == a_run.first) { continue; }
      LocalRtModel m = build_local_rt_model(d_run.second.donors, a_run.second.donors);
      m.n_anchors = lrt_n_anchors_;
      const double g = m.median_anchor_gap();
      if (g > 0.0) { gaps.push_back(g); }
      models.push_back(std::move(m));
    }
  }
  if (gaps.empty())
  {
    OPENMS_LOG_WARN << "[LOCAL-RT] auto: no shared anchors found; keeping default windows."
                    << std::endl;
    return;
  }
  std::sort(gaps.begin(), gaps.end());
  const double gap_med = gaps[gaps.size() / 2];
  const double locality = gap_med * static_cast<double>(lrt_n_anchors_) * MARGIN;

  // Sample the local RT-shift sigma within that locality across all pairs.
  std::vector<double> sigmas;
  for (auto& m : models)
  {
    m.locality_s = locality;
    for (double a : m.anchor_rt)
    {
      if (auto s = m.local_sigma_at(a); s.has_value() && *s > 0.0) { sigmas.push_back(*s); }
    }
  }

  const bool fwhm_used = host_fwhm_ > 0.0;
  double cap = 0.0, floor = 0.0;
  double p25 = 0, p50 = 0, p75 = 0, p90 = 0, p95 = 0;
  if (sigmas.size() >= 2)
  {
    std::sort(sigmas.begin(), sigmas.end());
    auto pct = [&](double p) {
      return sigmas[std::min(sigmas.size() - 1, static_cast<std::size_t>(p * sigmas.size()))];
    };
    p25 = pct(0.25); p50 = pct(0.50); p75 = pct(0.75); p90 = pct(0.90); p95 = pct(0.95);
    cap = std::max(c * p90, 0.25 * gap_med);         // sigma-driven, anchor-spacing guard
    // Floor = minimum window. With FWHM it should reflect chromatographic peak/apex
    // width, NOT alignment scatter (which already enters via c*sigma and the cap);
    // otherwise locally well-aligned donors get a needlessly broad minimum. Without
    // FWHM, fall back to a sigma-based floor (P25, not P10: avoid too-perfect tails).
    floor = fwhm_used ? std::max(1.5 * host_fwhm_, p25)
                      : std::max(c * p25, 0.10 * cap);
  }
  else  // too few sigma samples -> density-only fallback
  {
    cap = std::max(0.25 * gap_med, 1.0);
    floor = fwhm_used ? 1.5 * host_fwhm_ : 0.10 * cap;
  }
  if (fwhm_used) { cap = std::max(cap, K_FWHM * host_fwhm_); }  // physical lower guard on cap
  floor = std::min(floor, cap);

  lrt_cap_ = cap;
  lrt_floor_ = floor;
  lrt_locality_ = locality;
  lrt_fallback_half_ = floor;
  OPENMS_LOG_INFO << "[LOCAL-RT] auto windows (" << (fwhm_used ? "FWHM-aware" : "RT-residual only")
                  << "): cap=" << cap << "s floor=" << floor << "s locality=" << locality
                  << "s | local-sigma P25/50/75/90/95=" << p25 << "/" << p50 << "/" << p75 << "/"
                  << p90 << "/" << p95 << "s gap_med=" << gap_med << "s fwhm=" << host_fwhm_
                  << "s n=" << sigmas.size() << std::endl;
}

/******************************************************************************/
void Impl::partition_features(const std::vector<FeatureMap>& maps, RunMap& runs)
{
  std::size_t num_maps = maps.size();

  auto logger = ProgressLogger();
  logger.setLogType(ProgressLogger::CMD);
  logger.startProgress(0, num_maps, "building maps");

  for (std::size_t i = 0; i < num_maps; ++i)
  {
    logger.setProgress(i);

    const FeatureMap& map(maps[i]);
    std::string file_name(path_from_feature_map(map));
    Run& run = get_run_from_file_name(runs, file_name);

    for (auto& feature : map)
    {
      // We need to strip off the const-ness of the Feature to
      // prepare it for use.  We don't really care about it being
      // const, only that we have a reference/pointer.  But a
      // FeatureMap doesn't let us have a non-const reference
      // iterator.  So, we have to cheat and cast it.
      prepare_feature(const_cast<Feature&>(feature));
      run.insert(feature, i);
    }
  }

  logger.endProgress();
}

/******************************************************************************/
void Impl::link_donors_and_acceptors(const RunStatistics& stats,
                                     const DonorMap& donors,
                                     AcceptorMap& acceptors)
{
  for (auto& donor : donors.storage)
  {
    for (auto window = initial_window(*donor); window.has_value();
         window = next_window(window))
    {
      const auto target = find_acceptor_for(stats, acceptors, *donor, *window);

      if (target.has_value())
      {
        auto acceptor = target->second;
        acceptor->push_back(target->first, donor, DonorType::Target);
      }

      const auto random_donor = find_random_donor(donors, *donor, *window);
      match_t decoy = {}; // No std::optional::and_then in C++ 20 :(

      if (random_donor.has_value())
      {
        decoy = find_acceptor_for(stats, acceptors, *donor, *window,
                                  (*random_donor)->feature.getRT());

        if (decoy.has_value())
        {
          auto acceptor = decoy->second;
          acceptor->push_back(decoy->first, donor, DonorType::Decoy);
        }
      }

      if (target.has_value() || decoy.has_value())
      {
        break; // We can stop searching.
      }
    }
  }
}

/******************************************************************************/
Impl::match_t Impl::find_acceptor_for(const RunStatistics& stats,
                                      const AcceptorMap& acceptors,
                                      const Donor& donor,
                                      const Window& window,
                                      const std::optional<double> rt_override)
{
  // FIXME: What does the paper say about finding acceptor peaks?
  // In FlashLFQ they do some weird envelope cutting and don't use
  // the actual found peak (FindIndividualAcceptorPeak).
  const double anchor_rt = rt_override.value_or(donor.feature.getRT());

  // [LOCAL-RT] re-center the search on the locally predicted acceptor RT
  // and gate to the local residual scatter. Unlike the global baseline, the window
  // is ALWAYS bounded (<= cap_s; tight fallback when local anchors are sparse) --
  // there is no fall-through to the +/-1206 s global window. The score's
  // rt_diff_error is measured against the predicted center. Baseline path (env
  // unset) is byte-identical: local_on stays false, no re-center, no gate.
  double center_rt = anchor_rt;
  double rt_gate = std::numeric_limits<double>::infinity();  // no extra gate (baseline)
  bool local_on = false;
  if (use_local_rt_ && local_rt_)
  {
    const LocalRtModel::Pred p = local_rt_->predict(anchor_rt);
    center_rt = p.center;
    rt_gate = p.half_window;
    local_on = true;
    if (p.supported) { ++lrt_supported_; } else { ++lrt_fallback_; }
  }
  const std::optional<double> eff_override =
    local_on ? std::optional<double>(center_rt) : rt_override;

  const AcceptorMap::grid_center_t center(center_rt, donor.feature.getMZ());
  const AcceptorMap::grid_index_t index
    = acceptors.grid.cellIndexAtClusterCenter(center);

  return acceptors.nearby<match_t>(
    index, window.grid_neighbors, match_t {},
    [&](match_t best_match, AcceptorMap::value_type acceptor) -> match_t {
      if (acceptor->is_donor_compatible(donor, window, eff_override)
          && std::fabs(acceptor->feature.getRT() - center_rt) <= rt_gate)
      {
        Score score(stats.score(donor.feature, acceptor->feature, eff_override));

        if (! best_match.has_value()
            || best_match->first.mbr_score < score.mbr_score)
        {
          return std::make_pair(score, acceptor);
        }
      }

      return best_match;
    });
}

/******************************************************************************/
std::optional<std::shared_ptr<const Donor>>
Impl::find_random_donor(const DonorMap& donors,
                        const Donor& start,
                        const Window& window) const
{
  /// FIXME: Is this mass calculation good enough?
  const auto mass = [](const Donor& donor) -> double {
    return (donor.feature.getMZ() - Constants::PROTON_MASS_U)
           * donor.feature.getCharge();
  };

  /// FIXME: Is this sequence calculation good enough?
  /// FIXME: What is a "Base Sequence"?
  const auto base_seq = [](const Donor& donor) -> std::string {
    const auto hit(PipEcho::Util::feature_hit(donor.feature));
    if (! hit.has_value())
    {
      throw Exception::MissingInformation(__FILE__, __LINE__,
                                          OPENMS_PRETTY_FUNCTION,
                                          "donor feature missing peptide sequence");
    }
    return hit->getSequence().toUnmodifiedString();
  };

  const double mass_min_diff = 5.0 * Constants::PROTON_MASS_U;
  const double mass_max_diff = 11.0 * Constants::PROTON_MASS_U;
  // [LOCAL-RT] with tight local windows, requiring the decoy anchor
  // > 2*GLOBAL window away starves decoys; require only > ~2 local windows
  // (>= 120 s) -- far enough to avoid the same elution event, near enough to
  // keep decoys plentiful for FDR estimation.
  double rt_min_diff = window.rt_tol * 2.0; // FIXME: Is this okay?
  if (use_local_rt_ && local_rt_)
  {
    const LocalRtModel::Pred p = local_rt_->predict(start.feature.getRT());
    // Use the BASE (un-widened) local window for decoy-donor separation: tying it
    // to the widened window would push decoy donors farther away on each widening
    // pass, opposing the decoy generation widening is meant to achieve (Codex review).
    const double lw = p.half_window / std::max(local_rt_->widen, 1.0);
    rt_min_diff = std::max(2.0 * lw, 120.0);
  }

  const double start_mass = mass(start);
  const double start_rt = start.feature.getRT();
  const std::string start_seq = base_seq(start);

  std::vector<std::shared_ptr<const Donor>> matching_donors;

  auto explore
    = [&](double distance, std::shared_ptr<const Donor> donor) -> double {
    const double mass_diff = std::fabs(start_mass - mass(*donor));
    const double rt_diff = std::fabs(start_rt - donor->feature.getRT());

    if (rt_diff > rt_min_diff && mass_diff > mass_min_diff
        && mass_diff < mass_max_diff && start_seq != base_seq(*donor))
    {
      matching_donors.push_back(donor);
    }

    return std::max(distance, mass_diff);
  };

  double cur_max_mass_diff = mass_max_diff;
  Int64 neighbors = window.grid_neighbors;

  const DonorMap::grid_center_t center(start.feature.getRT(),
                                       start.feature.getMZ());
  const DonorMap::grid_index_t index
    = donors.grid.cellIndexAtClusterCenter(center);

  do
  {
    const double distance = donors.nearby<double>(index, neighbors, 0, explore);

    // Adjust search params for next round.
    if (distance < cur_max_mass_diff)
    {
      // We're not exploring far enough.
      ++neighbors;
    }
    else
    {
      // We need to increase our tolerance.
      cur_max_mass_diff *= 10;
    }
  } while (matching_donors.empty() && cur_max_mass_diff < 1e5 && neighbors < 5);

  if (! matching_donors.empty())
  {
    std::uniform_int_distribution<std::size_t> pick(0,
                                                    matching_donors.size() - 1);
    return matching_donors[pick(rng_)];
  }

  return {};
}

/******************************************************************************/
void Impl::generate_consensus_map(RunMap& runs, ConsensusMap& consensus_map)
{
  // [LOCAL-RT] report how often the local window was used vs fell back.
  if (use_local_rt_)
  {
    const std::size_t tot = lrt_supported_ + lrt_fallback_;
    const double frac = tot ? (100.0 * lrt_fallback_ / tot) : 0.0;
    OPENMS_LOG_INFO << "[LOCAL-RT] find_acceptor_for queries: supported="
                    << lrt_supported_ << " fallback=" << lrt_fallback_ << " ("
                    << frac << "% fallback)" << std::endl;
  }

  using ident_val_t = std::set<FeatureRef, FeatureRefCmp>;
  using ident_map_t = std::map<std::string, ident_val_t>;
  struct transfer_info_t
  {
    IntList map_indices;
    DoubleList q_values;
  };
  using transfer_map_t = std::map<std::string, transfer_info_t>;

  std::size_t acceptors_seen {}, decoys_seen {};
  ident_map_t ident_map;
  transfer_map_t transfer_map;

  auto insert = [&](const std::string key, const FeatureRef& peak) {
    auto [place, inserted]
      = ident_map.insert(std::make_pair(key, ident_val_t {peak}));
    if (! inserted) place->second.insert(peak);
  };

  for (auto& run : runs)
  {
    // Place each donor into the identity map.
    for (auto& donor : run.second.donors.storage)
    {
      std::string key = feature_sequence_key(donor->feature);
      insert(key, *donor);
    }
  }

  auto view = runs | std::views::transform([](auto& run) {
                return run.second.acceptors.storage;
              })
              | std::views::join;

  std::vector<std::shared_ptr<Acceptor>> all_acceptors;

  //  FIXME
  // all_acceptors.reserve(std::ranges::size(view));
  // all_acceptors.insert(std::ranges::begin(view),
  // std::ranges::end(view));

  for (auto& acpt : view)
  {
    all_acceptors.push_back(acpt);
  }

  PipEcho::TransferFDRModel transfer_fdr(all_acceptors, min_decoys, max_training_points);
  const PipEcho::TransferFDRModel::group_t& acceptors = transfer_fdr.run(mbr_fdr);

  // Put each acceptor into the correct feature bucket.
  for (auto& acceptor : acceptors)
  {
    // Decoy transfers exist only to estimate the MBR FDR; they must not leak
    // into the quantitative consensus output.
    if (! acceptor->is_target())
    {
      ++decoys_seen;
      continue;
    }

    ++acceptors_seen;

    std::string key
      = feature_sequence_key(acceptor->donor_ident, acceptor->donor_charge);

    transfer_info_t& transfer_info = transfer_map[key];
    transfer_info.map_indices.push_back(
      static_cast<Int>(acceptor->acceptor.get().map_index));
    transfer_info.q_values.push_back(acceptor->q_value);

    insert(key, acceptor->acceptor);
  }

  OPENMS_LOG_INFO << "PIP-ECHO: added " << acceptors_seen
                  << " acceptors to the consensus map with " << decoys_seen
                  << " decoys seen" << std::endl;

  // We can now turn that IdentMap into a ConsensusMap.
  for (auto& group : ident_map)
  {
    ConsensusFeature consensus_feature;

    for (auto& peak : group.second)
    {
      // Note: This is apparently the correct way to add a feature
      // to a ConsensusFeature.  It does more work, such as copying
      // identifications, that the other insert methods don't do.
      consensus_feature.insert(peak.map_index, peak.feature);
    }

    auto transfer = transfer_map.find(group.first);
    if (transfer != transfer_map.end())
    {
      consensus_feature.setMetaValue("mbr_transfer_map_indices",
                                     transfer->second.map_indices);
      consensus_feature.setMetaValue("mbr_transfer_qvalues",
                                     transfer->second.q_values);
    }

    consensus_feature.computeConsensus();
    consensus_map.push_back(consensus_feature);
  }
}

/******************************************************************************/
std::string Impl::path_from_feature_map(const FeatureMap& map)
{
  StringList paths;
  map.getPrimaryMSRunPath(paths);

  if (paths.empty())
  {
    OPENMS_LOG_WARN << "Warning: feature map has no MS run path!" << std::endl;
    return "unknown";
  }
  else
  {
    if (paths.size() > 1)
    {
      OPENMS_LOG_WARN << "Warning: feature map has more than one MS run path, "
                      << "using first path: " << paths[0] << std::endl;
    }

    return paths[0];
  }
}

/******************************************************************************/
Run& Impl::get_run_from_file_name(RunMap& runs, const std::string& file_name)
{
  auto run_it = runs.find(file_name);

  if (run_it == runs.end())
  {
    const auto [it, status]
      = runs.try_emplace(file_name, rt_sec_max_window, mz_grid_center);

    assert(status);
    run_it = it;
  }

  return run_it->second;
  ;
}

/******************************************************************************/
std::optional<Window> Impl::initial_window(const Donor& donor)
{
  double mz_diff = mz_max_diff.mz_diff(donor.feature.getMZ());

  return Window {rt_sec_max_window, mz_diff, 1};
}

/******************************************************************************/
std::optional<Window> Impl::next_window(const std::optional<Window>& prev)
{
  if (prev && prev->rt_tol <= rt_sec_max_window && prev->grid_neighbors < 3)
  {
    Window next(*prev);
    next.grid_neighbors++;
    return next;
  }
  else
  {
    return {};
  }
}

} // namespace OpenMS::PipEcho
