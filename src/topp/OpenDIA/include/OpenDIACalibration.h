// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIARTMap.h>
#include <OpenDIAPeakPicker.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenDIALibrary.h>
#include <OpenDIAMzIndex.h>

#include <array>
#include <cstddef>
#include <functional>
#include <iosfwd>
#include <limits>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS::OpenDIA
{
  /**
    @brief Stage-1 map and windows of a calibration whose stage 3 had too few anchors.

    Used for the windowed second pass: a larger panel is collected only inside these
    windows around the stage-1 map, and stage 3 and 4 are repeated on it unchanged.
  */
  struct CalibrationWindow
  {
    RTMap map; ///< Stage-1 library RT to run RT map.
    double half_window_s = 0.0; ///< Stage-3 RT half window.
    std::size_t half_window_spectra = 0; ///< Co-elution half window in spectra.
    std::size_t stage1_kept = 0;
    double stage1_sigma_s = 0.0;
    double peak_fwhm_s = 0.0;
    /// Library RT range of the stage-1 anchors behind @p map. Outside it the map only extrapolates, so the windowed
    /// pass collects and searches those members from their window to the run edge instead of around the map.
    double support_low = -std::numeric_limits<double>::infinity();
    double support_high = std::numeric_limits<double>::infinity();
    /// Measured in the first pass (memory sizing only): median spectrum spacing per precursor and the MS2 RT range.
    double cycle_s = 0.0;
    double run_first_rt = 0.0, run_last_rt = 0.0;
  };

  /** @brief Fixed settings of the identification-free run calibration. */
  struct CalibrationOptions
  {
    std::size_t panel_targets = 10000; ///< Targets sampled evenly over library RT; each brings its DECOY_ partner.
    double search_ppm = 15.0; ///< Uncalibrated fragment search half width.
    double search_im = 0.05; ///< Uncalibrated 1/K0 search half width around the library mobility.
    std::size_t candidates = 10; ///< Local maxima of the summed trace evaluated per precursor.
    std::size_t seed_half_window = 4; ///< Co-elution half window (spectra) before the peak width is measured.
    // Stage 1: whole-run anchors by paired target-decoy competition (no RT information).
    std::size_t strict_min_fragments = 5; ///< Library fragments a panel precursor needs.
    double stage1_fdr = 0.05; ///< Paired target-decoy q-value cut for the whole-run anchors.
    std::size_t mode_bins = 20;
    double mode_bandwidth_s = 15.0;
    // Stage 2: identification-free latent peak/no-peak mixture.
    double presence_prior = 0.3; ///< Fixed presence prior pi (v0.13).
    double em_sigma_fraction = 0.07; ///< Fixed Gaussian RT width as a fraction of the run duration (v0.13).
    double score_beta = 0.8; ///< Weight of the clipped within-precursor score contrast.
    std::size_t em_knots = 11;
    std::size_t em_max_iterations = 200;
    // Stage 3: target-decoy controlled anchors inside the stage-2 window.
    std::size_t min_fragments = 4;
    double min_window_s = 30.0;
    double anchor_fdr = 0.01;
    // Tolerance multiples of the robust error widths.
    double rt_sigma_multiple = 4.0;
    double ms2_sigma_multiple = 5.0; ///< Unused since the mass tolerance is fitted (MassToleranceFit); kept for the JSON options record.
    double im_sigma_multiple = 4.0;
  };

  /** @brief One confident anchor: a panel target accepted at the anchor FDR. */
  struct CalibrationAnchor
  {
    std::string precursor_id;
    double library_rt = 0.0;
    double run_rt = 0.0;
    double library_im = 0.0;
    double observed_im = 0.0; ///< Median 1/K0 of the apex fragment centroids; 0 without mobility data.
    double correlation = 0.0;
    double cosine = 0.0;
    double q_value = 1.0;
    std::size_t fragments = 0;
    bool used = false; ///< Survived RT outlier rejection.
    // Raw-scan fragment traces around the anchor apex, input of the peak-boundary cap estimate.
    std::size_t anchor_scan = 0; ///< Index of the calibration apex in trace_rt.
    std::vector<double> trace_rt; ///< Spectrum RTs of the stored trace slice (at most +-64 spectra around the apex).
    std::vector<std::array<float, 6>> trace; ///< Fragment intensities per stored spectrum (unused fragments are 0).
  };

  /** @brief Data-driven fragment mass tolerance: signal (two-Gaussian scale mixture) plus uniform background, in ppm. */
  struct MassToleranceFit
  {
    double centre = 0.0;                  ///< Common centre of the signal components (the mass offset).
    std::array<double, 2> weights{};      ///< Shares of the narrow and the wide signal component.
    std::array<double, 2> sigmas{};       ///< Widths of the narrow and the wide signal component.
    double background_share = 0.0;        ///< Share of errors that are uniform background over the search window.
    double half_width = 0.0;              ///< Distance from the centre where signal density equals background density.
  };

  /** @brief Calibrated library-to-run coordinates and extraction tolerances. */
  struct CalibrationResult
  {
    // Panel and stage counts.
    std::size_t panel_targets = 0, panel_decoys = 0, pairs = 0;
    std::size_t stage1_targets = 0, stage1_decoys = 0, stage1_kept = 0;
    double stage1_sigma_s = 0.0;
    std::size_t em_iterations = 0;
    bool em_converged = false;
    double em_sigma_s = 0.0; ///< Fixed model width in seconds (em_sigma_fraction x run duration).
    double em_presence = 0.0; ///< Mean peak responsibility per hypothesis at convergence (diagnostic).
    std::vector<double> em_knots_library, em_knots_run;
    double stage3_half_window_s = 0.0;
    std::size_t stage3_extended_targets = 0; ///< Windowed pass: targets searched to the run edge (outside the stage-1 support).
    bool windowed = false; ///< Windowed second pass: stages 1 and 2 were skipped, stage 3 centred on the stage-1 map.
    std::size_t half_window_spectra = 0;
    double peak_fwhm_s = 0.0;
    std::size_t anchor_targets = 0, anchor_decoys = 0;
    // Final models.
    RTMap rt_map;
    double rt_p50_s = 0.0, rt_p95_s = 0.0, rt_p99_s = 0.0;
    /// RT error of an arbitrary precursor: the Gaussian width of a signal-plus-background fit to all stage-3 target
    /// winners' residuals to rt_map (background shape from the decoy winners). Not selected for small residuals, unlike
    /// the anchors behind rt_map.sigma. Equals rt_map.sigma when the fit is not possible.
    double rt_spread_s = 0.0;
    double rt_spread_offset_s = 0.0; ///< Offset of that Gaussian from rt_map (0 when not fitted).
    double rt_spread_share = 0.0; ///< Share of the target winners in the Gaussian (0 when not fitted).
    std::size_t rt_spread_points = 0; ///< Target winners in the fit (0 when not fitted).
    double rt_tolerance_s = 0.0; ///< Half width, rt_sigma_multiple x rt_spread_s.
    double ms2_offset_ppm = 0.0, ms2_sigma_ppm = 0.0, ms2_tolerance_ppm = 0.0;
    MassToleranceFit ms2_fit; ///< Model behind ms2_tolerance_ppm.
    std::size_t ms2_points = 0;
    BoundaryHalfWidthFit boundary_fit; ///< Peak boundary cap in half widths, estimated on the used anchors.
    bool has_im = false;
    double im_intercept = 0.0, im_slope = 1.0; ///< observed 1/K0 = intercept + slope x library 1/K0.
    double im_sigma = 0.0, im_tolerance = 0.0;
    std::size_t im_points = 0;
    std::vector<CalibrationAnchor> anchors;

    /// Library RT to run RT in seconds.
    double runRT(double library_rt) const { return rt_map(library_rt); }
    /// Library 1/K0 to calibrated run 1/K0.
    double runIM(double library_im) const { return im_intercept + im_slope * library_im; }
  };

  /**
    @brief Identification-free calibration of RT, fragment m/z and 1/K0 from the run itself.

    A panel of library targets, stratified over library RT, and their paired decoys
    is extracted across the whole run in one streaming pass: per MS2 spectrum that
    isolates the precursor, the most intense centroid within the wide m/z and 1/K0
    search windows of each of the top six library fragments. No identifications,
    iRT peptides or external RTs are used.

    Stage 1 (DIATraceExtractor concept): strict co-elution anchors over the whole
    run (fragment Pearson correlation x library cosine), filtered by the per-bin RT
    mode, give the initial map. Stage 2 (identification-free concept): a latent
    peak-versus-no-peak mixture over all candidate peaks of all panel targets,
    weighted by within-precursor score contrast and a Gaussian RT term around the
    current monotone map, refines the map and estimates its RT spread. Stage 3:
    targets and their paired decoys are rescored inside max(4 sigma, min_window)
    of that map; paired target-decoy competition selects anchors at anchor_fdr,
    ranking each member by the evidence of its co-elution score (-log of its
    upper-tail p-value among the panel's decoy scores) plus the log-likelihood
    ratio of its apex RT residual under a Gaussian of the stage-1 spread against a
    uniform over the window. The score alone stops separating targets from decoys
    when decoys match background as often as targets do (predicted decoys); a
    chance peak's RT is uniform over the window, a target's is near the map.
    Stage 4: the final RT map, fragment mass offset and 1/K0 model are fitted on
    the anchors, with tolerances as fixed multiples of their robust widths.

    The collector keeps the panel's traces in memory (panel x spectra x 6 x 12 bytes).

    @note This is a heuristic calibration, not an identification
    or FDR-controlled result. Anchor q-values only control the anchor set.
  */
  class CalibrationCollector final : public Interfaces::IMSDataConsumer
  {
  public:
    static constexpr std::size_t FRAGMENTS = 6;

    /**
      @brief Select the panel from @p library.
      @param[in] library Prepared library with library RT and DECOY_ naming available through @p traml_ids.
      @param[in] traml_ids PQP precursor ID to TRAML_ID map, used to pair decoys with targets.
      @param[in] options Fixed calibration settings.
      @param[in] window If set, collect each member only inside its stage-3 window around this stage-1 map
                        (plus a 30 s margin), and let calibrate() skip stages 1 and 2.
      @throws std::invalid_argument If no eligible target-decoy pairs exist or options are invalid.
    */
    CalibrationCollector(const OpenDIA::Library& library,
                         const std::vector<std::pair<std::string, std::string>>& traml_ids,
                         CalibrationOptions options = {}, std::optional<CalibrationWindow> window = std::nullopt);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override {}
    void setExpectedSize(std::size_t, std::size_t) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}

    /// Treat a spectrum ion-mobility array with a non-CCS unit as 1/K0 (VSSC).
    void setIonMobilityArrayAsVSSC(bool value) { ion_mobility_array_as_vssc_ = value; }
    /// Worker threads for the per-spectrum trace collection (results do not depend on it). Default 1.
    void setThreads(int threads) { threads_ = threads < 1 ? 1 : threads; }

    /**
      @brief Run stages 1-4 on the collected traces.
      @param[out] retry_window Optional: set when stage 1 succeeded but stage 3 found too few anchors (seed for a
                  windowed second pass); left unchanged otherwise.
      @throws std::runtime_error If a stage has too few anchors to fit.
    */
    CalibrationResult calibrate(std::optional<CalibrationWindow>* retry_window = nullptr) const;

    /// Number of panel precursors (targets plus decoys).
    std::size_t panelSize() const { return panel_.size(); }

  private:
    struct Member
    {
      std::string id;
      double precursor_mz = 0.0;
      double library_rt = 0.0;
      double library_im = 0.0;
      bool is_decoy = false;
      std::size_t partner = 0; ///< Index of the paired decoy/target.
      std::size_t fragment_count = 0;
      double window_low = -std::numeric_limits<double>::infinity(); ///< Collected run RT range (windowed pass only).
      double window_high = std::numeric_limits<double>::infinity();
      std::array<double, FRAGMENTS> fragment_mz{};
      std::array<double, FRAGMENTS> library_intensity{};
      // One entry per isolating spectrum, in RT order. im is empty when the member has no library 1/K0.
      std::vector<double> rt;
      std::vector<std::array<float, FRAGMENTS>> intensity, ppm, im; ///< ppm: error of the strongest centroid.
      bool closed = false; ///< Past window_high: traces complete and trimmed.
    };

    struct Candidate
    {
      std::size_t member = 0;
      std::size_t apex = 0; ///< Index into the member's spectra.
      double rt = 0.0;
      double correlation = 0.0;
      double cosine = 0.0;
      double score = -1.0;
      std::size_t fragments = 0;
    };

    std::vector<Candidate> candidates_(std::size_t member, std::size_t half_window, double rt_low, double rt_high,
                                       std::size_t max_candidates, std::size_t min_fragments) const;

    /// Stages 3 and 4 around @p centre (library RT to run RT); throws if too few anchors. Targets whose library RT
    /// is below @p support_low (above @p support_high) and their decoys are searched from the run start (to the run end).
    void finish_(CalibrationResult& result, const std::function<double(double)>& centre, std::size_t half_window,
                 double support_low = -std::numeric_limits<double>::infinity(),
                 double support_high = std::numeric_limits<double>::infinity()) const;

    CalibrationOptions options_;
    std::optional<CalibrationWindow> window_;
    std::vector<Member> panel_;
    std::vector<std::size_t> by_mz_;
    bool ion_mobility_array_as_vssc_ = false;
    bool saw_im_ = false;
    int threads_ = 1;
    SpectrumMzIndex mz_index_;
    double first_rt_ = 0.0, last_rt_ = 0.0;
    bool saw_spectrum_ = false;
  };

  /// Write @p result as JSON (units in keys) to @p out.
  void writeCalibrationJSON(const CalibrationResult& result, const CalibrationOptions& options, std::ostream& out);

  /** @brief Models read back from a calibration JSON written by writeCalibrationJSON(). */
  struct CalibrationRT
  {
    RTMap map; ///< Library RT to run RT (s).
    double tolerance_s = 0.0; ///< Half width of the RT window.
    double ms2_offset_ppm = 0.0; ///< Run fragment m/z = library m/z x (1 + offset x 1e-6).
    bool has_ms2 = false; ///< False when too few fragment errors were measured.
    double ms2_tolerance_ppm = 0.0; ///< Half width around the calibrated fragment m/z.
    bool has_im = false; ///< False without mobility data or too few anchors.
    double im_intercept = 0.0, im_slope = 1.0; ///< Run 1/K0 = intercept + slope x library 1/K0.
    double im_tolerance = 0.0; ///< Half width around the calibrated 1/K0.
    double boundary_half_widths = 0.0; ///< Peak boundary cap per side in half widths; 0 when absent (older files).
  };

  /**
    @brief Read the RT map and tolerance from a calibration JSON written by writeCalibrationJSON().
    @param[in] name Used in error messages.
    @throws std::invalid_argument If the JSON is not version 1, or the map or tolerance is invalid.
  */
  CalibrationRT readCalibrationRT(std::istream& in, const std::string& name);
}
