// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include <OpenDIARTMap.h>

#include <OpenMS/INTERFACES/IMSDataConsumer.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

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

    Used for an experimental windowed second pass: a larger panel is collected only inside these
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
    double ms2_sigma_multiple = 5.0;
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
    bool windowed = false; ///< Experimental windowed second pass: stages 1 and 2 were skipped, stage 3 centred on the stage-1 map.
    std::size_t half_window_spectra = 0;
    double peak_fwhm_s = 0.0;
    std::size_t anchor_targets = 0, anchor_decoys = 0;
    // Final models.
    RTMap rt_map;
    double rt_p50_s = 0.0, rt_p95_s = 0.0, rt_p99_s = 0.0;
    double rt_tolerance_s = 0.0; ///< Half width, rt_sigma_multiple x rt_map.sigma.
    double ms2_offset_ppm = 0.0, ms2_sigma_ppm = 0.0, ms2_tolerance_ppm = 0.0;
    std::size_t ms2_points = 0;
    bool has_im = false;
    double im_intercept = 0.0, im_slope = 1.0; ///< observed 1/K0 = intercept + slope x library 1/K0.
    double im_sigma = 0.0, im_tolerance = 0.0;
    std::size_t im_points = 0;
    std::vector<CalibrationAnchor> anchors;

    /// Library RT to run RT in seconds.
    double runRT(double library_rt) const { return rt_map(library_rt); }
    /// Library 1/K0 to calibrated run 1/K0.
    double runIM(double library_im) const { return im_intercept + im_slope * library_im; }
    /// Library m/z to calibrated MS2 m/z.
    double runMZ(double library_mz) const { return library_mz * (1.0 + ms2_offset_ppm * 1e-6); }
  };

  /**
    @brief Per panel member stage values, for diagnosing a calibration (also filled when it fails).

    Rows are in panel order. Values of a stage that was not reached stay NaN (q-values 1).
  */
  struct CalibrationDiagnostics
  {
    struct Row
    {
      std::string precursor_id;
      bool is_decoy = false;
      double library_rt = 0.0;
      std::size_t spectra = 0; ///< Isolating MS2 spectra collected.
      double stage1_rt = std::numeric_limits<double>::quiet_NaN(), stage1_score = std::numeric_limits<double>::quiet_NaN();
      double stage1_q = 1.0;
      double stage1_map_rt = std::numeric_limits<double>::quiet_NaN(); ///< Stage-1 map at the target library RT.
      double em_rt = std::numeric_limits<double>::quiet_NaN(); ///< Stage-2 map at the target library RT.
      double stage3_centre_rt = std::numeric_limits<double>::quiet_NaN(); ///< Centre of the stage-3 window.
      double stage3_rt = std::numeric_limits<double>::quiet_NaN(), stage3_score = std::numeric_limits<double>::quiet_NaN();
      double stage3_q = 1.0;
    };
    std::vector<Row> rows;
    double first_rt = 0.0, last_rt = 0.0;
    double cycle_s = 0.0; ///< Median spectrum spacing per precursor.
    double peak_fwhm_s = 0.0; ///< Median stage-1 anchor FWHM (0 if stage 1 failed).
    std::size_t half_window_spectra = 0; ///< Stage-2/3 co-elution half window (0 if stage 1 failed).
    std::string failure; ///< Empty if calibration succeeded.
    /// Set when stage 1 succeeded but stage 3 found too few anchors (seed for a windowed second pass).
    std::optional<CalibrationWindow> window;
  };

  /// Write @p diagnostics as TSV to @p out (header comment carries the run RT range and any failure).
  void writeCalibrationDiagnostics(const CalibrationDiagnostics& diagnostics, std::ostream& out);

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
    of that map; paired target-decoy competition selects anchors at anchor_fdr.
    Stage 4: the final RT map, fragment mass offset and 1/K0 model are fitted on
    the anchors, with tolerances as fixed multiples of their robust widths.

    The collector keeps the panel's traces in memory (panel x spectra x 6 x 12 bytes).

    @note This is an experimental heuristic calibration, not an identification
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
                        (plus a 30 s margin), and let calibrate() skip stages 1 and 2. Experimental.
      @throws std::invalid_argument If no eligible target-decoy pairs exist or options are invalid.
    */
    CalibrationCollector(const OpenSwath::LightTargetedExperiment& library,
                         const std::vector<std::pair<std::string, std::string>>& traml_ids,
                         CalibrationOptions options = {}, std::optional<CalibrationWindow> window = std::nullopt);

    void consumeSpectrum(MSSpectrum& spectrum) override;
    void consumeChromatogram(MSChromatogram&) override {}
    void setExpectedSize(std::size_t, std::size_t) override {}
    void setExperimentalSettings(const ExperimentalSettings&) override {}

    /// Treat a spectrum ion-mobility array with a non-CCS unit as 1/K0 (VSSC).
    void setIonMobilityArrayAsVSSC(bool value) { ion_mobility_array_as_vssc_ = value; }

    /**
      @brief Run stages 1-4 on the collected traces.
      @param[out] diagnostics Optional per-member stage values; filled up to the stage reached, also on failure.
      @throws std::runtime_error If a stage has too few anchors to fit.
    */
    CalibrationResult calibrate(CalibrationDiagnostics* diagnostics = nullptr) const;

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
      // One entry per isolating spectrum, in RT order.
      std::vector<double> rt;
      std::vector<std::array<float, FRAGMENTS>> intensity, ppm, im;
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

    /// Stages 3 and 4 around @p centre (library RT to run RT); throws if too few anchors.
    void finish_(CalibrationResult& result, const std::function<double(double)>& centre, std::size_t half_window,
                 CalibrationDiagnostics* diagnostics) const;

    CalibrationOptions options_;
    std::optional<CalibrationWindow> window_;
    std::vector<Member> panel_;
    std::vector<std::size_t> by_mz_;
    bool ion_mobility_array_as_vssc_ = false;
    bool saw_im_ = false;
    double first_rt_ = 0.0, last_rt_ = 0.0;
    bool saw_spectrum_ = false;
  };

  /// Write @p result as JSON (units in keys) to @p out.
  void writeCalibrationJSON(const CalibrationResult& result, const CalibrationOptions& options, std::ostream& out);

  /// Write anchors as TSV to @p out.
  void writeCalibrationAnchors(const CalibrationResult& result, std::ostream& out);

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
  };

  /**
    @brief Read the RT map and tolerance from a calibration JSON file.
    @throws std::invalid_argument If the file is missing, not version 1, or the map or tolerance is invalid.
  */
  CalibrationRT readCalibrationRT(const std::string& path);
}
