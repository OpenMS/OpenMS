// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

#include "OpenDIATraces.h"

#include <array>
#include <cstddef>
#include <string>
#include <unordered_map>
#include <vector>

namespace OpenMS::OpenDIA
{
  /// Feature order of the experimental second-stage score. Feature 0 is the main (initial ranking) score.
  enum class Stage2Feature : std::size_t
  {
    CNN_Z = 0,             ///< CNN logit standardised by its fold's decoy precursor winners (median, 1.4826 MAD).
    FRAGMENT_FRACTION,     ///< Share of library fragments with signal in the core window.
    FRAGMENT_COELUTION,    ///< Mean pairwise Pearson correlation of fragment traces in the core window.
    LIBRARY_CORRELATION,   ///< Pearson correlation of sqrt core fragment areas and sqrt library intensities.
    MS2_ABS_PPM,           ///< |intensity-weighted fragment ppm error - run offset|, capped.
    MS1_ABS_PPM,           ///< |intensity-weighted MS1 monoisotope ppm error - run offset|, capped.
    MS1_COELUTION,         ///< Pearson correlation of the MS1 monoisotope trace and the summed fragment trace.
    MS1_ISOTOPE,           ///< 1 if both the MS1 monoisotope and first isotope carry signal in the core window.
    COUNT
  };
  constexpr std::size_t STAGE2_FEATURES = static_cast<std::size_t>(Stage2Feature::COUNT);
  /// Column names, in Stage2Feature order.
  const std::array<std::string, STAGE2_FEATURES>& stage2FeatureNames();

  struct Stage2ExtractionOptions
  {
    double rt_half_width = 15.0;     ///< Tensor RT half width in seconds (TraceCollector).
    double core_half_width = 4.0;    ///< Half width around the apex used by every feature.
    double ms1_offset_ppm = 0.0;     ///< Run MS1 offset subtracted before taking the absolute error.
    double ms2_offset_ppm = 0.0;     ///< Run fragment offset subtracted before taking the absolute error.
    double ppm_cap = 20.0;           ///< Absolute errors are capped here; missing signal also gets the cap.
    /// Experimental refined peak picking for the boundary group: 0 = off; otherwise fragments with refinement r below this
    /// value (against the median unit-area shape of the first pick) are left out when the peak group is picked again.
    double refined_picking_min_r = 0.0;
  };

  /**
    @brief Extraction features 1..7 of one window (feature 0, the CNN score, is filled by the caller).

    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats, fragment-major as in the collector.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence from TraceCollector::enableEvidence().
    @param[in] options Window geometry and run mass offsets.
    @return All features with feature 0 set to zero. Undefined correlations are zero.
  */
  std::array<double, STAGE2_FEATURES> stage2Features(const float* tensor, const std::array<double, 6>& library_intensities,
                                                     const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                     const Stage2ExtractionOptions& options);

  /**
    @brief Which feature groups enter the second-stage discriminant. Feature 0 is always the standardised CNN logit.

    Row order: cnn_z, then the seven extraction features, then the search features, then the isotope features, then the
    boundary features, then the b/y ion counts, then the shape features, then the refine features, each group only if selected.
  */
  struct Stage2FeatureSet
  {
    bool extraction = true;   ///< Stage2Feature columns 1..7 (tensor and MS1 evidence; the original -stage2 set).
    bool search = false;      ///< The four search peak-group features (needs the OpenDIA -features_out table).
    bool isotope = false;     ///< The MS1 isotope-envelope features (Stage2IsotopeFeature).
    bool boundary = false;    ///< Experimental peak-group boundary features (Stage2BoundaryFeature).
    bool byion = false;       ///< Experimental apex-spectrum b/y ion counts (Stage2ByIonFeature; needs TraceCollector::enableByIonCounts()).
    bool shape = false;       ///< Experimental lagged cross-correlation, library agreement and signal features over the picked peak (Stage2ShapeFeature).
    bool refine = false;      ///< Experimental interference-robust transition refinement features (Stage2RefineFeature).
  };
  /**
    @brief Parse a '+'-separated list of groups, e.g. "extraction", "search+isotope"; "extraction_search" means "extraction+search".
    @throws std::invalid_argument On an unknown or repeated group name.
  */
  Stage2FeatureSet parseStage2FeatureSet(const std::string& name);
  /// Column names of a feature set, in row order.
  std::vector<std::string> stage2FeatureSetNames(const Stage2FeatureSet& set);

  /// MS1 isotope-envelope features, following the OpenSWATH MS1 isotope scores with three extra isotopes.
  enum class Stage2IsotopeFeature : std::size_t
  {
    ISOTOPE_CORRELATION = 0, ///< Pearson of core M+0..M+3 intensities and the averagine pattern at the precursor mass.
    ISOTOPE_OVERLAP,         ///< Core M-1/z intensity over M+0 intensity, capped at 10 (a larger peak before the monoisotope).
    ISOTOPE_COELUTION,       ///< Mean Pearson of the M+1 and M+2 traces with the M+0 trace over the core MS1 bins.
    COUNT
  };
  constexpr std::size_t STAGE2_ISOTOPE_FEATURES = static_cast<std::size_t>(Stage2IsotopeFeature::COUNT);
  /// Isotope column names, in Stage2IsotopeFeature order.
  const std::array<std::string, STAGE2_ISOTOPE_FEATURES>& stage2IsotopeFeatureNames();
  /**
    @brief Isotope-envelope features of one window from its MS1 evidence.

    @param[in] evidence Window evidence from TraceCollector::enableEvidence().
    @param[in] precursor_mz Library precursor m/z.
    @param[in] charge Precursor charge; 0 or less gives all-zero features.
    @param[in] options Window geometry (only the RT and core half widths are used).
    @return Features in Stage2IsotopeFeature order; undefined correlations are zero.
  */
  std::array<double, STAGE2_ISOTOPE_FEATURES> stage2IsotopeFeatures(const TraceEvidence& evidence, double precursor_mz, int charge,
                                                                    const Stage2ExtractionOptions& options);

  /// Experimental peak-group boundary features: the window is re-picked with pickPeakGroup and scored over [left, right].
  enum class Stage2BoundaryFeature : std::size_t
  {
    LOG_FWHM = 0,          ///< log1p of the interpolated full width at half maximum in seconds.
    LOG_WIDTH,             ///< log1p of the picked base width in seconds.
    BOUNDARY_COELUTION,    ///< Mean pairwise Pearson of fragment traces over the picked peak (0 if fewer than 3 scans).
    BOUNDARY_LIBRARY_CORRELATION, ///< Pearson of sqrt fragment areas over the picked peak and sqrt library intensities.
    APEX_AGREEMENT,        ///< Share of library fragments whose maximum in the peak is within one scan of the consensus apex.
    COUNT
  };
  constexpr std::size_t STAGE2_BOUNDARY_FEATURES = static_cast<std::size_t>(Stage2BoundaryFeature::COUNT);
  /// Boundary column names, in Stage2BoundaryFeature order.
  const std::array<std::string, STAGE2_BOUNDARY_FEATURES>& stage2BoundaryFeatureNames();
  /**
    @brief Boundary features of one window, identical for targets and decoys and free of training.

    The ion-mobility-summed fragment traces are restricted to the RT bins that received an isolating MS2 scan (bins without a
    scan are not zero signal), the peak group around the window apex (centre bin) is picked with pickPeakGroup defaults, and
    width, co-elution, library correlation and apex agreement are measured over the picked boundaries.

    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats, fragment-major as in the collector.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @param[in] options Window geometry (only rt_half_width is used).
    @return Features in Stage2BoundaryFeature order; all zero if fewer than 3 scan bins or no signal.
  */
  std::array<double, STAGE2_BOUNDARY_FEATURES> stage2BoundaryFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                      const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                                                                      const Stage2ExtractionOptions& options);

  /// Experimental b/y ion evidence of the isolating MS2 spectrum nearest the window apex (all theoretical ions, not only the library's six).
  enum class Stage2ByIonFeature : std::size_t
  {
    BYION_B = 0,  ///< Number of matched singly charged b ions (TraceEvidence::byion_b).
    BYION_Y,      ///< Number of matched singly charged y ions (TraceEvidence::byion_y).
    COUNT
  };
  constexpr std::size_t STAGE2_BYION_FEATURES = static_cast<std::size_t>(Stage2ByIonFeature::COUNT);
  /// b/y ion column names, in Stage2ByIonFeature order.
  const std::array<std::string, STAGE2_BYION_FEATURES>& stage2ByIonFeatureNames();
  /**
    @brief b/y ion count features of one window; zero if no isolating MS2 spectrum fell inside the core window.
    @param[in] evidence Window evidence collected with TraceCollector::enableByIonCounts().
    @return Features in Stage2ByIonFeature order.
  */
  std::array<double, STAGE2_BYION_FEATURES> stage2ByIonFeatures(const TraceEvidence& evidence);

  /**
    @brief Experimental peak-shape and library-agreement features over the peak group picked as for the boundary group.

    Pairs are all pairs of library fragments; the lagged cross-correlation uses lags up to 3 scan bins on mean-centred
    traces normalised by their zero-lag norms. Library weights are sqrt(library intensity), normalised to sum 1.
  */
  enum class Stage2ShapeFeature : std::size_t
  {
    XCORR_COELUTION = 0,   ///< Mean |lag| of the best cross-correlation over fragment pairs (scan bins).
    XCORR_SHAPE,           ///< Mean best cross-correlation over fragment pairs.
    XCORR_COELUTION_W,     ///< As XCORR_COELUTION, pairs weighted by the product of library weights.
    XCORR_SHAPE_W,         ///< As XCORR_SHAPE, pairs weighted by the product of library weights.
    XCORR_SHAPE_MIN,       ///< Worst best cross-correlation over fragment pairs.
    LIBRARY_DOTPRODUCT,    ///< Cosine of sqrt fragment areas and sqrt library intensities.
    LIBRARY_MANHATTAN,     ///< L1 distance of the sum-normalised sqrt areas and sqrt library intensities.
    LIBRARY_RMSD,          ///< RMS difference of the same normalised vectors.
    LIBRARY_SPECTRAL_ANGLE,///< 1 - 2 arccos(dot product) / pi.
    LOG_SN,                ///< Mean over fragments of log1p(apex intensity / median scan-bin intensity of the window).
    INTENSITY_SHARE,       ///< Summed fragment area inside the peak over the summed fragment signal of the window.
    SIGNAL_FRAGMENTS,      ///< Fragments whose maximum inside the peak exceeds 3 x their window median.
    MIN_FRAGMENT_SHARE,    ///< Smallest fragment area over the summed area inside the peak.
    COUNT
  };
  constexpr std::size_t STAGE2_SHAPE_FEATURES = static_cast<std::size_t>(Stage2ShapeFeature::COUNT);
  /// Shape column names, in Stage2ShapeFeature order.
  const std::array<std::string, STAGE2_SHAPE_FEATURES>& stage2ShapeFeatureNames();
  /**
    @brief Shape features of one window, identical for targets and decoys and free of training.
    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @return Features in Stage2ShapeFeature order; all zero if no peak is picked or fewer than two fragments are present.
  */
  std::array<double, STAGE2_SHAPE_FEATURES> stage2ShapeFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence);

  /**
    @brief Experimental interference-robust features from automated transition refinement.

    After EncyclopeDIA (Searle et al., Nat Commun 9:5128, 2018, Supplementary Figure 4): over the peak group picked as
    for the boundary group, widened by REFINE_PAD_BINS scan bins on each side, every fragment trace is smoothed once with
    [1 2 1]/4 and normalised to unit area. The consensus peak shape is the median normalised intensity at every scan bin
    over the fragments with signal. A fragment whose normalised trace has Pearson r >= 0.9 with the consensus is
    'quantitative', r >= 0.75 is 'usable', lower r is treated as interfered. The median is robust to a minority of
    interfered fragments, so these features count only fragments that share one elution shape, not all signal at the
    fragment m/z values. The best-reference correlation follows the DIA-NN idea of one reference elution profile: the
    fragment with the highest summed correlation to the others is the reference.
  */
  enum class Stage2RefineFeature : std::size_t
  {
    QUANTITATIVE = 0,      ///< Library fragments with r >= 0.9 against the consensus.
    USABLE,                ///< Library fragments with r >= 0.75 against the consensus.
    MEAN_R,                ///< Mean r against the consensus over library fragments (no signal counts as 0).
    LIBRARY_WEIGHT,        ///< sqrt(library intensity) share of the quantitative fragments.
    QUANT_AREA_SHARE,      ///< Smoothed area of the quantitative fragments over the smoothed area of all fragments.
    BEST_REFERENCE_R,      ///< Mean r of the other fragments with the best reference fragment (no signal counts as 0).
    COUNT
  };
  constexpr std::size_t STAGE2_REFINE_FEATURES = static_cast<std::size_t>(Stage2RefineFeature::COUNT);
  /// Scan bins added on each side of the picked peak for the refinement correlations.
  constexpr std::size_t REFINE_PAD_BINS = 2;
  /// Refine column names, in Stage2RefineFeature order.
  const std::array<std::string, STAGE2_REFINE_FEATURES>& stage2RefineFeatureNames();
  /**
    @brief Per-fragment Pearson r against the median consensus shape, in tensor channel order.

    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @return r per channel; padding channels, fragments without signal, and windows with no picked peak, fewer than 5
      scan bins in the range or fewer than 3 fragments with signal give 0 for every channel.
  */
  std::array<double, 6> refineFragmentCorrelations(const float* tensor, const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence);
  /**
    @brief Experimental interference mask: zero every fragment channel whose refinement r is below @p min_r.

    Uses refineFragmentCorrelations(); a window it cannot refine is left unchanged. A masked channel keeps its m/z, so
    later features count it as a library fragment without signal.
    @param[in,out] tensor One window of TraceCollector::SAMPLE_SIZE floats.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @param[in] min_r Channels with signal and r < min_r are zeroed.
    @return Number of channels zeroed.
  */
  std::size_t maskInterferedFragments(float* tensor, const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence, double min_r);
  /**
    @brief Refine features of one window, identical for targets and decoys and free of training.
    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats.
    @param[in] library_intensities Library intensity per fragment channel.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @return Features in Stage2RefineFeature order; all zero under the conditions of refineFragmentCorrelations().
  */
  std::array<double, STAGE2_REFINE_FEATURES> stage2RefineFeatures(const float* tensor, const std::array<double, 6>& library_intensities,
                                                                  const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence);

  /**
    @brief Experimental quantification of one window, naive and with automated transition refinement.

    Naive: the peak group picked as for the boundary group, area summed over all fragments. Refined (EncyclopeDIA): the
    peak group re-picked on the fragments with r >= @p pick_min_r, refinement r recomputed over the refined peak, and the
    area summed over the fragments with r >= @p quant_min_r. Areas are raw (unsmoothed) IM-summed scan-bin sums.
  */
  struct Stage2Quant
  {
    bool found = false;                     ///< False if no peak group was picked.
    double naive_left = 0.0, naive_right = 0.0;     ///< Naive boundaries, seconds relative to the window centre (bin centres).
    double refined_left = 0.0, refined_right = 0.0; ///< Refined boundaries, same units.
    double naive_area = 0.0;                ///< All fragments over the naive boundaries.
    double refined_area = 0.0;              ///< Fragments with r >= quant_min_r over the refined boundaries.
    std::size_t quantitative = 0;           ///< Number of fragments in refined_area.
    std::array<double, 6> fragment_area{};  ///< Per channel, over the refined boundaries (0 for padding).
    std::array<double, 6> fragment_r{};     ///< Per channel refinement r over the refined boundaries (0 without signal).
  };
  /**
    @brief Quantify one window (see Stage2Quant).
    @param[in] tensor One window of TraceCollector::SAMPLE_SIZE floats.
    @param[in] fragment_mz Fragment m/z per channel; zero marks padding.
    @param[in] evidence Window evidence (only ms2_mask is used).
    @param[in] options Window geometry (only rt_half_width is used).
    @param[in] pick_min_r Refinement cut for re-picking the peak group.
    @param[in] quant_min_r Refinement cut for the quantitative fragments.
    @return The quantities; all zero if no peak group is picked.
  */
  Stage2Quant stage2Quantify(const float* tensor, const std::array<double, 6>& fragment_mz, const TraceEvidence& evidence,
                             const Stage2ExtractionOptions& options, double pick_min_r = 0.75, double quant_min_r = 0.9);

  /**
    @brief Peak-group evidence of one search candidate, read from the OpenDIA -features_out table.

    These come from the candidate search, not from the CNN tensor, and are computed identically for targets and decoys.
  */
  struct SearchPeakFeatures
  {
    std::string precursor_id;
    double rt_apex = 0.0;
    double transition_support = 0.0;   ///< Fragments that formed the peak group.
    double matched_scans = 0.0;        ///< Scans that contributed points.
    double peak_width_seconds = 0.0;
    double abs_rt_error_seconds = 0.0; ///< |apex - (calibrated) library RT|; 0 without a library RT.
  };
  /// Number of search columns appended by the SEARCH sets.
  constexpr std::size_t STAGE2_SEARCH_FEATURES = 4;
  /// Search column names: transition_support, log_matched_scans, log_peak_width, abs_library_rt_error.
  const std::array<std::string, STAGE2_SEARCH_FEATURES>& stage2SearchFeatureNames();
  /// The four search columns of one candidate: support, log1p(scans), log1p(width), |RT error|.
  std::array<double, STAGE2_SEARCH_FEATURES> stage2SearchFeatures(const SearchPeakFeatures& peak);

  /**
    @brief Read the OpenDIA -features_out candidate table, keyed by candidate_row (the window sample_id).

    @param[in] path TSV with candidate_row, precursor_id, rt_apex, peak_width_seconds, transition_support,
      matched_scans, has_library_rt and apex_rt_error_seconds columns.
    @return Map from candidate_row to its peak-group features.
    @throws std::invalid_argument If the file is missing, a column is absent, a value is not finite or a row is duplicated.
  */
  std::unordered_map<std::string, SearchPeakFeatures> readSearchPeakFeatures(const std::string& path);

  /// One window entering the second-stage discriminant. Every row of one call must have the same feature count.
  struct Stage2Row
  {
    std::string precursor_id;
    bool is_decoy = false;
    int fold = 0;
    std::vector<double> features;
  };

  struct Stage2Options
  {
    int iterations = 3;            ///< Semi-supervised LDA rounds per fold.
    double initial_fdr = 0.15;     ///< Positive cut in round 1, ranked by feature 0 (PyProphet ss_initial_fdr).
    double iteration_fdr = 0.05;   ///< Positive cut in later rounds (PyProphet ss_iteration_fdr).
    double ridge = 1e-3;           ///< Added to the diagonal of the standardised within-class covariance.
    std::size_t min_positives = 10;
  };

  /// Fitted weights of one held-out fold, on standardised features.
  struct Stage2FoldReport
  {
    int fold = 0;
    std::size_t positives = 0;
    std::size_t negatives = 0;
    int iterations = 0;
    std::vector<double> weights;
  };

  /**
    @brief Cross-fitted semi-supervised linear discriminant (PyProphet-style) over the rows.

    For each fold k, the discriminant is trained only on rows of the other folds:
    rank by feature 0, take precursor winners, use target winners at decoy-count q <= initial_fdr as
    positives and every decoy winner as negatives, fit a ridge LDA on standardised features, rescore and
    repeat with iteration_fdr. The final weights score fold k. Scores of different folds are therefore
    on different scales; standardise per fold before pooling (aggregatePairedCNNQValues pooled).

    @param[in] rows Windows with folds 0..folds-1.
    @param[in] folds Number of folds.
    @param[in] options Iteration settings.
    @param[out] report One entry per fold.
    @return One score per row, in input order.
    @throws std::invalid_argument If folds, features (non-finite, or differing counts) or a training set are invalid (too few positives or negatives).
  */
  std::vector<double> crossfitStage2LDA(const std::vector<Stage2Row>& rows, int folds, const Stage2Options& options,
                                        std::vector<Stage2FoldReport>& report);
}
