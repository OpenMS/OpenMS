// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/FEATUREFINDER/FFIDAlgoExternalIDHandler.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>

#include <vector>
#include <fstream>
#include <map>

namespace OpenMS {
/**
  @brief ID-guided MS1 feature finder; the algorithm behind @ref TOPP_FeatureFinderIdentification.

  Uses peptide identifications as targets to drive chromatogram extraction from MS1 data,
  builds a per-peptide assay library on the fly, scores extracted peak groups with
  @ref MRMFeatureFinderScoring, and emits a @ref FeatureMap of quantified features.

  @section FFid_inputs Inputs

  Two ID layers are accepted:
    - @em Internal IDs (@p peptides / @p proteins) — the high-confidence anchors used to
      define candidate elution windows and to train the SVM classifier (when classification
      is enabled).
    - @em External IDs (@p peptides_ext / @p proteins_ext) — optional, lower-confidence
      transfer IDs (e.g. from match-between-runs) used to extend the search space without
      contributing to SVM training. When external lists are empty, the algorithm skips
      machine learning and FDR estimation entirely.

  Optional @p seeds from an upstream untargeted feature finder can also be added.

  @section FFid_faims FAIMS handling

  If the input @ref PeakMap carries multiple FAIMS compensation voltages, @ref run splits
  the data via @ref IMDataConverter::splitByFAIMSCV, runs the per-CV pipeline in a fresh
  algorithm instance per group (peptide IDs filtered by @c FAIMS_CV through
  @ref FAIMSHelper::filterPeptidesByFAIMSCV), annotates the resulting features with the
  @c FAIMS_CV meta value, and merges everything back into the caller's @p features.
  IDs without a @c FAIMS_CV are processed against @em every CV group for backward
  compatibility. In the multi-FAIMS case, @ref getLibrary returns an empty library because
  each CV group has its own assay library.

  @section FFid_sideeffects ID-list transformations inside run()

  Both ID lists are taken by value (see the run() signature), so the @em caller's
  containers are never observably modified. The transformations described here happen
  on the function's local copies and influence which IDs reach @p features:

  - The local copy of @p peptides is shrunk to the best hit per identification.
  - FFid meta values are added to each surviving hit before it is written into
    @p features.
  - If @p seeds carries IDs, those are appended to the local copy of @p peptides
    before scoring.

  @section FFid_primary Primary MS-run path

  After feature detection the @ref FeatureMap's @c primaryMSRunPath is set to the
  @c ms_data_'s recorded run path; if that is not a valid / readable mzML, @p spectra_file
  is used as a fallback annotation (see @ref MSExperiment overload of
  @ref FeatureMap::setPrimaryMSRunPath).

  @ingroup FeatureFinder
*/
class OPENMS_DLLAPI FeatureFinderIdentificationAlgorithm :
  public DefaultParamHandler
{
public:
  /// Default constructor; installs the FFid parameters (see class docs)
  FeatureFinderIdentificationAlgorithm();

  /**
    @brief Run the FFid pipeline; for FAIMS data this dispatches one run per CV group and merges results.

    See the class brief for the role of external IDs, the seeds list, the FAIMS handling,
    the side effects on the input ID lists, and the primary-MS-run-path fallback semantics.

    @param[in]  peptides      Internal peptide IDs (taken by value). The local copy is shrunk to best hit per identification and FFid meta values are added before being written into @p features; the caller's container is not modified.
    @param[in]  proteins      Protein IDs corresponding to @p peptides.
    @param[in]  peptides_ext  External peptide IDs, optional (may be empty); taken by value with the same shrink-and-annotate treatment as @p peptides.
    @param[in]  proteins_ext  Protein IDs corresponding to @p peptides_ext.
    @param[out] features      Quantified feature map; pre-existing contents are cleared for FAIMS data and replaced.
    @param[in]  seeds         Optional pre-detected features from an untargeted feature finder.
    @param[in]  spectra_file  Source mzML path used as a fallback for @c primaryMSRunPath annotation when the MSExperiment's own path isn't usable.
  */
  void run(
    PeptideIdentificationList peptides,
    const std::vector<ProteinIdentification>& proteins,
    PeptideIdentificationList peptides_ext,
    std::vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds = FeatureMap(),
    const String& spectra_file = ""
    );

  /// Re-score / filter an existing candidate @ref FeatureMap in place using the configured classifier and quality cutoffs (entry point used to re-process candidates exported via the @c candidates_out parameter).
  void runOnCandidates(FeatureMap& features);

  /**
    @brief Mutable access to the cached MS1 data.
    @return Reference to the cached @ref PeakMap (already pre-processed by @ref run when called after run()).
  */
  PeakMap& getMSData();
  /**
    @brief Read-only access to the cached MS1 data.
    @return Const reference to the cached @ref PeakMap.
  */
  const PeakMap& getMSData() const;

  /**
    @brief Copy the MS data into the algorithm instance; useful from pyOpenMS where moving is awkward.
    @param[in] ms_data Source @ref PeakMap; copied into the internal cache.
  */
  void setMSData(const PeakMap& ms_data); // for pyOpenMS
  /**
    @brief Move the MS data into the algorithm instance (no copy).
    @param[in] ms_data Source @ref PeakMap; moved into the internal cache. @ref getMSData will later expose the cached (possibly preprocessed) copy, not the original.
  */
  void setMSData(PeakMap&& ms_data); // moves peak data and saves the copy. Note that getMSData() will give back a processed/modified version.

  /**
    @brief Mutable access to the accumulated extracted chromatograms (XICs) from the last @ref run.
    @return Reference to the accumulated XIC @ref PeakMap.
  */
  PeakMap& getChromatograms();
  /**
    @brief Read-only access to the accumulated extracted chromatograms.
    @return Const reference to the accumulated XIC @ref PeakMap.
  */
  const PeakMap& getChromatograms() const;

  /**
    @brief Mutable access to the progress logger used by the algorithm.
    @return Reference to the @ref ProgressLogger (configure log type, etc.).
  */
  ProgressLogger& getProgressLogger();
  /**
    @brief Read-only access to the progress logger.
    @return Const reference to the @ref ProgressLogger.
  */
  const ProgressLogger& getProgressLogger() const;

  /**
    @brief Mutable access to the assay library used / produced by the last @ref run.
    @return Reference to the @ref TargetedExperiment library.
    @note For multi-FAIMS inputs this returns an empty library because each CV group has its own internal library.
  */
  TargetedExperiment& getLibrary();
  /**
    @brief Read-only access to the assay library used / produced by the last @ref run.
    @return Const reference to the @ref TargetedExperiment library (empty for multi-FAIMS inputs — see the non-const overload).
  */
  const TargetedExperiment& getLibrary() const;

protected:

  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  /// mapping: RT (not necessarily unique) -> pointer to peptide
  typedef std::multimap<double, PeptideIdentification*> RTMap;
  /// mapping: charge -> internal/external: (RT -> pointer to peptide)
  typedef std::map<Int, std::pair<RTMap, RTMap> > ChargeMap;
  /// mapping: sequence -> charge -> internal/external ID information
  typedef std::map<AASequence, ChargeMap> PeptideMap;
  /// mapping: peptide ref. -> int./ext.: (RT -> pointer to peptide)
  typedef std::map<String, std::pair<RTMap, RTMap> > PeptideRefRTMap;

  PeptideMap peptide_map_;

  Size n_internal_peps_; ///< number of internal peptide
  Size n_external_peps_; ///< number of external peptides

  Size batch_size_; ///< nr of peptides to use at the same time during chromatogram extraction
  double rt_window_; ///< RT window width
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?

  double mapping_tolerance_; ///< RT tolerance for mapping IDs to features

  double isotope_pmin_; ///< min. isotope probability for peptide assay
  Size n_isotopes_; ///< number of isotopes for peptide assay

  double rt_quantile_;

  double peak_width_;
  double min_peak_width_;
  double signal_to_noise_;

  String elution_model_;

  // SVM related parameters
  double svm_min_prob_;
  StringList svm_predictor_names_;
  String svm_xval_out_;
  double svm_quality_cutoff;
  Size svm_n_parts_; ///< number of partitions for SVM cross-validation
  Size svm_n_samples_; ///< number of samples for SVM training

  // output file (before filtering)
  String candidates_out_;

  Size debug_level_;

  void updateMembers_() override;

  /// region in RT in which a peptide elutes:
  struct RTRegion
  {
    double start, end;
    ChargeMap ids; ///< internal/external peptide IDs (per charge) in this region
  };

  /**
   * @brief Ion mobility statistics for a peptide in a specific RT region and charge state
   *
   * This structure stores statistical measures of ion mobility values collected from
   * peptide identifications within a single RT region. These statistics are used for:
   * - Chromatogram extraction with appropriate IM windows (using median)
   * - Feature annotation for quality control (all three values)
   * - Detecting potential mis-identifications (large min-max spread may indicate issues)
   *
   * All values default to -1.0 to indicate missing/unavailable IM data.
   */
  struct IMStats
  {
    double median = -1.0; ///< Median IM value (robust central tendency, used for extraction)
    double min = -1.0;    ///< Minimum IM value (lower bound of IM distribution)
    double max = -1.0;    ///< Maximum IM value (upper bound of IM distribution)
  };

  /// predicate for filtering features by overall quality:
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.getOverallQuality() == 0.0;
    }
  } feature_filter_quality_;

  /// predicate for filtering features by assigned peptides:
  struct FeatureFilterPeptides
  {
    bool operator()(const Feature& feature)
    {
      return feature.getPeptideIdentifications().empty();
    }
  } feature_filter_peptides_;

  /// comparison functor for (unassigned) peptide IDs
  struct PeptideCompare
  {
    bool operator()(const PeptideIdentification& p1,
                    const PeptideIdentification& p2)
    {
      const String& seq1 = p1.getHits()[0].getSequence().toString();
      const String& seq2 = p2.getHits()[0].getSequence().toString();
      if (seq1 == seq2)
      {
        Int charge1 = p1.getHits()[0].getCharge();
        Int charge2 = p2.getHits()[0].getCharge();
        if (charge1 == charge2)
        {
          return p1.getRT() < p2.getRT();
        }
        return charge1 < charge2;
      }
      return seq1 < seq2;
    }
  } peptide_compare_;

  /// comparison functor for features
  struct FeatureCompare
  {
    bool operator()(const Feature& f1, const Feature& f2)
    {
      const String& ref1 = f1.getMetaValue("PeptideRef");
      const String& ref2 = f2.getMetaValue("PeptideRef");
      if (ref1 == ref2)
      {
        return f1.getRT() < f2.getRT();
      }
      return ref1 < ref2;
    }
  } feature_compare_;

  PeakMap ms_data_; ///< input LC-MS data
  PeakMap chrom_data_; ///< accumulated chromatograms (XICs)
  TargetedExperiment library_; ///< assays for peptides (cleared per chunk during processing)
  TargetedExperiment output_library_; ///< accumulated assays for output (populated from library_ before clearing)

  bool quantify_decoys_;
  double add_mass_offset_peptides_{0.0}; ///< non-zero if for every feature an additional offset features should be extracted
  double seed_apex_rt_tolerance_{5.0}; ///< max allowed RT deviation (s) between seed apex and detected feature apex; seed-derived features exceeding this are removed (0 = disabled)
  bool use_psm_cutoff_;
  double psm_score_cutoff_;
  PeptideIdentificationList unassignedIDs_;

  const double seed_rt_window_ = 60.0; ///< extraction window used for seeds (smaller than rt_window_ as we know the exact apex positions)

  /// SVM probability -> number of pos./neg. features (for FDR calculation):
  std::map<double, std::pair<Size, Size> > svm_probs_internal_;
  /// SVM probabilities for "external" features (for FDR calculation):
  std::multiset<double> svm_probs_external_;
  Size n_internal_features_; ///< internal feature counter (for FDR calculation)
  Size n_external_features_; ///< external feature counter (for FDR calculation)
  /// TransformationDescription trafo_; // RT transformation (to range 0-1)
  std::map<String, double> isotope_probs_; ///< isotope probabilities of transitions
  /**
   * @brief Ion mobility statistics per peptide reference (peptide sequence/charge:region)
   *
   * Maps from full peptide reference (e.g., "PEPTIDE/2:1") to IM statistics.
   * Populated during createAssayLibrary_() and used during annotateFeatures_()
   * to add IM_median, IM_min, and IM_max meta-values to features.
   */
  std::map<String, IMStats> im_stats_;

  /**
   * @brief Global ion mobility statistics from all peptide identifications
   *
   * Calculated from peptide identifications BEFORE seeds are added (ensuring we only
   * learn from real IDs with IM annotation). Provides context for the typical IM range
   * in the dataset.
   */
  IMStats global_im_stats_;

  MRMFeatureFinderScoring feat_finder_; ///< OpenSWATH feature finder
  Internal::FFIDAlgoExternalIDHandler external_id_handler_; ///< Handler for external peptide IDs

  ProgressLogger prog_log_;

  /// generate transitions (isotopic traces) for a peptide ion and add them to the library:
  void generateTransitions_(const String& peptide_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist);

  void addPeptideRT_(TargetedExperiment::Peptide& peptide, double rt) const;

  /// get regions in which peptide eludes (ideally only one) by clustering RT elution times
  void getRTRegions_(ChargeMap& peptide_data, std::vector<RTRegion>& rt_regions, bool clear_IDs = true) const;

  /**
   * @brief Calculate ion mobility statistics for peptide identifications in an RT region
   *
   * Computes median, min, and max IM values from peptide identifications within
   * the given RT region (across all charge states). Individual IDs lacking IM
   * annotation are skipped (with warning), and statistics are calculated from the
   * remaining IDs with valid IM data. The median is used for robust central tendency
   * estimation and is more resistant to outliers than the mean.
   *
   * Seeds from untargeted feature finders may or may not have an IM meta value set,
   * depending on the feature finder. If IM is annotated on the seed, it is used for
   * targeted extraction. If not, the seed is extracted across the full IM range
   * (ChromatogramExtractor disables IM filtering when ion_mobility < 0).
   *
   * Note: RT region boundaries are determined from ALL IDs (including those without IM),
   * so this only affects IM statistics calculation, not RT extraction.
   *
   * @param[in] r RT region containing peptide identifications grouped by charge state
   * @return IMStats structure with median/min/max, or {-1, -1, -1} if no valid IM data
   *
   * @see IMStats for details on the returned structure
   */
  IMStats getRTRegionIMStats_(const RTRegion& r);

  /**
   * @brief Calculate global IM statistics from MS data and peptide identifications
   *
   * Uses MSExperiment::getMinMobility()/getMaxMobility() to get the full IM range
   * from raw data (min/max), and calculates median from peptide identifications
   * for robust central tendency. Must be called BEFORE addSeeds_() to ensure
   * global statistics are based only on identified peptides.
   *
   * Seeds may or may not have IM annotation depending on the feature finder.
   * Seeds with IM annotation use their own IM value; seeds without IM are
   * extracted across the full IM range of the dataset.
   */
  void calculateGlobalIMStats_();

  void annotateFeaturesFinalizeAssay_(
    FeatureMap& features,
    std::map<Size, std::vector<PeptideIdentification*> >& feat_ids,
    RTMap& rt_internal);

  /// annotate identified features with m/z, isotope probabilities, etc.
  void annotateFeatures_(FeatureMap& features, PeptideRefRTMap& ref_rt_map);

  void ensureConvexHulls_(Feature& feature) const;

  void postProcess_(FeatureMap& features, bool with_external_ids);

  /// Helper functions for run()
  void validateSVMParameters_() const;
  void initializeFeatureFinder_();
  double calculateRTWindow_(double rt_uncertainty) const;
  void removeSeedPseudoIDs_(FeatureMap& features);

  /// Helper function to check if a peptide hit is a seed pseudo-ID
  static bool isSeedPseudoHit_(const PeptideHit& hit);

  /// Calculate RT bounds with optional tolerance expansion
  std::pair<double, double> calculateRTBounds_(double rt_min, double rt_max) const;

  /// some statistics on detected features
  void statistics_(const FeatureMap& features) const;

  /// creates an assay library out of the peptide sequences and their RT elution windows
  /// the PeptideMap is mutable since we clear it on-the-go
  /// @p clear_IDs set to false to keep IDs in internal charge maps (only needed for debugging purposes)
  void createAssayLibrary_(const PeptideMap::iterator& begin, const PeptideMap::iterator& end, PeptideRefRTMap& ref_rt_map, bool clear_IDs = true);

  /// CAUTION: This method stores a pointer to the given @p peptide reference in internals
  /// Make sure it stays valid until destruction of the class.
  /// @todo find better solution
  void addPeptideToMap_(PeptideIdentification& peptide,
    PeptideMap& peptide_map,
    bool external = false);

  void filterFeatures_(FeatureMap& features, bool classified);

  /// Core processing logic for a single (non-FAIMS or single FAIMS group) dataset
  /// Called by run() either directly or for each FAIMS CV group
  void runSingleGroup_(
    PeptideIdentificationList peptides,
    const std::vector<ProteinIdentification>& proteins,
    PeptideIdentificationList peptides_ext,
    std::vector<ProteinIdentification> proteins_ext,
    FeatureMap& features,
    const FeatureMap& seeds,
    const String& spectra_file);

  // seeds for untargeted extraction
  Size addSeeds_(PeptideIdentificationList& peptides, const FeatureMap& seeds);

  // quant. decoys
  Size addOffsetPeptides_(PeptideIdentificationList& peptides, double offset);

  /// Chunks an iterator range (allowing advance and distance) into batches of size batch_size.
  /// Last batch might be smaller.
  template <typename It>
  std::vector<std::pair<It,It>>
  chunk_(It range_from, It range_to, const std::ptrdiff_t batch_size)
  {
    /* Aliases, to make the rest of the code more readable. */
    using std::vector;
    using std::pair;
    using std::make_pair;
    using std::distance;
    using diff_t = std::ptrdiff_t;

    /* Total item number and batch_size size. */
    const diff_t total {distance(range_from, range_to)};
    const diff_t num {total / batch_size};

    vector<pair<It,It>> chunks(num);

    It batch_end {range_from};

    /* Use the 'generate' algorithm to create batches. */
    std::generate(begin(chunks), end(chunks), [&batch_end, batch_size]()
    {
      It batch_start {batch_end };

      std::advance(batch_end, batch_size);
      return make_pair(batch_start, batch_end);
    });

    /* The last batch_size's end must always be 'range_to'. */
    if (chunks.empty())
    {
      chunks.emplace_back(range_from, range_to);
    }
    else
    {
      chunks.back().second = range_to;
    }

    return chunks;
  }
}; // namespace OpenMS
} // namespace OpenMS
 