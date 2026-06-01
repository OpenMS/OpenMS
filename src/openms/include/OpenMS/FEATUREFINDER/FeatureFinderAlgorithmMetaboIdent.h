// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CHEMISTRY/AdductInfo.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>

#include <map>
#include <vector>

namespace OpenMS
{
class IsotopeDistribution;
  
class OPENMS_DLLAPI FeatureFinderAlgorithmMetaboIdent :
  public DefaultParamHandler
{
public:
  /**
    @brief One target compound for the @ref FeatureFinderAlgorithmMetaboIdent assay library.

    Plain data container. No validation is performed in the constructor or accessors —
    the contract below is enforced when the compound is consumed by
    @ref FeatureFinderAlgorithmMetaboIdent::run via the internal
    @c addTargetToLibrary_ helper. Compounds that violate the contract are silently
    dropped with an @c OPENMS_LOG_ERROR message at run time:
      - Either @c mass @c > @c 0 or a non-empty @c formula is required (one of the two must
        be present so the m/z can be derived).
      - @c rts must be non-empty.
      - For each entry in @c charges, the charge state must be non-zero.
      - @c rt_ranges (when provided) and @c ion_mobilities (when provided) must have
        either @c 1 entry (broadcast to every RT) or exactly @c rts.size() entries.

    The struct stores the per-target data verbatim; the downstream library builder fills in
    defaults (e.g. RT range @c 0.0, ion mobility sentinel @c -1.0 for "no IM filter") and
    computes the isotope distribution from @c formula when @c iso_distrib is empty or its
    first entry is @c 0.
  */
  struct OPENMS_DLLAPI FeatureFinderMetaboIdentCompound
  {
    /**
      @brief Construct a target compound by aggregating its identification data.

      @param[in] _name             Target name (used as a tag in the resulting assay library).
      @param[in] _formula          Molecular formula (e.g. @c "C12H22O11"); empty when only @p _mass is provided.
      @param[in] _mass             Monoisotopic mass; @c <= @c 0 means "derive from @p _formula".
      @param[in] _charges          Charge states to register for this target; each non-zero entry produces one library transition.
      @param[in] _rts              Expected retention times (one library entry per RT). Must be non-empty when the compound is fed into @ref FeatureFinderAlgorithmMetaboIdent::run.
      @param[in] _rt_ranges        RT tolerance per RT (parallel to @p _rts); a one-element vector is broadcast, an empty vector uses the algorithm default.
      @param[in] _iso_distrib      Pre-computed isotope-intensity distribution; empty (or first entry @c 0) triggers computation from @p _formula at library-build time.
      @param[in] _ion_mobilities   Optional ion-mobility values (parallel to @p _rts); a one-element vector is broadcast, an empty vector disables IM filtering for this target.
      @param[in] _adduct           Adduct string (e.g. @c "M+H;1+", @c "M+Na;1+"); empty lets downstream logic infer defaults from charge polarity.
    */
    FeatureFinderMetaboIdentCompound(const String& _name,
        const String& _formula,
        double _mass,
        const std::vector<int>& _charges,
        const std::vector<double>& _rts,
        const std::vector<double>& _rt_ranges,
        const std::vector<double>& _iso_distrib,
        const std::vector<double>& _ion_mobilities = {},
        const String& _adduct = ""):
      name_(_name),
      formula_(_formula),
      mass_(_mass),
      charges_(_charges),
      rts_(_rts),
      rt_ranges_(_rt_ranges),
      iso_distrib_(_iso_distrib),
      ion_mobilities_(_ion_mobilities),
      adduct_(_adduct)
      {
      }

    private:
      String name_;                          ///< Target name; used as the @c name meta value and in the auto-generated library ID.
      String formula_;                       ///< Molecular formula; empty if the caller only provides @c mass_.
      double mass_;                          ///< Monoisotopic mass; @c <= @c 0 means "compute from @c formula_".
      std::vector<int> charges_;             ///< Charge states to register; entries @c == @c 0 are dropped at library-build time.
      std::vector<double> rts_;              ///< Expected retention times; one library transition per entry.
      std::vector<double> rt_ranges_;        ///< RT tolerance per @c rts_ entry; one-element broadcasts, empty falls back to the algorithm default.
      std::vector<double> iso_distrib_;      ///< Optional pre-computed isotope-intensity distribution; empty (or @c [0]) triggers computation from @c formula_.
      std::vector<double> ion_mobilities_;   ///< Optional ion-mobility values per @c rts_ entry; one-element broadcasts, empty disables IM filtering.
      String adduct_;                        ///< Adduct string (e.g. @c "M+H;1+", @c "M+Na;1+", @c "M-H;1-"). Empty lets downstream logic infer defaults from charge polarity.

    public:
      /// Return the target name.
      const String& getName() const {
        return name_;
      }

      /// Return the molecular formula (may be empty if only @c mass was provided).
      const String& getFormula() const {
        return formula_;
      }

      /// Return the monoisotopic mass; @c <= @c 0 means "compute from @c formula at library-build time".
      double getMass() const {
        return mass_;
      }

      /// Return the list of charge states to register for this target.
      const std::vector<int>& getCharges() const {
        return charges_;
      }

      /// Return the expected retention times (one library transition per entry).
      const std::vector<double>& getRTs() const {
        return rts_;
      }

      /// Return the RT tolerance per @c rts_ entry (or a single broadcast value, or empty for the algorithm default).
      const std::vector<double>& getRTRanges() const {
        return rt_ranges_;
      }

      /// Return the pre-computed isotope-intensity distribution (or empty / @c [0] for runtime computation from @c formula).
      const std::vector<double>& getIsotopeDistribution() const {
        return iso_distrib_;
      }

      /// Return the optional ion-mobility values per @c rts_ entry (or a single broadcast value, or empty to disable IM filtering).
      const std::vector<double>& getIonMobilities() const {
        return ion_mobilities_;
      }

      const String& getAdduct() const {
        return adduct_;
      }
  };

  /// default constructor
  FeatureFinderAlgorithmMetaboIdent();

  /// @brief perform targeted feature extraction of compounds from @p metaboIdentTable and stores them in @p features.
  /// If @p spectra_file is provided it will be used as a fall-back to setPrimaryMSRunPath
  /// in the feature map in case a proper primaryMSRunPath is not annotated in the MSExperiment.
  /// If there are no MS1 scans in the MSData return @p features unchanged.
  ///
  /// FAIMS data is handled automatically: if the MS data contains multiple FAIMS compensation
  /// voltages, each CV group is processed independently and results are combined with
  /// FAIMS_CV annotation on features. For multi-FAIMS data, getLibrary() returns an empty
  /// library since each FAIMS group has its own assay library.
  void run(const std::vector<FeatureFinderMetaboIdentCompound>& metaboIdentTable, FeatureMap& features, const String& spectra_file = "");

  /// @brief Retrieve chromatograms (empty if run was not executed)
  PeakMap& getMSData() { return ms_data_; }
  const PeakMap& getMSData() const { return ms_data_; }

  /// @brief Set spectra
  void setMSData(const PeakMap& m); // needed because pyOpenMS can't wrap the non-const reference version

  void setMSData(PeakMap&& m); // moves peak data and saves the copy. Note that getMSData() will give back a processed/modified version.

  /// @brief Retrieve chromatograms (empty if run was not executed)
  const PeakMap& getChromatograms() const { return chrom_data_; }
  PeakMap& getChromatograms() { return chrom_data_; }

  /// @brief Retrieve the assay library (e.g., to store as TraML, empty if run was not executed)
  const TargetedExperiment& getLibrary() const { return library_; }
  
  /// @brief Retrieve deviations between provided coordinates and extracted ones (e.g., to store as TrafoXML or for plotting)
  const TransformationDescription& getTransformations() const { return trafo_; }

  /// @brief Retrieve number of features with shared identifications
  size_t getNShared() const  { return n_shared_; }

  static String prettyPrintCompound(const TargetedExperiment::Compound& compound);
protected:

  /// Boundaries for a mass trace in a feature
  struct MassTraceBounds
  {
    Size sub_index;
    double rt_min, rt_max, mz_min, mz_max;
  };

  /// Boundaries for all mass traces per feature
  typedef std::map<UInt64, std::vector<MassTraceBounds> > FeatureBoundsMap;

  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  typedef std::vector<Feature*> FeatureGroup; ///< group of (overlapping) features

  /// Predicate for filtering features by overall quality
  struct FeatureFilterQuality
  {
    bool operator()(const Feature& feature)
    {
      return feature.metaValueExists("FFMetId_remove");
    }
  } feature_filter_;

  /// Comparison functor for features
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


  void extractTransformations_(const FeatureMap& features);

  /// Add a target (from the input file) to the assay library
  void addTargetToLibrary_(const String& name, const String& formula,
                           double mass, const std::vector<Int>& charges,
                           const std::vector<double>& rts,
                           std::vector<double> rt_ranges,
                           const std::vector<double>& iso_distrib,
                           const std::vector<double>& ion_mobilities = {},
                           const String& adduct = "");

  /// Add "peptide" identifications with information about targets to features
  Size addTargetAnnotations_(FeatureMap& features);

  void addTargetRT_(TargetedExperiment::Compound& target, double rt);

  /// Calculate mass-to-charge ratio from mass and charge
  double calculateMZ_(double mass, Int charge) const;

  void generateTransitions_(const String& target_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist);

  void annotateFeatures_(FeatureMap& features);

  void ensureConvexHulls_(Feature& feature) const;

  void selectFeaturesFromCandidates_(FeatureMap& features);

  /// Core processing logic for a single (non-FAIMS or single FAIMS group) dataset
  /// Called by run() either directly or for each FAIMS CV group
  void runSingleGroup_(const std::vector<FeatureFinderMetaboIdentCompound>& metaboIdentTable,
                       FeatureMap& features,
                       const String& spectra_file);

  double rt_window_; ///< RT window width
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?
  double im_window_; ///< Ion mobility window width (0 = disabled)

  double isotope_pmin_; ///< min. isotope probability for peptide assay
  Size n_isotopes_; ///< number of isotopes for peptide assay

  double peak_width_;
  double min_peak_width_;
  double signal_to_noise_;

  String elution_model_;

  // output file (before filtering)
  String candidates_out_;

  Size debug_level_;

  void updateMembers_() override;

  PeakMap ms_data_; ///< input LC-MS data
  PeakMap chrom_data_; ///< accumulated chromatograms (XICs)

  MRMFeatureFinderScoring feat_finder_; ///< OpenSWATH feature finder

  TargetedExperiment library_; ///< accumulated assays for targets
  
  TransformationDescription trafo_;
  
  CoarseIsotopePatternGenerator iso_gen_; ///< isotope pattern generator
  std::map<String, double> isotope_probs_; ///< isotope probabilities of transitions
  std::map<String, double> target_rts_; ///< RTs of targets (assays)
  
  size_t n_shared_ = 0;
};

} // namespace OpenMS
