// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERIDENTIFICATIONALGORITHM_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERIDENTIFICATIONALGORITHM_H

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureFinderScoring.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmPickedHelperStructs.h>

#include <vector>
#include <fstream>
#include <map>

namespace OpenMS
{
  class IsotopeDistribution;


class OPENMS_DLLAPI FeatureFinderIdentificationAlgorithm :
  public DefaultParamHandler
{
public:
  /// default constructor
  FeatureFinderIdentificationAlgorithm();

  /**
     @brief Run feature detection

     @param features Output feature map
     @param id_data Primary ("internal") identifications as targets for feature detection
     @param id_data_ext Additional ("external") identifications as targets for feature detection
     @param spectra_file Fall-back value for setting @p primaryMSRunPath in the output (by default set based on the MS data being processed)

     External IDs (@p id_data_ext) may be empty, in which case no machine learning or FDR estimation will be performed.
  */
  void run(FeatureMap& features,
           IdentificationData& id_data,
           IdentificationData& id_data_ext,
           const String& spectra_file = "");

  /// Convert seeds to an IdentificationData representation
  void convertSeeds(const FeatureMap& seeds, IdentificationData& id_data,
                    Size n_overlap_traces = 6);

  // void runOnCandidates(FeatureMap& features);

  PeakMap& getMSData();
  const PeakMap& getMSData() const;

  /// @brief set the MS data used for feature detection
  void setMSData(const PeakMap& ms_data); // for pyOpenMS
  void setMSData(PeakMap&& ms_data); // moves peak data and saves the copy. Note that getMSData() will give back a processed/modified version.

  PeakMap& getChromatograms();
  const PeakMap& getChromatograms() const;

  ProgressLogger& getProgressLogger();
  const ProgressLogger& getProgressLogger() const;

  TargetedExperiment& getLibrary();
  const TargetedExperiment& getLibrary() const;

protected:
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  // aggregate all search hits (internal and external) grouped by molecule (e.g.
  // peptide) and charge state, ordered by RT:
  /// mapping: RT (not necessarily unique) -> reference to search hit
  typedef std::multimap<double, IdentificationData::ObservationMatchRef> RTMap;
  /// mapping: charge -> internal/external: (RT -> ref. to search hit)
  typedef std::map<Int, std::pair<RTMap, RTMap>> ChargeMap;

  struct TargetData
  {
    IdentificationData::IdentifiedMolecule molecule;
    IdentificationData::AdductOpt adduct;
    ChargeMap hits_by_charge;
  };
  /// mapping: target ion ID -> associated data
  typedef std::map<String, TargetData> TargetMap;

  /// region in RT in which a target elutes:
  struct RTRegion
  {
    double start, end;
    ChargeMap ids; ///< internal/external IDs (per charge) in this region
  };

  /// comparison functor for features
  struct FeatureCompare
  {
    bool operator()(const Feature& f1, const Feature& f2)
    {
      const String& ref1 = f1.getMetaValue("CompoundRef");
      const String& ref2 = f2.getMetaValue("CompoundRef");
      if (ref1 == ref2)
      {
        return f1.getRT() < f2.getRT();
      }
      return ref1 < ref2;
    }
  } feature_compare_;

  TargetMap target_map_; ///< aggregated IDs for each identified molecule

  Size n_internal_targets_; ///< number of internal target molecules
  Size n_external_targets_; ///< number of external target molecules
  Size n_seed_targets_; ///< number of targets derived from seeds

  Size batch_size_; ///< number of target molecules to consider together during chromatogram extraction
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?
  double rt_window_; ///< RT window width (for "proper" IDs)
  double rt_window_seeds_; ///< RT window width for seeds

  double mapping_tolerance_; ///< RT tolerance for mapping IDs to features

  Size n_isotopes_; ///< number of isotopes for assay
  bool max_isotopes_; ///< consider most abundant isotopes?
  double isotope_pmin_; ///< min. isotope probability for assay

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

  PeakMap ms_data_; ///< input LC-MS data
  PeakMap chrom_data_; ///< accumulated chromatograms (XICs)
  TargetedExperiment library_; ///< accumulated assays for targets (one chunk)
  TargetedExperiment combined_library_; ///< accumulated assays for targets (all chunks)

  bool quantify_decoys_;

  /// SVM probability -> number of pos./neg. features (for FDR calculation):
  std::map<double, std::pair<Size, Size> > svm_probs_internal_;
  /// SVM probabilities for "external" features (for FDR calculation):
  std::multiset<double> svm_probs_external_;
  Size n_internal_features_; ///< internal feature counter (for FDR calculation)
  Size n_external_features_; ///< external feature counter (for FDR calculation)
  /// TransformationDescription trafo_; // RT transformation (to range 0-1)
  TransformationDescription trafo_external_; ///< transform. to external RT scale
  std::map<String, double> isotope_probs_; ///< isotope probabilities of transitions
  MRMFeatureFinderScoring feat_finder_; ///< OpenSWATH feature finder

  ProgressLogger prog_log_;

  void updateMembers_() override;

  /// generate transitions (isotopic traces) for an ion and add them to the library:
  void generateTransitions_(const String& target_id, double target_mass,
                            Int charge, const IsotopeDistribution& iso_dist);

  void addTargetRT_(TargetedExperiment::Compound& target, double rt) const;

  /// get regions in which target elutes (ideally only one) by clustering RT elution times
  void makeRTRegions_(const ChargeMap& charge_data, std::vector<RTRegion>& rt_regions,
                      bool is_seed = false) const;

  /// annotate identified features with m/z, isotope probabilities, etc.
  void annotateFeatures_(FeatureMap& features);

  void annotateFeaturesOneTarget_(FeatureMap& features, const String& target_id,
                                  Int charge, const std::vector<Size>& indexes);

  void ensureConvexHulls_(Feature& feature);

  void postProcess_(FeatureMap& features, bool with_external_ids);

  /// print some statistics on detected features
  void statistics_(const FeatureMap& features, bool with_external_ids) const;

  /*!
    @brief Creates an assay library given target molecule information

    @p MoleculeMap will be (partially) cleared and thus has to be mutable.
  */
  void createAssayLibrary_(TargetMap::iterator begin, TargetMap::iterator end);

  void addMatchToTargetMap_(IdentificationData::ObservationMatchRef ref, bool external = false);

  void checkNumObservations_(Size n_pos, Size n_neg, const String& note = "") const;

  void getUnbiasedSample_(const std::multimap<double, std::pair<Size, bool>>& valid_obs,
                          std::map<Size, Int>& training_labels);

  void getRandomSample_(std::map<Size, Int>& training_labels);

  void classifyFeatures_(FeatureMap& features);

  void filterFeaturesFinalizeAssay_(Feature& best_feature, double best_quality,
                                    const double quality_cutoff, const String& target_id);

  void filterFeatures_(FeatureMap& features, bool classified);

  void calculateFDR_(FeatureMap& features);

  static std::pair<String, Int> extractTargetID_(const Feature& feature, bool extract_charge = false);

  /// Chunks an iterator range (allowing advance and distance) into batches of size @p batch_size.
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
};

} // namespace OpenMS

#endif
