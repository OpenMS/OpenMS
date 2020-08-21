// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

  /// Main method for actual FeatureFinder
  /// External IDs (@p id_data_ext) may be empty, in which case no machine learning or FDR estimation will be performed.
  void run(
    FeatureMap& features,
    IdentificationData& id_data,
    IdentificationData& id_data_ext);

  /// Convert seeds to an IdentificationData representation
  void convertSeeds(const FeatureMap& seeds, IdentificationData& id_data,
                    Size n_overlap_traces = 6);

  // void runOnCandidates(FeatureMap& features);

  PeakMap& getMSData() { return ms_data_; }
  const PeakMap& getMSData() const { return ms_data_; }

  PeakMap& getChromatograms() { return chrom_data_; }
  const PeakMap& getChromatograms() const { return chrom_data_; }

  ProgressLogger& getProgressLogger() { return prog_log_; }
  const ProgressLogger& getProgressLogger() const { return prog_log_; }

  // @TODO: how does this work if the library is cleared between chunks?
  TargetedExperiment& getLibrary() { return combined_library_; }
  const TargetedExperiment& getLibrary() const { return combined_library_; }

protected:
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTrace MassTrace;
  typedef FeatureFinderAlgorithmPickedHelperStructs::MassTraces MassTraces;

  typedef std::pair<IdentificationData::IdentifiedMolecule,
                    boost::optional<IdentificationData::AdductRef>> AdductedID;

  // aggregate all search hits (internal and external) grouped by molecule (e.g.
  // peptide) and charge state, ordered by RT:
  /// mapping: RT (not necessarily unique) -> reference to search hit
  typedef std::multimap<double, IdentificationData::QueryMatchRef> RTMap;
  /// mapping: charge -> internal/external: (RT -> ref. to search hit)
  typedef std::map<Int, std::pair<RTMap, RTMap>> ChargeMap;
  /// mapping: sequence (with adduct?) -> charge -> internal/external ID information
  typedef std::map<AdductedID, ChargeMap> MoleculeMap;

  /// mapping: target ion ID -> pos. in @p molecule_map_, charge
  typedef std::map<String, std::pair<MoleculeMap::iterator, Int>> TargetIonRTs;

  // need to map from a MoleculeQueryMatch to the corresponding exported
  // PeptideIdentification, so generate a look-up table:
  typedef std::tuple<double, double, String, String> PepIDKey; ///< (RT, m/z, molecule, adduct)
  typedef std::multimap<PepIDKey, const PeptideIdentification*> PepIDLookup;
  // @TODO: why does this crash when a reference is used instead of the pointer?

  /// region in RT in which a target elutes:
  struct RTRegion
  {
    double start, end;
    ChargeMap ids; ///< internal/external IDs (per charge) in this region
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
      const String& ref1 = f1.getMetaValue("CompoundRef");
      const String& ref2 = f2.getMetaValue("CompoundRef");
      if (ref1 == ref2)
      {
        return f1.getRT() < f2.getRT();
      }
      return ref1 < ref2;
    }
  } feature_compare_;

  MoleculeMap molecule_map_; ///< aggregated IDs for each identified molecule
  TargetIonRTs target_ion_rts_; ///< reference into @p molecule_map_ for each target ion ID
  PepIDLookup pep_id_lookup_; ///< mapping to PeptideIdentifications

  Size n_internal_targets_; ///< number of internal target molecules
  Size n_external_targets_; ///< number of external target molecules
  Size n_seed_targets_; ///< number of targets derived from seeds

  Size batch_size_; ///< number of target molecules to consider together during chromatogram extraction
  double rt_window_; ///< RT window width
  double mz_window_; ///< m/z window width
  bool mz_window_ppm_; ///< m/z window width is given in PPM (not Da)?

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
  void generateTransitions_(const String& target_id, double mz, Int charge,
                            const IsotopeDistribution& iso_dist);

  void addTargetRT_(TargetedExperiment::Compound& target, double rt) const;

  /// get regions in which target elutes (ideally only one) by clustering RT elution times
  void getRTRegions_(ChargeMap& charge_data, std::vector<RTRegion>& rt_regions) const;

  /// annotate identified features with m/z, isotope probabilities, etc.
  void annotateFeatures_(FeatureMap& features);

  void annotateFeaturesOneTarget_(FeatureMap& features, const String& target_id,
                                  const std::vector<Size>& indexes);

  void ensureConvexHulls_(Feature& feature);

  void postProcess_(FeatureMap& features, bool with_external_ids);

  /// some statistics on detected features
  void statistics_(const FeatureMap& features) const;

  /*!
    @brief Creates an assay library given target molecule information

    @p MoleculeMap will be (partially) cleared and thus has to be mutable.
  */
  void createAssayLibrary_(MoleculeMap::iterator begin, MoleculeMap::iterator end);

  void addTargetMolecule_(IdentificationData::QueryMatchRef ref, bool external = false);

  void checkNumObservations_(Size n_pos, Size n_neg, const String& note = "") const;

  void getUnbiasedSample_(const std::multimap<double, std::pair<Size, bool> >& valid_obs,
                          std::map<Size, Int>& training_labels);

  void getRandomSample_(std::map<Size, Int>& training_labels);

  void classifyFeatures_(FeatureMap& features);

  void filterFeaturesFinalizeAssay_(Feature& best_feature, double best_quality,
                                    const double quality_cutoff);

  void filterFeatures_(FeatureMap& features, bool classified);

  void calculateFDR_(FeatureMap& features);

  /// Look up peptide IDs based on given keys and store the results
  void lookUpPeptideIDs_(const std::set<PepIDKey> pep_id_keys,
                         std::vector<PeptideIdentification>& output);

  std::pair<String, Int> extractTargetID_(const Feature& feature,
                                          bool extract_charge = false);

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
