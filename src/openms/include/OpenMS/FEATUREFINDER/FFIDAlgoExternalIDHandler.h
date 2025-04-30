// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/Types.h>

#include <map>
#include <set>
#include <vector>

namespace OpenMS
{
  /**
   * @brief Class for handling external peptide identifications in feature finding
   *
   * This class encapsulates all functionality related to external peptide IDs in the
   * feature finding process, including storage, RT transformation, and feature annotation.
   */
  class OPENMS_DLLAPI FFIDAlgoExternalIDHandler
  {
  public:
    /// RTMap for external data structure storage
    typedef std::multimap<double, PeptideIdentification*> ExternalRTMap;
    
    /// Charge to External RTMap mapping
    typedef std::map<Int, ExternalRTMap> ExternalChargeMap;
    
    /// Sequence to External Charge Map mapping
    typedef std::map<AASequence, ExternalChargeMap> ExternalPeptideMap;

    /// Default constructor
    FFIDAlgoExternalIDHandler();

    /// Reset the handler's state
    void reset();

    /// Add an external peptide to the handler's map
    void addExternalPeptide(PeptideIdentification& peptide);
    
    /// Process external peptide IDs
    void processExternalPeptides(std::vector<PeptideIdentification>& peptides_ext);
    
    /// Align internal and external IDs to estimate RT shifts and return RT uncertainty
    double alignInternalAndExternalIDs(
        const std::vector<PeptideIdentification>& peptides_internal,
        const std::vector<PeptideIdentification>& peptides_external,
        double rt_quantile);
    
    /// Transform RT from internal to external scale
    double transformRT(double rt) const;
    
    /// Check if we have RT transformation data
    bool hasRTTransformation() const;
    
    /// Get the RT transformation
    const TransformationDescription& getRTTransformation() const;
        
    /// Classify features using SVM
    void classifyFeaturesWithSVM(FeatureMap& features, const Param& param);

    /// Filter classified features
    void filterClassifiedFeatures(FeatureMap& features, double quality_cutoff);

    /// Calculate FDR for classified features
    void calculateFDR(FeatureMap& features);
    
    /// Get SVM probabilities for internal features
    const std::map<double, std::pair<Size, Size> >& getSVMProbsInternal() const;
     
  private:
    /// Add external peptide to charge map (merged version for compatibility)
    void addExternalPeptideToMap_(PeptideIdentification& peptide,
                               std::map<AASequence,
                               std::map<Int, std::pair<std::multimap<double, PeptideIdentification*>,
                                                      std::multimap<double, PeptideIdentification*>>>>& peptide_map);
    
    /// Fill an external RTMap from our data for a specific peptide and charge
    bool fillExternalRTMap_(const AASequence& sequence, Int charge,
                         std::multimap<double, PeptideIdentification*>& rt_map);
    
    /// Check and set feature class based on external data
    void annotateFeatureWithExternalIDs_(Feature& feature);
  
    /// Initialize SVM parameters
    void initSVMParameters_(const Param& param);

    /// Finalize assay features
    void finalizeAssayFeatures_(Feature& best_feature, double best_quality, double quality_cutoff);

    /// Get random sample for SVM training
    void getRandomSample_(std::map<Size, double>& training_labels);

    /// Check observation counts for SVM
    void checkNumObservations_(Size n_pos, Size n_neg, const String& note = "") const;

    /// Get unbiased sample for SVM training
    void getUnbiasedSample_(const std::multimap<double, std::pair<Size, bool> >& valid_obs,
                          std::map<Size, double>& training_labels);

    /// Add dummy peptide identification from external data
    void addDummyPeptideID_(Feature& feature, const PeptideIdentification* ext_id);
    
    /// Handle external feature probability
    void handleExternalFeature_(Feature& feature, double prob_positive, double quality_cutoff);
    
    /// Adjust FDR calculation for external features
    void adjustFDRForExternalFeatures_(std::vector<double>& fdr_probs,
                                    std::vector<double>& fdr_qvalues,
                                    Size n_internal_features);

    /// External peptide storage
    ExternalPeptideMap external_peptide_map_;
    
    /// RT transformation description
    TransformationDescription rt_transformation_;
    
    /// Number of external peptides
    Size n_external_peptides_;
    
    /// Number of external features
    Size n_external_features_;
    
    /// SVM probabilities for external features
    std::multiset<double> svm_probs_external_;

    /// SVM probabilities for internal features
    std::map<double, std::pair<Size, Size> > svm_probs_internal_;
    
    /// SVM number of parts for cross-validation
    Size svm_n_parts_;
    
    /// SVM number of samples for training
    Size svm_n_samples_;
    
    /// SVM minimum probability threshold
    double svm_min_prob_;

    /// SVM quality cutoff
    double svm_quality_cutoff;

    /// SVM predictor names
    std::vector<String> svm_predictor_names_;

    /// SVM cross-validation output file
    String svm_xval_out_;

    /// Debug level
    Int debug_level_;

    /// Number of internal features
    Size n_internal_features_;
  };

} // namespace OpenMS
