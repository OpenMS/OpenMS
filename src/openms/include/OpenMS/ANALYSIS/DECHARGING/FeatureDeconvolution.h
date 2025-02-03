// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#pragma once

// OpenMS
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <map>

namespace OpenMS
{

  class Compomer;

  /**
    @brief An algorithm to decharge features (i.e. as found by FeatureFinder).

    @htmlinclude OpenMS_FeatureDeconvolution.parameters

    @ingroup Analysis
  */

  class OPENMS_DLLAPI FeatureDeconvolution :
    public DefaultParamHandler
  {
public:

    enum CHARGEMODE {QFROMFEATURE = 1, QHEURISTIC, QALL};

    typedef DPosition<2> ClusterPointType;
    typedef Feature::CoordinateType CoordinateType;
    typedef ILPDCWrapper::PairsType PairsType;

    /** @name Constructors and Destructor s
    */
    //@{
    /// default constructor
    FeatureDeconvolution();

    /// Copy constructor
    FeatureDeconvolution(const FeatureDeconvolution& source);

    /// Assignment operator
    FeatureDeconvolution& operator=(const FeatureDeconvolution& source);

    /// destructor
    ~FeatureDeconvolution() override;
    //@}

    /**
      @brief Compute a zero-charge feature map from a set of charged features

      Find putative ChargePairs, then score them and hand over to ILP.

      @param fm_in  Input feature-map
      @param fm_out Output feature-map (sorted by position and augmented with user params)
      @param cons_map   [out] Output of grouped features belonging to a charge group
      @param cons_map_p [out] Output of paired features connected by an edge
    */
    void compute(const FeatureMap& fm_in, FeatureMap& fm_out, ConsensusMap& cons_map, ConsensusMap& cons_map_p);

protected:

    void updateMembers_() override;

    /**
      @brief 1-sided Compomer for a feature

      Holds information on an explicit (with H+) 1-sided Compomer of a feature.
    */
    struct CmpInfo_;

    /*
      @brief test for obviously wrong parameter settings and warn

      Currently supports the following scenarios:
      * If the lower charge bound is too high, consensus features with gapped, even charges will occur (e.g. 30,32,34), instead of the true (15,16,17)
      When more than 5% of the cf's look like this, we report it.

    */
    void checkSolution_(const ConsensusMap& cons_map) const;

    /// test if "simple" edges have alternative
    /// (more difficult explanation) supported by neighboring edges
    /// e.g. (.)   -> (H+) might be augmented to
    ///      (Na+) -> (H+Na+)
    void inferMoreEdges_(PairsType& edges, std::map<Size, std::set<CmpInfo_> >& feature_adducts);

    /// A function mostly for debugging
    void printEdgesOfConnectedFeatures_(Size idx_1, Size idx_2, const PairsType& feature_relation);

    /**
      @brief returns true if the intensity filter was passed or switched off

      Filter for adding an edge only when the two features connected by it, fulfill the
      intensity criterion.
    **/
    inline bool intensityFilterPassed_(const Int q1, const Int q2, const Compomer& cmp, const Feature& f1, const Feature& f2) const;

    /**
      @brief determines if we should test a putative feature charge

      Answer query given the internal status of @em q_try.
      Features with q<=0 always return true.
    **/
    bool chargeTestworthy_(const Int feature_charge, const Int putative_charge, const bool other_unchanged) const;

    /// List of adducts used to explain mass differences
    MassExplainer::AdductsType potential_adducts_;
    /// labeling table
    std::map<Size, String> map_label_;
    /// labeling table inverse
    std::map<String, Size> map_label_inverse_;
    /// status of intensity filter for edges
    bool enable_intensity_filter_;
    /// status of charge discovery
    CHARGEMODE q_try_;
    /// amount of debug information displayed
    Int verbose_level_;

  };
} // namespace OpenMS
