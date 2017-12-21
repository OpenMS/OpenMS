// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H
#define OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H

// OpenMS
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/ANALYSIS/DECHARGING/ILPDCWrapper.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/MassExplainer.h>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

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

    typedef FeatureMap FeatureMapType;
    typedef Feature FeatureType;
    typedef DPosition<2> ClusterPointType;
    typedef FeatureMapType::FeatureType::CoordinateType CoordinateType;
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
    void compute(const FeatureMapType& fm_in, FeatureMapType& fm_out, ConsensusMap& cons_map, ConsensusMap& cons_map_p);

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
    void inferMoreEdges_(PairsType& edges, Map<Size, std::set<CmpInfo_> >& feature_adducts);

    /// A function mostly for debugging
    void printEdgesOfConnectedFeatures_(Size idx_1, Size idx_2, const PairsType& feature_relation);

    /**
      @brief returns true if the intensity filter was passed or switched off

      Filter for adding an edge only when the two features connected by it, fulfill the
      intensity criterion.
    **/
    inline bool intensityFilterPassed_(const Int q1, const Int q2, const Compomer& cmp, const FeatureType& f1, const FeatureType& f2);

    /**
      @brief determines if we should test a putative feature charge

      Answer query given the internal status of @em q_try.
      Features with q<=0 always return true.
    **/
    bool chargeTestworthy_(const Int feature_charge, const Int putative_charge, const bool other_unchanged) const;

    /// List of adducts used to explain mass differences
    MassExplainer::AdductsType potential_adducts_;
    /// labeling table
    Map<Size, String> map_label_;
    /// labeling table inverse
    Map<String, Size> map_label_inverse_;
    /// status of intensity filter for edges
    bool enable_intensity_filter_;
    /// status of charge discovery
    CHARGEMODE q_try_;
    /// amount of debug information displayed
    Int verbose_level_;

  };
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DECHARGING_FEATUREDECONVOLUTION_H
