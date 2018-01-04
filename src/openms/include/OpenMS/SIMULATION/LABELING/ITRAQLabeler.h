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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_ITRAQLABELER_H
#define OPENMS_SIMULATION_LABELING_ITRAQLABELER_H

#include <OpenMS/ANALYSIS/QUANTITATION/ItraqConstants.h>
#include <OpenMS/DATASTRUCTURES/Utils/MatrixUtils.h>
#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>
#include <OpenMS/SIMULATION/SimTypes.h>

namespace OpenMS
{

  /**
    @brief Simulate iTRAQ experiments

    Adds reporter ion intensities to MS/MS scans.
    Supports custom channel allocation and isotope matrices.

        @htmlinclude OpenMS_ITRAQLabeler.parameters
  */
  class OPENMS_DLLAPI ITRAQLabeler :
    public BaseLabeler
  {
public:

    typedef ItraqConstants::ChannelInfo ChannelInfo;
    typedef ItraqConstants::ChannelMapType ChannelMapType;
    typedef ItraqConstants::IsotopeMatrices IsotopeMatrices;

    /// default constructor
    ITRAQLabeler();

    /// destructor
    ~ITRAQLabeler() override;

    /// create new object (needed by Factory)
    static BaseLabeler* create()
    {
      return new ITRAQLabeler();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "itraq";
    }

    // redeclaration of virtual methods
    void preCheck(Param& param) const override;

    void setUpHook(SimTypes::FeatureMapSimVector& /* channels */) override;

    void postDigestHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRTHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postDetectabilityHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postIonizationHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRawMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRawTandemMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */, SimTypes::MSSimExperiment& /* simulated map */) override;

protected:

    /**
      @brief Modify the first peptide hit of the feature with a modification at @p pos

    */
    void addModificationToPeptideHit_(Feature& feature, const String& modification, const Size& pos) const;

    /**
      @brief tag a feature with iTRAQ modifications

      This might produce several features, due to incomplete labeling of "Y" residues.
      The resulting features are exact copies, except for their intensity and modification state.

    */
    void labelPeptide_(const Feature& feature, SimTypes::FeatureMapSim& result) const;

    Feature mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, std::map<AASequence, Feature>& unlabeled_features_index) const;

    /// Synchronize members with param class
    void updateMembers_() override;

    // get the closest RT profile factor of a feature for a given RT
    double getRTProfileIntensity_(const Feature& f, const double MS2_RT_time) const;

    /// convert meta information from feature into intensity values for iTRAQ
    EigenMatrixXdPtr getItraqIntensity_(const Feature& f, const double MS2_RT_time) const;


    // Members:

    /// set to either ItraqConstants::FOURPLEX or ItraqConstants::EIGHTPLEX
    Int itraq_type_;

    /// map the channel-name (e.g. 114) onto its description and the centroid mass
    /// the channel-name is also the id-string in the mapList section of the ConsensusMap
    ChannelMapType channel_map_;

    /// Matrices with isotope correction values (one for each plex-type)
    IsotopeMatrices isotope_corrections_;

    /// efficiency of "Y" labeling
    double y_labeling_efficiency_;

  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_ITRAQLabeler_H
