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
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

namespace OpenMS
{


/**
    @brief Simulate SILAC experiments

    Add modified features to MS1 scans.

    @htmlinclude OpenMS_SILACLabeler.parameters
  */
  class OPENMS_DLLAPI SILACLabeler :
    public BaseLabeler
  {
public:

    /// default constructor
    SILACLabeler();

    /// destructor
    ~SILACLabeler() override;

    /// create new object (needed by Factory)
    static BaseLabeler* create()
    {
      return new SILACLabeler();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
      return "SILAC";
    }

    // redeclaration of virtual methods
    void preCheck(Param& /* param */) const override;

    void setUpHook(SimTypes::FeatureMapSimVector& /* channels */) override;

    void postDigestHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRTHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postDetectabilityHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postIonizationHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRawMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */) override;

    void postRawTandemMSHook(SimTypes::FeatureMapSimVector& /* features_to_simulate */, SimTypes::MSSimExperiment& /* simulated map */) override;

protected:
    void addModificationToPeptideHit_(Feature& feature, const String& modification) const;

    Feature mergeFeatures_(Feature& labeled_channel_feature, const String& unmodified_sequence, Map<String, Feature>& unlabeled_features_index, Int index_channel_id, Int labeled_channel_id) const;

    Feature mergeAllChannelFeatures_(Feature& heavy_channel_feature, const String& unmodified_feature_sequence, Map<String, Feature>& light_channel_feature_index, Map<String, Feature>& medium_channel_feature_index) const;

    String medium_channel_lysine_label_;
    String medium_channel_arginine_label_;

    String heavy_channel_lysine_label_;
    String heavy_channel_arginine_label_;

    bool canModificationBeApplied_(const String& modification_id, const String& aa) const;

    void applyLabelToProteinHit_(SimTypes::FeatureMapSim& protein_hit, const String& arginine_label, const String& lysine_label) const;

    void updateMembers_() override;

    String getUnmodifiedSequence_(const Feature& feature, const String& arginine_label, const String& lysine_label) const;

private:
    /// Map ID for the light/unlabeled channel
    static const int LIGHT_FEATURE_MAPID_;
    /// Map ID for the medium labeled channel
    static const int MEDIUM_FEATURE_MAPID_;
    /// Map ID for the heavy labeled channel
    static const int HEAVY_FEATURE_MAPID_;
  };

} // namespace OpenMS

