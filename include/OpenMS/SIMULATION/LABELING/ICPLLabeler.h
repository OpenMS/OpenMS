// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#ifndef OPENMS_SIMULATION_LABELING_ICPLLABELER_H
#define OPENMS_SIMULATION_LABELING_ICPLLABELER_H

#include <OpenMS/SIMULATION/LABELING/BaseLabeler.h>

namespace OpenMS
{

  /**
    @brief Simulate ICPL experiments

    Add modified features to MS1 scans. Highly Experimental.

    @htmlinclude OpenMS_ICPLLabeler.parameters
  */
  class OPENMS_DLLAPI ICPLLabeler
    : public BaseLabeler
  {
  public:

    /// default constructor
    ICPLLabeler();

    /// destructor
    virtual ~ICPLLabeler();

    /// create new object (needed by Factory)
    static BaseLabeler* create()
    {
        return new ICPLLabeler();
    }

    /// name of the model (needed by Factory)
    static const String getProductName()
    {
        return "ICPL";
    }

    // redeclaration of virtual methods
    void preCheck(Param &param) const;

    void setUpHook(FeatureMapSimVector & /* channels */);

    void postDigestHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRTHook(FeatureMapSimVector & /* features_to_simulate */);

    void postDetectabilityHook(FeatureMapSimVector & /* features_to_simulate */);

    void postIonizationHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRawMSHook(FeatureMapSimVector & /* features_to_simulate */);

    void postRawTandemMSHook(FeatureMapSimVector & /* features_to_simulate */, MSSimExperiment & /* simulated map */);

  protected:
    void addModificationToPeptideHit_(Feature& feature, const String& modification) const;

    void addLabelToProteinHits_(FeatureMapSim& features, const String& label) const;

    Feature mergeFeatures_(Feature& labeled_channel_feature, const AASequence& unmodified_sequence, Map<String, Feature>& unlabeled_features_index) const;

    String light_channel_label_;
    String medium_channel_label_;
    String heavy_channel_label_;

    void updateMembers_();

    String getUnmodifiedAASequence_(const Feature& feature, const String& label) const;
  };
} // namespace OpenMS

#endif //#ifndef OPENMS_SIMULATION_LABELING_ICPLLABELER_H
