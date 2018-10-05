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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESELECTOR_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESELECTOR_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>

namespace OpenMS
{
  class OPENMS_DLLAPI MRMFeatureSelector :
    public DefaultParamHandler
  {
    friend class MRMFeatureSelectorQMIP;
    friend class MRMFeatureSelectorScore;
public:
    MRMFeatureSelector();
    virtual ~MRMFeatureSelector();

    virtual std::vector<String> optimize(
      std::vector<std::pair<double, String>> time_to_name, 
      std::map< String, std::vector<Feature> > feature_name_map
    )=0;
    FeatureMap select_MRMFeature(FeatureMap& features);
    virtual double make_score(Feature& feature)=0;

    String remove_spaces(String str);

    void setNNThreshold(const double& nn_threshold);
    double getNNThreshold() const;

    void setLocalityWeight(const bool& locality_weight);
    bool getLocalityWeight() const;

    void setSelectTransitionGroup(const bool& select_transition_group);
    bool getSelectTransitionGroup() const;

    void setSegmentWindowLength(const double& segment_window_length);
    double getSegmentWindowLength() const;

    void setSegmentStepLength(const double& segment_step_length);
    double getSegmentStepLength() const;

    void setSelectHighestCount(const bool& select_highest_count);
    bool getSelectHighestCount() const;

    void setVariableType(const String& variable_type);
    String getVariableType() const;

    void setOptimalThreshold(const double& optimal_threshold);
    double getOptimalThreshold() const;

    void getDefaultParameters(Param& params);

protected:
    void updateMembers_(); /// overridden function from DefaultParamHandler to keep members up to date, when a parameter is changed

private:
    double nn_threshold_;
    bool   locality_weight_;
    bool   select_transition_group_;
    double segment_window_length_;
    double segment_step_length_;
    bool   select_highest_count_;
    String variable_type_;
    double optimal_threshold_;
    Int _addVariable(LPWrapper& problem, String& name, bool bounded, double obj);
    void _addConstraint(LPWrapper& problem, size_t size, Int *indices_array, double *values_array, String name, double lb, double ub, LPWrapper::Type param);
  };

  class OPENMS_DLLAPI MRMFeatureSelectorQMIP : public MRMFeatureSelector
  {
public:
    std::vector<String> optimize(
      std::vector<std::pair<double, String>> time_to_name, 
      std::map< String, std::vector<Feature> > feature_name_map
    );
    double make_score(Feature& feature);
  };

  class OPENMS_DLLAPI MRMFeatureSelectorScore : public MRMFeatureSelector
  {
public:
    std::vector<String> optimize(
      std::vector<std::pair<double, String>> time_to_name, 
      std::map< String, std::vector<Feature> > feature_name_map
    );
    double make_score(Feature& feature);
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESELECTOR_H
