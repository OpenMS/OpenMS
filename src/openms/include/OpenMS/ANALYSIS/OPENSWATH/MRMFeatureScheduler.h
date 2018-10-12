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

#ifndef OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESCHEDULER_H
#define OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESCHEDULER_H

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <map>
#include <vector>

namespace OpenMS
{
  class OPENMS_DLLAPI MRMFeatureScheduler :
    public DefaultParamHandler
  {
public:
    MRMFeatureScheduler();
    virtual ~MRMFeatureScheduler();

    void schedule_MRMFeatures(MRMFeatureSelector& feature_selector, const FeatureMap& features, FeatureMap& output_features);
    void schedule_MRMFeaturesQMIP(const FeatureMap& features, FeatureMap& output_features);
    void schedule_MRMFeatures_score(const FeatureMap& features, FeatureMap& output_features);

    void setNNThresholds(const std::vector<double>& nn_thresholds);
    void setLocalityWeights(const std::vector<String>& locality_weights);
    void setSelectTransitionGroups(const std::vector<String>& select_transition_groups);
    void setSegmentWindowLengths(const std::vector<double>& segment_window_lengths);
    void setSegmentStepLengths(const std::vector<double>& segment_step_lengths);
    void setSelectHighestCounts(const std::vector<String>& select_highest_counts);
    void setVariableTypes(const std::vector<String>& variable_types);
    void setOptimalThresholds(const std::vector<double>& optimal_thresholds);

private:
    std::vector<double> nn_thresholds_;
    std::vector<String>   locality_weights_;
    std::vector<String>   select_transition_groups_;
    std::vector<double> segment_window_lengths_;
    std::vector<double> segment_step_lengths_;
    std::vector<String>   select_highest_counts_;
    std::vector<String> variable_types_;
    std::vector<double> optimal_thresholds_;
  };
}

#endif // OPENMS_ANALYSIS_OPENSWATH_MRMFEATURESCHEDULER_H
