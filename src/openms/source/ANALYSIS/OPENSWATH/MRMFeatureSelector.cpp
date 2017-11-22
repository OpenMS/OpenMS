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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>

namespace OpenMS
{
  MRMFeatureSelector::MRMFeatureSelector() :
    DefaultParamHandler("MRMFeatureSelector")
  {
    getDefaultParameters(defaults_);
    defaultsToParam_(); // write defaults into Param object param_
  }

  MRMFeatureSelector::~MRMFeatureSelector() {}

  void MRMFeatureSelector::optimize_Tr() {}

  void MRMFeatureSelector::optimize_score() {}

  void MRMFeatureSelector::select_MRMFeature_qmip(
    FeatureMap& features,
    TargetedExperiment& targeted_exp
  )
  {}

  void MRMFeatureSelector::select_MRMFeature_score() {}

  double MRMFeatureSelector::make_score(
    Feature& feature,
    String& metaValue,
    double& weight // maybe others
  )
  {}

  void MRMFeatureSelector::setNNThreshold(const double& nn_threshold)
  {
    nn_threshold_ = nn_threshold;
  }

  double MRMFeatureSelector::getNNThreshold() const
  {
    return nn_threshold_;
  }

  void MRMFeatureSelector::setLocalityWeight(const bool& locality_weight)
  {
    locality_weight_ = locality_weight;
  }

  bool MRMFeatureSelector::getLocalityWeight() const
  {
    return locality_weight_;
  }

  void MRMFeatureSelector::setSelectTransitionGroup(const bool& select_transition_group)
  {
    select_transition_group_ = select_transition_group;
  }

  bool MRMFeatureSelector::getSelectTransitionGroup() const
  {
    return select_transition_group_;
  }

  void MRMFeatureSelector::setSegmentWindowLength(const double& segment_window_length)
  {
    segment_window_length_ = segment_window_length;
  }

  double MRMFeatureSelector::getSegmentWindowLength() const
  {
    return segment_window_length_;
  }

  void MRMFeatureSelector::setSegmentStepLength(const double& segment_step_length)
  {
    segment_step_length_ = segment_step_length;
  }

  double MRMFeatureSelector::getSegmentStepLength() const
  {
    return segment_step_length_;
  }

  void MRMFeatureSelector::setSelectHighestCount(const bool& select_highest_count)
  {
    select_highest_count_ = select_highest_count;
  }

  bool MRMFeatureSelector::getSelectHighestCount() const
  {
    return select_highest_count_;
  }

  void MRMFeatureSelector::setVariableType(const String& variable_type)
  {
    variable_type_ = variable_type;
  }

  String MRMFeatureSelector::getVariableType() const
  {
    return variable_type_;
  }

  void MRMFeatureSelector::setOptimalThreshold(const double& optimal_threshold)
  {
    optimal_threshold_ = optimal_threshold;
  }

  double MRMFeatureSelector::getOptimalThreshold() const
  {
    return optimal_threshold_;
  }

  void MRMFeatureSelector::getDefaultParameters(Param& params)
  {
    params.clear();
    // TODO Adjust defaults
    // TODO set limits on parameters
    params.setValue("nn_threshold", 4.0);
    params.setValue("locality_weight", "false");
    params.setValue("select_transition_group", "true");
    params.setValue("segment_window_length", 8.0);
    params.setValue("segment_step_length", 4.0);
    params.setValue("select_highest_count", "false");
    params.setValue("variable_type", "continuous");
    params.setValue("optimal_threshold", 0.5);
  }

  void MRMFeatureSelector::updateMembers_()
  {
    nn_threshold_ = (double)param_.getValue("nn_threshold");
    locality_weight_ = param_.getValue("locality_weight").toBool();
    select_transition_group_ = param_.getValue("select_transition_group").toBool();
    segment_window_length_ = (double)param_.getValue("segment_window_length");
    segment_step_length_ = (double)param_.getValue("segment_step_length");
    select_highest_count_ = param_.getValue("select_highest_count").toBool();
    variable_type_ = (String)param_.getValue("variable_type");
    optimal_threshold_ = (double)param_.getValue("optimal_threshold");
  }
}
