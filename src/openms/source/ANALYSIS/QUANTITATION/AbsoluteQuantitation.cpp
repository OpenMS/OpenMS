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
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>

//Kernal classes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

//Quantitation classes
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationStandards.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>


//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>
#include <math.h>

namespace OpenMS
{
  
  AbsoluteQuantitation::AbsoluteQuantitation() :
  DefaultParamHandler("AbsoluteQuantitation")
  {
    //todo:  see MRMTransitionGroupPicker.cpp
  }
  
  AbsoluteQuantitation::~AbsoluteQuantitation()
  {
  }

  double AbsoluteQuantitation::calculateRatio(Feature & component_1, Feature & component_2, std::string & feature_name)
  {
    double ratio = 0.0;
    if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
    {
      double feature_1 = component_1.getMetaValue(feature_name);
      double feature_2 = component_2.getMetaValue(feature_name);
      ratio = feature_1/feature_2;
    } 

    return ratio;
  }
  
  double AbsoluteQuantitation::calculateBias(double & actual_concentration, double & calculated_concentration)
  {
    double bias = 0.0;
    bias = fabs(actual_concentration - calculated_concentration)/actual_concentration*100;
    return bias;
  }
  
  void AbsoluteQuantitation::fitCalibration(std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    std::string & feature_name,
    std::string & transformation_model,
    Param & transformation_model_params)
  {
    // extract out the calibration points
    TransformationModel::DataPoints data;
    for (int i = 0; i != component_concentrations.size()-1; i++){
      data.push_back(make_pair(component_concentrations[i].actual_concentration, component_concentrations[i].feature.getMetaValue(feature_name)));
    }

    // fit the data to the model
    AbsoluteQuantitationMethod aqm;
    auto tm = aqm.getTransformationModel(transformation_model);
    tm(data,transformation_model_params);

    // store the information about the fit

  }
  
  double AbsoluteQuantitation::applyCalibration(Feature & component,
    Feature & IS_component,
    std::string & feature_name,
    std::string & transformation_model,
    Param & transformation_model_params)
  {
    // calculate the ratio
    double ratio = calculateRatio(component, IS_component, feature_name);

    // calculate the absolute concentration
    double calculated_concentration = 0.0;
    TransformationModel::DataPoints empty;
    AbsoluteQuantitationMethod aqm;
    auto tm = aqm.getTransformationModel(transformation_model);
    tm(empty,transformation_model_params);
    calculated_concentration = tm.apply(ratio);

    return calculated_concentration;
  }

} // namespace

