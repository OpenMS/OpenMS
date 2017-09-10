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
    else if (component_1.metaValueExists(feature_name))
    {
      double feature_1 = component_1.getMetaValue(feature_name);
      ratio = feature_1;
    } 

    return ratio;
  }
  
  double AbsoluteQuantitation::calculateBias(double & actual_concentration, double & calculated_concentration)
  {
    double bias = 0.0;
    bias = fabs(actual_concentration - calculated_concentration)/actual_concentration*100;
    return bias;
  }
  
  void AbsoluteQuantitation::fitCalibration(
    std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    std::string & feature_name,
    std::string & transformation_model,
    Param & transformation_model_params)
  {
    // extract out the calibration points
    TransformationModel::DataPoints data;
    for (size_t i = 0; i < component_concentrations.size(); i++){
      data.push_back(std::make_pair(component_concentrations[i].actual_concentration, component_concentrations[i].feature.getMetaValue(feature_name)));
    }

    // fit the data to the model
    AbsoluteQuantitationMethod aqm;
    Param params = aqm.fitTransformationModel(transformation_model, data, transformation_model_params);

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
    AbsoluteQuantitationMethod aqm;
    calculated_concentration = aqm.evaluateTransformationModel(
      transformation_model, ratio, transformation_model_params);

    return calculated_concentration;
  }

  void AbsoluteQuantitation::quantifyComponents(std::vector<FeatureMap>& unknowns)
  {
    //Potential Optimizations: create a map for each unknown FeatureMap
    // to reduce multiple loops

    // initialize all quant_method variables
    std::map<std::string,AbsoluteQuantitationMethod>::iterator quant_methods_it;
    std::string quant_component_name; //i.e., quant_method transition_id
    std::string quant_IS_component_name; //i.e., quant_method internal standard transition_id
    std::string quant_feature_name; //i.e., quant_method peak_apex_int or peak_area
    std::string transformation_model;
    Param transformation_model_params;

    // initialize all unknown variables
    std::string component_name; //i.e., transition_id
    std::string IS_component_name; //i.e., internal standard transition_id
    std::string component_group_name; //i.e., peptideRef
    double calculated_concentration;
    std::string concentration_units;

    // initalize all other variables
    bool IS_found;
    Feature empty_feature;
    size_t sub_it, is_sub_it; // keep sub-feature and IS sub_feature in the function scope
    FeatureMap::iterator is_feature_it; // keep the IS feature in the function scope

    // iterate through the unknowns
    for (size_t i = 0; i < unknowns.size(); i++)
    {      

      // iterate through each component_group/feature     
      for (FeatureMap::iterator feature_it = unknowns[i].begin(); feature_it != unknowns[i].end(); ++feature_it)
      {
        component_group_name = (std::string)feature_it->getMetaValue("PeptideRef");
        Feature unknowns_quant_feature;

        // iterate through each component/sub-feature
        for (sub_it = 0; sub_it < feature_it->getSubordinates().size(); ++sub_it)
        {
          component_name = (std::string)feature_it->getSubordinates()[sub_it].getMetaValue("native_id");
          quant_methods_it = quant_methods_.find(component_name);

          // apply the calibration curve to components that are in the quant_method
          if (quant_methods_it != quant_methods_.end())
          {
            quant_methods_it->second.getComponentISFeatureNames(quant_component_name,quant_IS_component_name,quant_feature_name);
            if (quant_IS_component_name != "")
            {
              // look up the internal standard for the component
              IS_found = false;
              // Optimization: 90% of the IS will be in the same component_group/feature
              for (is_sub_it = 0; is_sub_it < feature_it->getSubordinates().size(); ++is_sub_it)
              {
                IS_component_name = (std::string)feature_it->getSubordinates()[is_sub_it].getMetaValue("native_id");              
                if (quant_IS_component_name == IS_component_name)
                {
                  IS_found = true;
                  break;
                }
              }
              if (!IS_found)
              {// expand IS search to all components                
                // iterate through each component_group/feature     
                for (is_feature_it = unknowns[i].begin(); is_feature_it != unknowns[i].end(); ++is_feature_it)
                {
                  //iterate through each component/sub-feature
                  for (is_sub_it = 0; is_sub_it < is_feature_it->getSubordinates().size(); ++is_sub_it)
                  {
                    IS_component_name = (std::string)is_feature_it->getSubordinates()[is_sub_it].getMetaValue("native_id");                   
                    if (quant_IS_component_name == IS_component_name)
                    {
                      IS_found = true;
                      break;
                    }
                  }
                  if (IS_found)
                  {
                    break;
                  }
                }
              }
              if (IS_found)
              {                
                quant_methods_it->second.getTransformationModel(transformation_model,transformation_model_params);
                calculated_concentration = applyCalibration(
                  feature_it->getSubordinates()[sub_it],
                  is_feature_it->getSubordinates()[is_sub_it],
                  quant_feature_name,transformation_model,transformation_model_params);
              }
              else 
              {                
                LOG_INFO << "Component " << component_name << " IS " << quant_IS_component_name << " was not found.";
                LOG_INFO << "No concentration will be calculated.";
              }
            }
            else
            {
              quant_methods_it->second.getTransformationModel(transformation_model,transformation_model_params);
              calculated_concentration = applyCalibration(
                feature_it->getSubordinates()[sub_it],
                empty_feature,
                quant_feature_name,transformation_model,transformation_model_params);
            }

            // add new metadata (calculated_concentration, concentration_units) to the component
            feature_it->getSubordinates()[sub_it].setMetaValue("calculated_concentration",calculated_concentration);
            quant_methods_it->second.getConcentrationUnits(concentration_units);
            feature_it->getSubordinates()[sub_it].setMetaValue("concentration_units",concentration_units);
            // calculate the bias?
          }
          else 
          {
            LOG_INFO << "Component " << component_name << " does not have a quantitation method.";
            LOG_INFO << "No concentration will be calculated.";
            feature_it->getSubordinates()[sub_it].setMetaValue("calculated_concentration","");
            feature_it->getSubordinates()[sub_it].setMetaValue("concentration_units","");
          }
        }
      }
    }
  }

} // namespace

