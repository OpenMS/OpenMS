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
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

//Math classes
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <math.h>
#include <algorithm>

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

  void AbsoluteQuantitation::setQuantMethods(std::vector<AbsoluteQuantitationMethod>& quant_methods)
  {
    quant_methods_.clear();
    for (size_t i = 0; i < quant_methods.size(); i++)
    {
      String component_name = quant_methods[i].getComponentName();
      quant_methods_[component_name] = quant_methods[i];
    }
  }

  double AbsoluteQuantitation::calculateRatio(const Feature & component_1, const Feature & component_2, const String & feature_name)
  {
    double ratio = 0.0;
    if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
    {
      double feature_1 = component_1.getMetaValue(feature_name);
      double feature_2 = component_2.getMetaValue(feature_name);
      // std::cout <<  "ratio = " << ratio << "." << std::endl;
      ratio = feature_1/feature_2;
    } 
    else if (component_1.metaValueExists(feature_name))
    {
      LOG_INFO << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << ".";
      // std::cout <<  "Warning: no IS found for component " << component_1.getMetaValue("native_id") << "." << std::endl;
      double feature_1 = component_1.getMetaValue(feature_name);
      ratio = feature_1;
    } 
    else
    {
      LOG_INFO << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
      // std::cout << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << "." << std::endl;
    }

    return ratio;
  }
  
  double AbsoluteQuantitation::calculateBias(const double & actual_concentration, const double & calculated_concentration)
  {
    double bias = fabs(actual_concentration - calculated_concentration)/actual_concentration*100;
    return bias;
  }
  
  Param AbsoluteQuantitation::fitCalibration(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params)
  {
    // extract out the calibration points
    TransformationModel::DataPoints data;
    TransformationModel::DataPoint point;
    for (size_t i = 0; i < component_concentrations.size(); i++){
      point.first = component_concentrations[i].actual_concentration;
      point.second = component_concentrations[i].feature.getMetaValue(feature_name);
      data.push_back(point);
    }

    // fit the data to the model
    AbsoluteQuantitationMethod aqm;
    Param params = aqm.fitTransformationModel(transformation_model, data, transformation_model_params);

    // store the information about the fit
    return params;
  }
  
  double AbsoluteQuantitation::applyCalibration(const Feature & component,
    const Feature & IS_component,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params)
  {
    // calculate the ratio
    double ratio = calculateRatio(component, IS_component, feature_name);

    // calculate the absolute concentration
    AbsoluteQuantitationMethod aqm;
    double calculated_concentration = aqm.evaluateTransformationModel(
      transformation_model, ratio, transformation_model_params);

    // check for less than zero
    if (calculated_concentration < 0.0)
    {
      calculated_concentration = 0.0;
    }

    return calculated_concentration;
  }

  void AbsoluteQuantitation::quantifyComponents(FeatureMap& unknowns)
  {
    //Potential Optimizations: create a map for each unknown FeatureMap
    // to reduce multiple loops

    // initalize all other variables
    Feature empty_feature;
    size_t IS_component_it, IS_component_group_it;

    // // iterate through the unknowns
    // for (size_t i = 0; i < unknowns.size(); i++)
    // {      

    // iterate through each component_group/feature     
    for (size_t feature_it = 0; feature_it < unknowns.size(); ++feature_it)
    {
      String component_group_name = (String)unknowns[feature_it].getMetaValue("PeptideRef");
      Feature unknowns_quant_feature;

      // iterate through each component/sub-feature
      for (size_t sub_it = 0; sub_it < unknowns[feature_it].getSubordinates().size(); ++sub_it)
      {
        String component_name = (String)unknowns[feature_it].getSubordinates()[sub_it].getMetaValue("native_id");  

        // apply the calibration curve to components that are in the quant_method
        if (quant_methods_.count(component_name)>0)
        { 
          double calculated_concentration = 0.0;    
          std::map<String,AbsoluteQuantitationMethod>::iterator quant_methods_it = quant_methods_.find(component_name);
          String quant_component_name = quant_methods_it->second.getComponentName();
          String quant_IS_component_name = quant_methods_it->second.getISName();
          String quant_feature_name = quant_methods_it->second.getFeatureName();
          if (quant_IS_component_name != "")
          {
            // look up the internal standard for the component
            bool IS_found = false;
            // Optimization: 90% of the IS will be in the same component_group/feature
            for (size_t is_sub_it = 0; is_sub_it < unknowns[feature_it].getSubordinates().size(); ++is_sub_it)
            {
              String IS_component_name = (String)unknowns[feature_it].getSubordinates()[is_sub_it].getMetaValue("native_id");              
              if (quant_IS_component_name == IS_component_name)
              {
                IS_found = true;
                IS_component_group_it = feature_it;
                IS_component_it = is_sub_it;
                break;
              }
            }
            if (!IS_found)
            {// expand IS search to all components                
              // iterate through each component_group/feature     
              for (size_t is_feature_it = 0; is_feature_it < unknowns.size(); ++is_feature_it)
              {
                //iterate through each component/sub-feature
                for (size_t is_sub_it = 0; is_sub_it < unknowns[is_feature_it].getSubordinates().size(); ++is_sub_it)
                {
                  String IS_component_name = (String)unknowns[is_feature_it].getSubordinates()[is_sub_it].getMetaValue("native_id");                   
                  if (quant_IS_component_name == IS_component_name)
                  {
                    IS_found = true;
                    IS_component_group_it = is_feature_it;
                    IS_component_it = is_sub_it;
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
              String transformation_model = quant_methods_it->second.getTransformationModel();
              Param transformation_model_params = quant_methods_it->second.getTransformationModelParams();
              calculated_concentration = applyCalibration(
                unknowns[feature_it].getSubordinates()[sub_it],
                unknowns[IS_component_group_it].getSubordinates()[IS_component_it],
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
            String transformation_model = quant_methods_it->second.getTransformationModel();
            Param transformation_model_params = quant_methods_it->second.getTransformationModelParams();
            calculated_concentration = applyCalibration(
              unknowns[feature_it].getSubordinates()[sub_it],
              empty_feature,
              quant_feature_name,transformation_model,transformation_model_params);
          }

          // add new metadata (calculated_concentration, concentration_units) to the component
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("calculated_concentration",calculated_concentration);
          String concentration_units = quant_methods_it->second.getConcentrationUnits();
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("concentration_units",concentration_units);
          // calculate the bias?
        }
        else 
        {
          LOG_INFO << "Component " << component_name << " does not have a quantitation method.";
          LOG_INFO << "No concentration will be calculated.";
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("calculated_concentration","");
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("concentration_units","");
        }
      }
    }
    // }
  }

  void AbsoluteQuantitation::findIS_()
  {
    //TODO: possible refactor the method to include a seperate function to find the IS
  }

  /** TODO: this method is incomplete
   * 1. interface with MRMRTNormalizer
   * 2. make tests
   */
  void AbsoluteQuantitation::optimizeCalibrationCurveBruteForce(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params,
    Param & optimized_params)
  {
    
    std::vector<double> biases;
    double r2 = 0.0;
    bool bias_check;

    //TODO use internal params
    size_t min_points = 4;
    double max_bias = 30.0;
    double min_r2 = 0.9; 
    // size_t max_outliers = 1;  // not used currently

    std::vector<AbsoluteQuantitationStandards::featureConcentration>::const_iterator component_start_it;
    std::vector<AbsoluteQuantitationStandards::featureConcentration>::const_iterator component_end_it;

    // sort from min to max concentration
    std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sorted = component_concentrations;
    std::sort(component_concentrations_sorted.begin(), component_concentrations_sorted.end(),
      [](AbsoluteQuantitationStandards::featureConcentration lhs, AbsoluteQuantitationStandards::featureConcentration rhs)
      {
        return lhs.actual_concentration < rhs.actual_concentration;    
      }
    );

    // loop from all points to min_points
    for (size_t n_points = component_concentrations.size(); n_points >= min_points; --n_points)
    {
      size_t n_loops = component_concentrations.size() - n_points;
      for (size_t  component_it = 0; component_it < n_loops; ++component_it)
      {
        // TODO: support for outliers
        // // loop through max_outliers
        // for (size_t n_outliers = 0; n_outliers <= max_outliers; ++n_outliers)
        // {

          // extract out components
          component_start_it = component_concentrations.begin() + component_it;
          component_end_it = component_concentrations.begin() + component_it + n_points;
          const std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sub(component_start_it, component_end_it);

          // fit the model
          optimized_params = fitCalibration(component_concentrations_sub,
            feature_name,
            transformation_model,
            transformation_model_params);

          // calculate the R2 and bias
          calculateBiasAndR2(
            component_concentrations,
            feature_name,
            transformation_model,
            transformation_model_params,
            biases,
            r2);
          
          // check R2 and biases
          bias_check = true;
          for (size_t bias_it = 0; bias_it != biases.size(); --bias_it)
          {
            if (biases[bias_it] > max_bias)
            {
              bias_check = false;
            }
          }
          if (bias_check && r2 > min_r2)
          {
            LOG_INFO << "Valid calibration found for " << component_concentrations[0].feature.getMetaValue("native_id") << " .";
            return;
          }   
        // }   
      }
    }
    LOG_INFO << "Valid calibration not found for " << component_concentrations[0].feature.getMetaValue("native_id") << " .";
  }
  
  void AbsoluteQuantitation::calculateBiasAndR2(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params,
    std::vector<double> & biases,
    double & r2_value)
  {
    
    // extract out the calibration points
    std::vector<double> concentration_ratios, feature_amounts_ratios;
    for (size_t i = 0; i < component_concentrations.size(); ++i)
    {

      // calculate the actual and calculated concentration ratios
      double calculated_concentration_ratio = applyCalibration(component_concentrations[i].feature,
        component_concentrations[i].IS_feature,
        feature_name,
        transformation_model,
        transformation_model_params);

      double actual_concentration_ratio = component_concentrations[i].actual_concentration/component_concentrations[i].IS_actual_concentration;
      concentration_ratios.push_back(component_concentrations[i].actual_concentration);

      // calculate the bias
      double bias = calculateBias(actual_concentration_ratio, calculated_concentration_ratio);
      biases.push_back(bias);

      // extract out the feature amount ratios
      double feature_amount_ratio = calculateRatio(component_concentrations[i].feature,
        component_concentrations[i].IS_feature,
        feature_name);
      feature_amounts_ratios.push_back(feature_amount_ratio);
      
      //DEBUG
      // std::cout << "calculated_concentration_ratio = " << calculated_concentration_ratio << "." << std::endl;
      // std::cout << "actual_concentration_ratio = " << actual_concentration_ratio << "." << std::endl;
      // std::cout << "bias = " << bias << "." << std::endl;
      // std::cout << "feature_amount = " << (String)component_concentrations[i].feature.getMetaValue(feature_name) << "." << std::endl;
      // std::cout << "IS_feature_amount = " << (String)component_concentrations[i].IS_feature.getMetaValue(feature_name) << "." << std::endl;
      // std::cout << "feature_amount_ratio = " << bias << "." << std::endl;
    }

    // calculate the R2 (R2 = Pearson_R^2)
    double r_value = Math::pearsonCorrelationCoefficient(
      concentration_ratios.begin(), concentration_ratios.begin() + concentration_ratios.size(),
      feature_amounts_ratios.begin(), feature_amounts_ratios.begin() + feature_amounts_ratios.size()
    ); 
    r2_value = r_value*r_value;

    //DEBUG
    // std::cout << "r_value = " << r_value << "." << std::endl;
    // std::cout << "r2_value = " << r2_value << "." << std::endl;

  }

} // namespace

