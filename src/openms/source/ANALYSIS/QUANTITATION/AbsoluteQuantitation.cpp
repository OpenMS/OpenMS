// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitation.h>

//Kernal classes
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>

#include <OpenMS/CONCEPT/LogStream.h>

//OpenSWATH classes
#include <OpenMS/ANALYSIS/OPENSWATH/MRMRTNormalizer.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

//Quantitation classes
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>

//Math classes
#include <OpenMS/MATH/StatisticFunctions.h>

//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <cmath>
#include <numeric>
#include <boost/math/special_functions/erf.hpp>
#include <algorithm>

namespace OpenMS
{

  AbsoluteQuantitation::AbsoluteQuantitation() :
  DefaultParamHandler("AbsoluteQuantitation")
  {
    defaults_.setValue("min_points", 4, "The minimum number of calibrator points.");

    defaults_.setValue("max_bias", 30.0, "The maximum percent bias of any point in the calibration curve.");

    defaults_.setValue("min_correlation_coefficient", 0.9, "The minimum correlation coefficient value of the calibration curve.");

    defaults_.setValue("max_iters", 100, "The maximum number of iterations to find an optimal set of calibration curve points and parameters.");

    defaults_.setValue("outlier_detection_method", "iter_jackknife", "Outlier detection method to find and remove bad calibration points.");
    defaults_.setValidStrings("outlier_detection_method", {"iter_jackknife","iter_residual"});

    defaults_.setValue("use_chauvenet", "true", "Whether to only remove outliers that fulfill Chauvenet's criterion for outliers (otherwise it will remove any outlier candidate regardless of the criterion).");
    defaults_.setValidStrings("use_chauvenet", {"true","false"});

    defaults_.setValue("optimization_method", "iterative", "Calibrator optimization method to find the best set of calibration points for each method.");
    defaults_.setValidStrings("optimization_method", {"iterative"});

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void AbsoluteQuantitation::updateMembers_()
  {
    min_points_ = (size_t)param_.getValue("min_points");
    max_bias_ = (double)param_.getValue("max_bias");
    min_correlation_coefficient_ = (double)param_.getValue("min_correlation_coefficient");
    max_iters_ = (size_t)param_.getValue("max_iters");
    outlier_detection_method_ = param_.getValue("outlier_detection_method").toString();
    use_chauvenet_ = (bool)param_.getValue("use_chauvenet").toBool();
    optimization_method_ = param_.getValue("optimization_method").toString();
  }

  AbsoluteQuantitation::~AbsoluteQuantitation() = default;

  void AbsoluteQuantitation::setQuantMethods(std::vector<AbsoluteQuantitationMethod>& quant_methods)
  {
    quant_methods_.clear();
    for (size_t i = 0; i < quant_methods.size(); i++)
    {
      String component_name = quant_methods[i].getComponentName();
      quant_methods_[component_name] = quant_methods[i];
    }
  }

  std::vector<AbsoluteQuantitationMethod> AbsoluteQuantitation::getQuantMethods()
  {
    std::vector<AbsoluteQuantitationMethod> quant_methods;
    for (auto const& quant_method : quant_methods_)
    {
      quant_methods.push_back(quant_method.second);
    }
    return quant_methods;
  }

  std::map<String, AbsoluteQuantitationMethod> AbsoluteQuantitation::getQuantMethodsAsMap()
  {
    return quant_methods_;
  }

  double AbsoluteQuantitation::calculateRatio(const Feature & component_1, const Feature & component_2, const String & feature_name)
  {
    double ratio = 0.0;
    // member feature_name access
    if (feature_name == "intensity")
    {
      if (component_1.metaValueExists("native_id") && component_2.metaValueExists("native_id"))
      {
        const double feature_1 = component_1.getIntensity();
        const double feature_2 = component_2.getIntensity();
        ratio = feature_1 / feature_2;
      }
      else if (component_1.metaValueExists("native_id"))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << ".";
        const double feature_1 = component_1.getIntensity();
        ratio = feature_1;
      }
    }
    // metaValue feature_name access
    else
    {
      if (component_1.metaValueExists(feature_name) && component_2.metaValueExists(feature_name))
      {
        const double feature_1 = component_1.getMetaValue(feature_name);
        const double feature_2 = component_2.getMetaValue(feature_name);
        ratio = feature_1/feature_2;
      }
      else if (component_1.metaValueExists(feature_name))
      {
        OPENMS_LOG_DEBUG << "Warning: no IS found for component " << component_1.getMetaValue("native_id") << ".";
        const double feature_1 = component_1.getMetaValue(feature_name);
        ratio = feature_1;
      }
      else
      {
        OPENMS_LOG_DEBUG << "Feature metaValue " << feature_name << " not found for components " << component_1.getMetaValue("native_id") << " and " << component_2.getMetaValue("native_id") << ".";
      }
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
    for (size_t i = 0; i < component_concentrations.size(); i++)
    {
      point.first = component_concentrations[i].actual_concentration / component_concentrations[i].IS_actual_concentration / component_concentrations[i].dilution_factor; // adjust based on the dilution factor
      double ratio = calculateRatio(component_concentrations[i].feature, component_concentrations[i].IS_feature,feature_name);
      point.second = ratio;
      data.push_back(point);
    }

    // fit the data to the model
    TransformationDescription tmd(data);
    // tmd.setDataPoints(data);
    tmd.fitModel(transformation_model, transformation_model_params);
    Param params = tmd.getModelParameters();
    // AbsoluteQuantitationMethod aqm;
    // Param params = aqm.fitTransformationModel(transformation_model, data, transformation_model_params);

    // store the information about the fit
    return params;
  }

  void AbsoluteQuantitation::calculateBiasAndR(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params,
    std::vector<double> & biases,
    double & correlation_coefficient)
  {
    // reset biases
    biases.clear();

    // extract out the calibration points
    std::vector<double> concentration_ratios, feature_amounts_ratios;
    TransformationModel::DataPoints data;
    TransformationModel::DataPoint point;
    for (size_t i = 0; i < component_concentrations.size(); ++i)
    {

      // calculate the actual and calculated concentration ratios
      double calculated_concentration_ratio = applyCalibration(component_concentrations[i].feature,
        component_concentrations[i].IS_feature,
        feature_name,
        transformation_model,
        transformation_model_params);

      double actual_concentration_ratio = component_concentrations[i].actual_concentration/
        component_concentrations[i].IS_actual_concentration / component_concentrations[i].dilution_factor;
      concentration_ratios.push_back(component_concentrations[i].actual_concentration);

      // extract out the feature amount ratios
      double feature_amount_ratio = calculateRatio(component_concentrations[i].feature,
        component_concentrations[i].IS_feature,
        feature_name);
      feature_amounts_ratios.push_back(feature_amount_ratio);

      // calculate the bias
      double bias = calculateBias(actual_concentration_ratio, calculated_concentration_ratio);
      biases.push_back(bias);

      point.first = actual_concentration_ratio;
      point.second = feature_amount_ratio;
      data.push_back(point);
    }

    // apply weighting to the feature amounts and actual concentration ratios
    TransformationModel tm(data, transformation_model_params);
    tm.weightData(data);
    std::vector<double> concentration_ratios_weighted, feature_amounts_ratios_weighted;
    for (size_t i = 0; i < data.size(); ++i)
    {
      concentration_ratios_weighted.push_back(data[i].first);
      feature_amounts_ratios_weighted.push_back(data[i].second);
    }

    // calculate the R2 (R2 = Pearson_R^2)
    correlation_coefficient = Math::pearsonCorrelationCoefficient(
      concentration_ratios_weighted.begin(), concentration_ratios_weighted.begin() + concentration_ratios_weighted.size(),
      feature_amounts_ratios_weighted.begin(), feature_amounts_ratios_weighted.begin() + feature_amounts_ratios_weighted.size()
    );
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
    TransformationModel::DataPoints data;
    TransformationDescription tmd(data);
    // tmd.setDataPoints(data);
    tmd.fitModel(transformation_model, transformation_model_params);
    tmd.invert();
    double calculated_concentration = tmd.apply(ratio);

    // AbsoluteQuantitationMethod aqm;
    // double calculated_concentration = aqm.evaluateTransformationModel(
    //   transformation_model, ratio, transformation_model_params);

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

    // initialize all other variables
    Feature empty_feature;
    size_t IS_component_it(0), IS_component_group_it(0);

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
          if (!quant_IS_component_name.empty())
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
              OPENMS_LOG_INFO << "Component " << component_name << " IS " << quant_IS_component_name << " was not found.";
              OPENMS_LOG_INFO << "No concentration will be calculated.\n";
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
          OPENMS_LOG_INFO << "Component " << component_name << " does not have a quantitation method.";
          OPENMS_LOG_INFO << "No concentration will be calculated.\n";
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("calculated_concentration","");
          unknowns[feature_it].getSubordinates()[sub_it].setMetaValue("concentration_units","");
        }
      }
    }
    // }
  }

  bool AbsoluteQuantitation::optimizeCalibrationCurveIterative(
    std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params,
    Param & optimized_params)
  {

    // sort from min to max concentration
    std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sorted = component_concentrations;
    std::sort(component_concentrations_sorted.begin(), component_concentrations_sorted.end(),
      [](const AbsoluteQuantitationStandards::featureConcentration& lhs, const AbsoluteQuantitationStandards::featureConcentration& rhs)
      {
        return lhs.actual_concentration < rhs.actual_concentration; //ascending order
      }
    );

    // indices of component_concentrations
    std::vector<size_t> component_concentrations_sorted_indices;// loop from all points to min_points
    for (size_t index = 0; index < component_concentrations_sorted.size(); ++index)
    {
      component_concentrations_sorted_indices.push_back(index);
    }

    // starting parameters
    optimized_params = transformation_model_params;

    // for (size_t n_iters = 0; n_iters < max_iters_; ++n_iters)
    for (size_t n_iters = 0; n_iters < component_concentrations_sorted.size(); ++n_iters)
    {

      // extract out components
      const std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sub = extractComponents_(
        component_concentrations_sorted, component_concentrations_sorted_indices);

      // check if the min number of calibration points has been broken
      if (component_concentrations_sorted_indices.size() < min_points_)
      {
        OPENMS_LOG_INFO << "No optimal calibration found for " << component_concentrations_sub[0].feature.getMetaValue("native_id") << " .";
        return false;  //no optimal calibration found
      }

      // fit the model
      optimized_params = fitCalibration(component_concentrations_sub,
        feature_name,
        transformation_model,
        optimized_params);

      // calculate the R2 and bias
      std::vector<double> biases; // not needed (method parameters)
      double correlation_coefficient = 0.0; // not needed (method parameters)
      calculateBiasAndR(
        component_concentrations_sub,
        feature_name,
        transformation_model,
        optimized_params,
        biases,
        correlation_coefficient);

      // check R2 and biases
      bool bias_check = true;
      for (size_t bias_it = 0; bias_it < biases.size(); ++bias_it)
      {
        if (biases[bias_it] > max_bias_)
        {
          bias_check = false;
        }
      }
      if (bias_check && correlation_coefficient > min_correlation_coefficient_)
      {
        OPENMS_LOG_INFO << "Valid calibration found for " << component_concentrations_sub[0].feature.getMetaValue("native_id") << " .";

        // copy over the final optimized points before exiting
        component_concentrations = component_concentrations_sub;
        return true;  //optimal calibration found
      }

      // R2 and biases check failed, determine potential outlier
      int pos;
      if (outlier_detection_method_ == "iter_jackknife")
      {
        // get candidate outlier: removal of which datapoint results in best rsq?
        pos = jackknifeOutlierCandidate_(
          component_concentrations_sub,
          feature_name,
          transformation_model,
          optimized_params);
      }
      else if (outlier_detection_method_ == "iter_residual")
      {
        // get candidate outlier: removal of datapoint with largest residual?
        pos = residualOutlierCandidate_(
          component_concentrations_sub,
          feature_name,
          transformation_model,
          optimized_params);
      }
      else
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          String("Method ") + outlier_detection_method_ + " is not a valid method for optimizeCalibrationCurveIterative");
      }

      // remove if residual is an outlier according to Chauvenet's criterion
      // or if testing is turned off
      if (!use_chauvenet_ || MRMRTNormalizer::chauvenet(biases, pos))
      {
        component_concentrations_sorted_indices.erase(component_concentrations_sorted_indices.begin() + pos);
      }
      else
      {
        return false;  //no optimal calibration found
      }
    }
    return false;  //no optimal calibration found
  }

  std::vector<AbsoluteQuantitationStandards::featureConcentration> AbsoluteQuantitation::extractComponents_(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
    const std::vector<size_t>& component_concentrations_indices)
  {
    std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_sub;
    for (size_t iter = 0; iter < component_concentrations_indices.size(); ++iter)
    {
      component_concentrations_sub.push_back(component_concentrations[component_concentrations_indices[iter]]);
    }
    return component_concentrations_sub;

  }

  int AbsoluteQuantitation::jackknifeOutlierCandidate_(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params)
  {
    // Returns candidate outlier: A linear regression and rsq is calculated for
    // the data points with one removed pair. The combination resulting in
    // highest rsq is considered corresponding to the outlier candidate. The
    // corresponding iterator position is then returned.
    std::vector<double> rsq_tmp;
    Param optimized_params = transformation_model_params;

    for (Size i = 0; i < component_concentrations.size(); i++)
    {
      std::vector<AbsoluteQuantitationStandards::featureConcentration> component_concentrations_tmp = component_concentrations;
      component_concentrations_tmp.erase(component_concentrations_tmp.begin() + i);

      // fit the model
      optimized_params = fitCalibration(component_concentrations_tmp,
        feature_name,
        transformation_model,
        optimized_params);

      // calculate the R2 and bias
      std::vector<double> biases;
      double correlation_coefficient = 0.0;
      calculateBiasAndR(
        component_concentrations_tmp,
        feature_name,
        transformation_model,
        optimized_params,
        biases,
        correlation_coefficient);

      rsq_tmp.push_back(correlation_coefficient);
    }
    return max_element(rsq_tmp.begin(), rsq_tmp.end()) - rsq_tmp.begin();
  }

  int AbsoluteQuantitation::residualOutlierCandidate_(
    const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
    const String & feature_name,
    const String & transformation_model,
    const Param & transformation_model_params)
  {
    // Returns candidate outlier: A linear regression and residuals are calculated for
    // the data points. The one with highest residual error is selected as the outlier candidate. The
    // corresponding iterator position is then returned.

    // fit the model
    Param optimized_params = fitCalibration(component_concentrations,
      feature_name,
      transformation_model,
      transformation_model_params);

    // calculate the R2 and bias
    std::vector<double> biases;
    double correlation_coefficient = 0.0;
    calculateBiasAndR(
      component_concentrations,
      feature_name,
      transformation_model,
      optimized_params,
      biases,
      correlation_coefficient);

    return max_element(biases.begin(), biases.end()) - biases.begin();
  }

  void AbsoluteQuantitation::optimizeCalibrationCurves(
    std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> & components_concentrations)
  {
    std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>>& cc = components_concentrations;
    for (std::pair<const String, AbsoluteQuantitationMethod>& quant_method : quant_methods_)
    {
      const String& component_name = quant_method.first;
      AbsoluteQuantitationMethod& component_aqm = quant_method.second;
      if (cc.count(component_name) && optimization_method_ == "iterative")
      {
        // optimize the calibration curve for the component
        Param optimized_params;
        bool optimal_calibration_found = optimizeCalibrationCurveIterative(
          cc[component_name],
          component_aqm.getFeatureName(),
          component_aqm.getTransformationModel(),
          component_aqm.getTransformationModelParams(),
          optimized_params);

        // order component concentrations and update the lloq and uloq
        std::vector<AbsoluteQuantitationStandards::featureConcentration>::const_iterator it;
        it = std::min_element(cc[component_name].begin(), cc[component_name].end(), [](
            const AbsoluteQuantitationStandards::featureConcentration& lhs,
            const AbsoluteQuantitationStandards::featureConcentration& rhs
          )
          {
            return lhs.actual_concentration < rhs.actual_concentration;
          }
        );
        component_aqm.setLLOQ(it->actual_concentration);
        it = std::max_element(cc[component_name].begin(), cc[component_name].end(), [](
            const AbsoluteQuantitationStandards::featureConcentration& lhs,
            const AbsoluteQuantitationStandards::featureConcentration& rhs
          )
          {
            return lhs.actual_concentration < rhs.actual_concentration;
          }
        );
        component_aqm.setULOQ(it->actual_concentration);

        if (optimal_calibration_found)
        {
          // calculate the R2 and bias
          std::vector<double> biases;
          double correlation_coefficient = 0.0;
          calculateBiasAndR(
            cc[component_name],
            component_aqm.getFeatureName(),
            component_aqm.getTransformationModel(),
            optimized_params,
            biases,
            correlation_coefficient);

          // record the updated information
          component_aqm.setCorrelationCoefficient(correlation_coefficient);
          component_aqm.setTransformationModelParams(optimized_params);
          component_aqm.setNPoints(cc[component_name].size());
        }
        else 
        {
          component_aqm.setCorrelationCoefficient(0.0);
          component_aqm.setNPoints(0);
          component_aqm.setLLOQ(0.0);
          component_aqm.setULOQ(0.0);
        }
      }
      else if (optimization_method_ != "iterative")
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          "Unsupported calibration curve optimization method '" + optimization_method_ + "'.");
      }
      else
      {
        OPENMS_LOG_DEBUG << "Warning: Standards not found for component " << component_name << ".";
      }
    }
  }

  void AbsoluteQuantitation::optimizeSingleCalibrationCurve(
    const String& component_name,
    std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations
  )
  {
    std::map<String, std::vector<AbsoluteQuantitationStandards::featureConcentration>> cc_map;
    cc_map.insert({component_name, component_concentrations});
    optimizeCalibrationCurves(cc_map);
    component_concentrations = cc_map.at(component_name);
  }
} // namespace
