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

#ifndef OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATION_H
#define OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATION_H

#include <OpenMS/config.h>

//Kernal classes
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MRMFeature.h>

//Analysis classes
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationDescription.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModel.h>
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

//Quantitation classes
#include <OpenMS/METADATA/AbsoluteQuantitationStandards.h>
#include <OpenMS/ANALYSIS/QUANTITATION/AbsoluteQuantitationMethod.h>


//Standard library
#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  /**
    @brief AbsoluteQuantitation is a class to support absolute or relative quantitation for targeted or untargeted
    quantitation workflows (e.g., Isotope Dilution Mass Spectrometry).

    Method:
    A transformation model where y = ratio (analyte/IS) corresponding to peak height or peak area and
      x = ratio (analyte/IS) corresponding to concentration is used to fit a series of runs
      with standards of known concentrations that span the detection range of the instrument.  
      The fitted transformation model can then be used to quantify the concentration of an analyte
      in an unknown sample given the analyte peak height or area, IS peak height or area, and 
      IS concentration.

    Terms:
    component: A protein, peptide, or compund fragment, transition, or whole species that is measured by e.g.,
      LC-MS, LC-MS/MS, GC-MS, GC-MS/MS, LC-MS-TOF, HPLC-UV, HPLC-IR, etc.
    calibration curve:  A series of standards that are used to correlate instrument measurements to
      actual concentrations
  */
  class OPENMS_DLLAPI AbsoluteQuantitation :
    public DefaultParamHandler
  {

public:    
    //@{
    /// Constructor
    AbsoluteQuantitation();

    /// Destructor
    ~AbsoluteQuantitation();
    //@}
 
    /**
      @brief quant_method setter.  A list of AbsoluteQuantitationMethod classes are given as input
        and a map is constructed based on their component_name member.

      @param quant_methods A list of AbsoluteQuantitationMethod classes
    */ 
    void setQuantMethods(std::vector<AbsoluteQuantitationMethod>& quant_methods);

 
    /**
      @brief quant_method getter.  A list of AbsoluteQuantitationMethod classes are returned.
    */ 
    std::vector<AbsoluteQuantitationMethod> getQuantMethods();
    std::map<String, AbsoluteQuantitationMethod> getQuantMethodsAsMap();
 
    /**
      @brief This function calculates the ratio between features.

      @param component_1 component of the numerator
      @param component_2 component of the denomenator
      @param feature_name name of the feature to calculate the ratio on
       e.g., peak_apex, peak_area

      @return The ratio.

      @exception Exception::UnableToFit
    */ 
    double calculateRatio(const Feature & component_1, const Feature & component_2, const String & feature_name);
                   
    /**
      @brief This function calculates the bias of the calibration.

      The bias is defined as the following:
        |actual_concentration - calculated_concentration|/actual_concentration * 100%
      This is in contrast to accuracy, which is defined as the following:
        calculated_concentration/actual_concentration * 100%

      @param actual_concentration the actual concentration of the component
      @param calculated_concentration the calibration curve back calculated concentration 
        of the component

      @return The bias.

      @exception Exception::UnableToFit
    */ 
    double calculateBias(const double & actual_concentration, const double & calculated_concentration);
            
    /**
      @brief This function fits the calibration points to the model.

      @param component_concentrations list of structures with features and concentrations
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model

      @returns updated Param object

      @exception Exception::UnableToFit
    */ 
    Param fitCalibration(const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params);
      
    /**
      @brief This function calculates the biases and the correlation coefficient of the calibration points.

      @param component_concentrations list of structures with features and concentrations
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model
      @param biases Vector of point biases
      @param correlation_coefficient Pearson's R

      @exception None
    */ 
    void calculateBiasAndR(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params,
      std::vector<double> & biases,
      double & correlation_coefficient);
      
    /**
      @brief This function optimizes the parameters of the calibration for a 
        given component iteratively.

      @param component_concentrations list of structures with features and concentrations.  
        The optimal points will be returned.
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model
      @param optimized_params optimized parameters

      @exception Exception::UnableToFit
    */ 
    void optimizeCalibrationCurveIterative(
      std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params,
      Param & optimized_params);
        
    /**
      @brief This function optimizes the parameters of the calibration for a 
        all components.

      @param components_concentrations An AbsoluteQuantitationStandards::components_to_concentrations type.
        Note that the method will update the list of featureConcentrations in place.  The resulting
        components_concentrations will reflect the optimal set of points for downstream QC/QA.

    */ 
    void optimizeCalibrationCurves(std::map<String,std::vector<AbsoluteQuantitationStandards::featureConcentration>> & components_concentrations);    

    /**
      @brief This function applies the calibration curve to the component.

      @param component the component to be quantified
      @param IS_component the internal standard (IS) of the component to be quantified.
        This can be an empty feature if there is no IS for the component.
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model

      @return The absolute concentration.

      @exception Exception::UnableToFit
    */ 
    double applyCalibration(const Feature & component,
      const Feature & IS_component,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params);     
      
    /**
      @brief This function applies the calibration curve to all components.

      An additional annotation for metaValue of "calculated_concentration" and "concentration_units"
        corresponding to the absolute concentration as back-calculated from the
        fitted calibration curve model and parameters will be added to each sub-feature. 
        It is assumed that all duplicate components have been removed.  If not,
        the function will quantify all components, but the first internal standard found will
        be used to calculate the ratio for the calculation.

      @param unknowns A FeatureMap to quantify.

    */ 
    void quantifyComponents(FeatureMap& unknowns);    

protected:
    /**
      @brief This function extractous out the components.

      @param component_concentrations list of structures with features and concentrations
      @param component_concentrations_indices indices to extract out

      @returns component_concentrations_sub sublist of structures with features and concentrations.

      @exception None
    */ 
    std::vector<AbsoluteQuantitationStandards::featureConcentration> extractComponents_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration> & component_concentrations,
      const std::vector<size_t>& component_concentrations_indices);
  
    /**
      @brief This function computes a candidate outlier point by iteratively
       leaving one point out to find the one which results in the maximum R^2
       of a first order linear regression of the remaining ones.

      @param component_concentrations list of structures with features and concentrations
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model

      @return The position of the candidate outlier point in component_concentrations.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    int jackknifeOutlierCandidate_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params);

    /**
      @brief This function computes a candidate outlier point by computing
       the residuals of all points to the linear fit and selecting the one with
       the largest deviation.

      @param component_concentrations list of structures with features and concentrations
      @param feature_name name of the feature to calculate the absolute concentration.
      @param transformation_model model used to fit the calibration points
      @param transformation_model_params parameters used by the transformation_model

      @return The position of the candidate outlier point in component_concentrations.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    int residualOutlierCandidate_(
      const std::vector<AbsoluteQuantitationStandards::featureConcentration>& component_concentrations,
      const String & feature_name,
      const String & transformation_model,
      const Param & transformation_model_params);
     
private:  
    /// Synchronize members with param class
    void updateMembers_();
    
    size_t min_points_;
    double max_bias_;
    double min_correlation_coefficient_; 
    size_t max_iters_;
    String outlier_detection_method_;
    bool use_chauvenet_;
    String optimization_method_;
    
    // members
    /// map between components and quantitation methods
    std::map<String, AbsoluteQuantitationMethod> quant_methods_;

  };

}
#endif // OPENMS_ANALYSIS_QUANTITATION_ABSOLUTEQUANTITATION_H

