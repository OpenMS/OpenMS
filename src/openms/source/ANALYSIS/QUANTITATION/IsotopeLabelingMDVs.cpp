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
// $Maintainer: Ahmed Khalil $
// $Authors: Ahmed Khalil $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/QUANTITATION/IsotopeLabelingMDVs.h>


namespace OpenMS
{

  IsotopeLabelingMDVs::IsotopeLabelingMDVs()
  {
  }

  IsotopeLabelingMDVs::~IsotopeLabelingMDVs()
  {
  }

  
  void IsotopeLabelingMDVs::isotopicCorrection(
    Feature& normalized_feature,
    Feature& corrected_feature,
    std::vector<std::vector<double>> correction_matrix)
  {
    // MDV_corrected = correction_matrix_inversed * MDV_observed (normalized_features)
    
    uint8_t correction_matrix_n = correction_matrix.size();
    std::vector<std::vector<double>> correction_matrix_inversed(correction_matrix_n, std::vector<double>(correction_matrix_n,0));
    
    // 1- correction matrix inversion
    inverseMatrix_(correction_matrix, correction_matrix_inversed);
    
    // 2- element-wise expansion with MDV_observed
    std::vector<Feature> normalized_feature_subordinates = normalized_feature.getSubordinates();
    std::vector<double> MDV_observed;
    for (auto it = normalized_feature_subordinates.begin(); it != normalized_feature_subordinates.end(); it++)
    {
      MDV_observed.push_back(it->getMetaValue("peak_apex_int"));
    }
    
    double corrected_value;
    corrected_feature = normalized_feature;
    for (size_t i = 0; i < correction_matrix_inversed.size(); ++i) {
      corrected_value = 0.0;
      for (size_t j = 0; j < correction_matrix_inversed[0].size(); ++j) {
        corrected_value += correction_matrix_inversed[i][j] * MDV_observed[j];
      }
      corrected_feature.getSubordinates().at(i).setIntensity(corrected_value);
    }
    
  }


  void IsotopeLabelingMDVs::calculateIsotopicPurity(
    Feature& normalized_featuremap,
    Feature& featuremap_with_isotopic_purity,
    std::vector<double>& experiment_data,
    std::string& isotopic_purity_name)
  {
    
    normalized_featuremap = featuremap_with_isotopic_purity;
    
    double experiment_data_peak = 0.0;
    if ( !experiment_data.empty() )
    {
      std::vector<double>::iterator max_it  = std::max_element(experiment_data.begin(), experiment_data.end());
      uint64_t experiment_data_peak_idx     = std::distance(experiment_data.begin(), max_it);
      experiment_data_peak                  = experiment_data[experiment_data_peak_idx];
      
      if ( experiment_data_peak_idx >= 1 && experiment_data_peak != 0.0)
      {
        double previous_experiment_data_peak  = experiment_data[experiment_data_peak_idx - 1];
        double isotopic_purity                = experiment_data_peak_idx / (experiment_data_peak_idx + ( previous_experiment_data_peak / experiment_data_peak));
        featuremap_with_isotopic_purity.setMetaValue(isotopic_purity_name, isotopic_purity);
      }
    }

  }

  
  void IsotopeLabelingMDVs::calculateMDVAccuracy(FeatureMap& normalized_featuremap, FeatureMap& featuremap_with_accuracy_info)
  {
    // MARK: TODO calculateMDVAccuracy
  }


  void IsotopeLabelingMDVs::calculateMDV(
    Feature& measured_feature,
    Feature& normalized_feature,
    const String& mass_intensity_type,
    const String& feature_name)
  {
    
    std::vector<Feature> measured_feature_subordinates = measured_feature.getSubordinates();
    normalized_feature = measured_feature;
    if (mass_intensity_type == "norm_max")
    {
      if (feature_name == "intensity")
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getIntensity());
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
      
        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if(measured_feature_max != 0.0) {
            normalized_feature.getSubordinates().at(i).setIntensity(normalized_feature.getSubordinates().at(i).getIntensity() /  measured_feature_max);
          }
        }
      }
      
      else
      {
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getMetaValue(feature_name));
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
          
        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if (measured_feature_max != 0.0) {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / measured_feature_max);
          }
        }
      }
    }
    
    else if (mass_intensity_type == "norm_sum")
    {
      if (feature_name == "intensity")
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += it->getIntensity();
        }
        
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_feature.setMetaValue((it - measured_feature_subordinates.begin()), (it->getIntensity() / feature_peak_apex_intensity_sum));
        }
      }
      
      else
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += (Peak2D::IntensityType)it->getMetaValue(feature_name);
        }
          
        for (size_t i = 0; i < normalized_feature.getSubordinates().size(); ++i)
        {
          if (feature_peak_apex_intensity_sum != 0.0)
          {
            normalized_feature.getSubordinates().at(i).setIntensity((OpenMS::Peak2D::IntensityType)measured_feature_subordinates.at(i).getMetaValue(feature_name) / feature_peak_apex_intensity_sum);
          }
        }
      }
    }
    
  }

  template<typename T>
  void IsotopeLabelingMDVs::inverseMatrix_(
    std::vector<std::vector<T>>& correction_matrix,
    std::vector<std::vector<T>>& correction_matrix_inversed)
  {
    
    uint16_t correction_matrix_n = correction_matrix.size() == correction_matrix[0].size() ? correction_matrix.size() : 0;

    // 1- get the inverse
    double **CM_temp, **CM_inv;
    double temp;
    int i, j, k;

    CM_temp = (double **)malloc(correction_matrix_n * sizeof(double *));
    for(i = 0; i<correction_matrix_n; ++i)
    {
      CM_temp[i] = (double *)malloc(correction_matrix_n * sizeof(double));
    }

    CM_inv=(double **)malloc(correction_matrix_n * sizeof(double *));
    for(i = 0; i<correction_matrix_n; ++i)
    {
      CM_inv[i] = (double *)malloc(correction_matrix_n * sizeof(double));
    }

    for(i = 0; i < correction_matrix_n; ++i) {
      for(j = 0; j < correction_matrix_n; ++j) {
        CM_temp[i][j] = correction_matrix[i][j];
      }
    }

    // initialize as identity matrix
    for(i = 0; i < correction_matrix_n; ++i) {
      for(j = 0; j < correction_matrix_n; ++j) {
        if(i == j)
          CM_inv[i][j]=1;
        else
          CM_inv[i][j]=0;
      }
    }
    
    // inversion routine
    for(k = 0; k < correction_matrix_n; ++k)
    {
      temp = CM_temp[k][k];
     
      for(j = 0; j < correction_matrix_n; ++j)
      {
        CM_temp[k][j] /= temp;
        CM_inv[k][j]  /= temp;
      }
      for(i = 0; i < correction_matrix_n; ++i)
      {
        temp = CM_temp[i][k];
        for(j = 0; j < correction_matrix_n; ++j)
        {
          if(i == k)
            break;
          CM_temp[i][j] -= CM_temp[k][j] * temp;
          CM_inv[i][j]  -= CM_inv[k][j] * temp;
        }
      }
    }
    
    for(i = 0; i < correction_matrix_n; ++i)
    {
      for(j = 0; j < correction_matrix_n; ++j){
        correction_matrix_inversed[i][j] = CM_inv[i][j];
      }
    }
  
  }
  
  
} // namespace
