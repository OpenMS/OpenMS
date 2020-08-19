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

  
  void IsotopeLabelingMDVs::isotopicCorrection(FeatureMap& normalized_featuremap, FeatureMap& corrected_featuremap)
  {
    // MARK: TODO isotopicCorrection
  }


  void IsotopeLabelingMDVs::calculateIsotopicPurity(FeatureMap& normalized_featuremap, FeatureMap& featuremap_with_isotopic_purity)
  {
    // MARK: TODO calculateIsotopicPurity
  }

  
  void IsotopeLabelingMDVs::calculateMDVAccuracy(FeatureMap& normalized_featuremap, FeatureMap& featuremap_with_accuracy_info)
  {
    // MARK: TODO calculateMDVAccuracy
  }


  void IsotopeLabelingMDVs::calculateMDV(
    Feature& measured_feature,
    Feature& normalized_featuremap,
    const String& mass_intensity_type,
    const String& feature_name)
  {
    // Lactatex +
    //          |-> Lactate_x17 : peak_apex_int = 3.61E+08
    //          |-> Lactate_x18 : peak_apex_int = 1.20E+04
    //          |-> Lactate_x19 : peak_apex_int = 1.02E+05
    //          |-> Lactate_x20 : peak_apex_int = 2.59E+04
    
    std::vector<Feature> measured_feature_subordinates = measured_feature.getSubordinates();
    
    if (mass_intensity_type == "norm_max")
    {
      if (feature_name == "intensity")
      {
        measured_feature_subordinates = measured_feature.getSubordinates(); // i.e. Lactate1
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          intensities_vec.push_back(it->getIntensity());
        }
        std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
        double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
      
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_featuremap.setMetaValue((it - measured_feature_subordinates.begin()), (it->getIntensity() / measured_feature_max));
        }
      }
      
      else
      {
        if (measured_feature.metaValueExists(feature_name))
        {
          measured_feature_subordinates = measured_feature.getSubordinates(); // i.e. Lactate1
          std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
          for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
          {
            intensities_vec.push_back(it->getMetaValue(feature_name));
          }
          std::vector<OpenMS::Peak2D::IntensityType>::iterator max_it = std::max_element(intensities_vec.begin(), intensities_vec.end());
          double measured_feature_max = intensities_vec[std::distance(intensities_vec.begin(), max_it)];
          
          for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
          {
            normalized_featuremap.setMetaValue((it - measured_feature_subordinates.begin()), ((OpenMS::Peak2D::IntensityType)it->getMetaValue(feature_name) / measured_feature_max));
          }
        }
      }
    }
    
    else if (mass_intensity_type == "norm_sum")
    {
      
      if (feature_name == "intensity")
      {
        OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
        measured_feature_subordinates = measured_feature.getSubordinates(); // i.e. Lactate1
        std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          feature_peak_apex_intensity_sum += it->getIntensity();
        }
        
        for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
        {
          normalized_featuremap.setMetaValue((it - measured_feature_subordinates.begin()), (it->getIntensity() / feature_peak_apex_intensity_sum));
        }
      }
      
      else
      {
        if (measured_feature.metaValueExists(feature_name))
        {
          OpenMS::Peak2D::IntensityType feature_peak_apex_intensity_sum = 0.0;
          measured_feature_subordinates = measured_feature.getSubordinates(); // i.e. Lactate1
          std::vector<OpenMS::Peak2D::IntensityType> intensities_vec;
          for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
          {
            feature_peak_apex_intensity_sum += (Peak2D::IntensityType)it->getMetaValue(feature_name);
          }
          
          for (auto it = measured_feature_subordinates.begin(); it != measured_feature_subordinates.end(); it++)
          {
            normalized_featuremap.setMetaValue((it - measured_feature_subordinates.begin()), ((OpenMS::Peak2D::IntensityType)it->getMetaValue(feature_name) / feature_peak_apex_intensity_sum));
          }
        }
      }
    }
    
  }

} // namespace
