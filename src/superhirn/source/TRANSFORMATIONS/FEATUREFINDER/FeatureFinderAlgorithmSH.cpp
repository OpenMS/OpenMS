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
// $Maintainer: Timo Sachsenberg $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSH.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSHCtrl.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <boost/shared_ptr.hpp>

namespace OpenMS
{
  FeatureFinderAlgorithmSH::FeatureFinderAlgorithmSH() :
    FeatureFinderAlgorithm()
  {
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:active", "false", "MS1 data centroid data");
    this->defaults_.setValidStrings("centroiding:active", ListUtils::create<String>("true,false"));
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:precursor_detection_scan_levels", ListUtils::create<Int>("1"), "Precursor detection scan levels");
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:max_inter_scan_distance", 0, "MS1 max inter scan distance"); // was 0.1
    this->defaults_.setMinInt("ms1:max_inter_scan_distance", 0); // Markus needs to clarify this parameter
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:tr_resolution", 0.01, "MS1 LC retention time resolution"); // seems to have no effect
    this->defaults_.setMinFloat("ms1:tr_resolution", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:intensity_threshold", 1000.0, "FT peak detect MS1 intensity min threshold");
    this->defaults_.setMinFloat("ms1:intensity_threshold", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:max_inter_scan_rt_distance", 0.1, "MS1 max inter scan distance"); // seems to have no effect
    this->defaults_.setMinFloat("ms1:max_inter_scan_rt_distance", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:min_nb_cluster_members", 4, "FT peak detect MS1 min nb peak members");
    this->defaults_.setMinInt("ms1:min_nb_cluster_members", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:detectable_isotope_factor", 0.05, "Detectable isotope factor");
    this->defaults_.setMinFloat("ms1:detectable_isotope_factor", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:intensity_cv", 0.9, "IntensityCV");
    this->defaults_.setMinFloat("ms1:intensity_cv", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:window_width", 5, "Centroid window width");
    this->defaults_.setMinInt("centroiding:window_width", 1);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:absolute_isotope_mass_precision", 0.01, "Absolute isotope mass precision (Da)");
    this->defaults_.setMinFloat("centroiding:absolute_isotope_mass_precision", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:relative_isotope_mass_precision", 10.0, "Relative isotope mass precision");
    this->defaults_.setMinFloat("centroiding:relative_isotope_mass_precision", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:minimal_peak_height", 0.0, "Minimal peak height");
    this->defaults_.setMinFloat("centroiding:minimal_peak_height", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("centroiding:min_ms_signal_intensity", 50.0, "Minimal Centroid MS Signal Intensity");
    this->defaults_.setMinFloat("centroiding:min_ms_signal_intensity", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:retention_time_tolerance", 0.5, "MS1 retention time tolerance (minutes)");
    this->defaults_.setMinFloat("ms1:retention_time_tolerance", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1:mz_tolerance", 0.0, "MS1 m/z tolerance (ppm)");
    this->defaults_.setMinFloat("ms1:mz_tolerance", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:active", "true", "Activation of MS1 feature merging post processing");
    this->defaults_.setValidStrings("ms1_feature_merger:active", ListUtils::create<String>("true,false"));
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:tr_resolution", 0.01, "MS1 LC retention time resolution");
    this->defaults_.setMinFloat("ms1_feature_merger:tr_resolution", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:initial_apex_tr_tolerance", 5.0, "Initial Apex Tr tolerance");
    this->defaults_.setMinFloat("ms1_feature_merger:initial_apex_tr_tolerance", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:feature_merging_tr_tolerance", 1.0, "MS1 feature Tr merging tolerance");
    this->defaults_.setMinFloat("ms1_feature_merger:feature_merging_tr_tolerance", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:intensity_variation_percentage", 25.0, "Percentage of intensity variation between LC border peaks");
    this->defaults_.setMinFloat("ms1_feature_merger:intensity_variation_percentage", 0.0);
    this->defaults_.setMaxFloat("ms1_feature_merger:intensity_variation_percentage", 100.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_merger:ppm_tolerance_for_mz_clustering", 10.0, "PPM value for the m/z clustering of merging candidates");
    this->defaults_.setMinFloat("ms1_feature_merger:ppm_tolerance_for_mz_clustering", 0.0);
    // ----------------------------------------------------------------------------------------------------
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:start_elution_window", 0.0, "start elution window (minutes)");
    this->defaults_.setMinFloat("ms1_feature_selection_options:start_elution_window", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:end_elution_window", 180.0, "end elution window (minutes)");
    this->defaults_.setMinFloat("ms1_feature_selection_options:end_elution_window", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:mz_range_min", 0.0, "MS1 feature mz range min");
    this->defaults_.setMinFloat("ms1_feature_selection_options:mz_range_min", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:mz_range_max", 2000.0, "MS1 feature mz range max");
    this->defaults_.setMinFloat("ms1_feature_selection_options:mz_range_max", 0.0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:chrg_range_min", 1, "MS1 feature CHRG range min");
    this->defaults_.setMinInt("ms1_feature_selection_options:chrg_range_min", 0);
    // ----------------------------------------------------------------------------------------------------
    this->defaults_.setValue("ms1_feature_selection_options:chrg_range_max", 5, "MS1 feature CHRG range max");
    this->defaults_.setMinInt("ms1_feature_selection_options:chrg_range_max", 0);

    this->check_defaults_ =  false;
  }

  unsigned int FeatureFinderAlgorithmSH::getNativeScanId(String native_id)
  {

    Size start_idx = 0;
    while (start_idx < native_id.length() && !isdigit(native_id[start_idx]))
    {
      ++start_idx;
    }
    if (start_idx == native_id.length())
    {
      std::cout << "Native id could not be determined: " << native_id;
      throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Cannot convert native id to unsigned integer");
    }

    Size end_idx = start_idx;
    while (isdigit(native_id[end_idx]))
    {
      ++end_idx;
    }

    return native_id.substr(start_idx, end_idx - start_idx).toInt();
  }

  void FeatureFinderAlgorithmSH::run()
  {
    std::cout << "SuperHirn feature extraction...\n";

    map_ = *(FeatureFinderAlgorithm::map_);

    FeatureFinderAlgorithmSHCtrl::MyMap dummyMap;
    FeatureFinderAlgorithmSHCtrl::Vec datavec;
    datavec.resize(map_.size());
    unsigned int scanId = 0;

    // Ordering by native IDs order by scan numbers
    // To achieve the exact same results as the original
    // superhirn does, this is necessary.
    // However, its is very experimental and will work
    // for all data since its based on string comparison.
    bool orderByNativeIds = false;

    // go through map, extract data and store it in a vector of RawData objects
    for (unsigned int s = 0; s < map_.size(); s++)
    {
      const SpectrumType& spectrum = map_[s];
      double rt = spectrum.getRT();

      if (orderByNativeIds)
      {
        scanId = getNativeScanId(spectrum.getNativeID());
        if (scanId == 0)
        {
          std::cout << "Order by native ids not working, turning it off.\n";
          orderByNativeIds = false;
          scanId = 1;
        }
      }
      else
      {
        scanId++;
      }

      std::vector<double> vmzvals;
      std::vector<double> vintvals;

      for (Size p = 0; p < spectrum.size(); ++p)
      {
        vmzvals.push_back(spectrum[p].getMZ());
        vintvals.push_back(spectrum[p].getIntensity());
      }

      //RawData* data = new RawData(vmzvals, vintvals);
      boost::shared_ptr<RawData> data_ptr(new RawData(vmzvals, vintvals));

      FeatureFinderAlgorithmSHCtrl::MyMap map_ptr(rt / 60, data_ptr);
//        m[rt/60.0] = data;
      unsigned int scanIndex = scanId - 1;
      datavec[scanIndex] = map_ptr;
    }

    // apply the SuperHirn FeatureFinder algorithm
    FeatureFinderAlgorithmSHCtrl ctrl;
    ctrl.initParams(this->param_);
    std::vector<Feature> thefeatures = ctrl.extractPeaks(datavec);

    for (unsigned int i = 0; i < thefeatures.size(); ++i)
    {
      features_->push_back(thefeatures[i]);
    }
  }

  FeatureFinderAlgorithm* FeatureFinderAlgorithmSH::create()
  {
    return new FeatureFinderAlgorithmSH();
  }

  const String FeatureFinderAlgorithmSH::getProductName()
  {
    return "superhirn";
  }

}
