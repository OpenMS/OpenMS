// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Florian Zeller $
// $Authors: Florian Zeller $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERALGORITHMSH_H

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithmSHCtrl.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>

namespace OpenMS
{
	/** 
		@brief The Superhirn FeatureFinderAlgorithm.
	
		@ingroup FeatureFinder
	*/
  template<class PeakType, class FeatureType> class FeatureFinderAlgorithmSH :
    public FeatureFinderAlgorithm<PeakType, FeatureType>,
    public FeatureFinderDefs
  {
    
  public:
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::MapType MapType; // MSExperiment
    typedef typename FeatureFinderAlgorithm<PeakType, FeatureType>::FeatureMapType FeatureMapType;
    typedef typename MapType::SpectrumType SpectrumType;
    
    using FeatureFinderAlgorithm<PeakType, FeatureType>::features_;
    
    FeatureFinderAlgorithmSH() : FeatureFinderAlgorithm<PeakType, FeatureType>()
    {
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:active", "false", "MS1 data centroid data");
      this->defaults_.setValidStrings("centroiding:active", StringList::create("true,false"));
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:precursor_detection_scan_levels", IntList::create(1), "Precursor detection scan levels");
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:max_inter_scan_distance", 0, "MS1 max inter scan distance"); // was 0.1
      this->defaults_.setMinInt(      "ms1:max_inter_scan_distance", 0); // Markus needs to clarify this parameter
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:tr_resolution", 0.01, "MS1 LC retention time resolution");  // seems to have no effect
      this->defaults_.setMinFloat(    "ms1:tr_resolution", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:intensity_threshold", 1000.0, "FT peak detect MS1 intensity min threshold");
      this->defaults_.setMinFloat(    "ms1:intensity_threshold", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:max_inter_scan_rt_distance", 0.1, "MS1 max inter scan distance"); // seems to have no effect
      this->defaults_.setMinFloat(    "ms1:max_inter_scan_rt_distance", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:min_nb_cluster_members", 4, "FT peak detect MS1 min nb peak members");
      this->defaults_.setMinInt(      "ms1:min_nb_cluster_members", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:detectable_isotope_factor", 0.05, "Detectable isotope factor");
      this->defaults_.setMinFloat(    "ms1:detectable_isotope_factor", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:intensity_cv", 0.9, "IntensityCV");
      this->defaults_.setMinFloat(    "ms1:intensity_cv", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:window_width", 5, "Centroid window width");
      this->defaults_.setMinInt(      "centroiding:window_width", 1);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:absolute_isotope_mass_precision", 0.01, "Absolute isotope mass precision (Da)");
      this->defaults_.setMinFloat(    "centroiding:absolute_isotope_mass_precision", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:relative_isotope_mass_precision", 10.0, "Relative isotope mass precision");
      this->defaults_.setMinFloat(    "centroiding:relative_isotope_mass_precision", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:minimal_peak_height", 0.0, "Minimal peak height");
      this->defaults_.setMinFloat(    "centroiding:minimal_peak_height", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "centroiding:min_ms_signal_intensity", 50.0, "Minimal Centroid MS Signal Intensity");
      this->defaults_.setMinFloat(    "centroiding:min_ms_signal_intensity", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:retention_time_tolerance", 0.5, "MS1 retention time tolerance (minutes)");
      this->defaults_.setMinFloat(    "ms1:retention_time_tolerance", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1:mz_tolerance", 0.0, "MS1 m/z tolerance (ppm)");
      this->defaults_.setMinFloat(    "ms1:mz_tolerance", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:active", "true", "Activation of MS1 feature merging post processing");
      this->defaults_.setValidStrings("ms1_feature_merger:active", StringList::create("true,false"));
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:tr_resolution", 0.01, "MS1 LC retention time resolution");
      this->defaults_.setMinFloat(    "ms1_feature_merger:tr_resolution", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:initial_apex_tr_tolerance", 5.0, "Initial Apex Tr tolerance");
      this->defaults_.setMinFloat(    "ms1_feature_merger:initial_apex_tr_tolerance", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:feature_merging_tr_tolerance", 1.0, "MS1 feature Tr merging tolerance");
      this->defaults_.setMinFloat(    "ms1_feature_merger:feature_merging_tr_tolerance", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:intensity_variation_percentage", 25.0, "Percentage of intensity variation between LC border peaks");
      this->defaults_.setMinFloat(    "ms1_feature_merger:intensity_variation_percentage", 0.0);
      this->defaults_.setMaxFloat(    "ms1_feature_merger:intensity_variation_percentage", 100.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_merger:ppm_tolerance_for_mz_clustering", 10.0, "PPM value for the m/z clustering of merging candidates");
      this->defaults_.setMinFloat(    "ms1_feature_merger:ppm_tolerance_for_mz_clustering", 0.0);
      // ----------------------------------------------------------------------------------------------------
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:start_elution_window", 0.0, "start elution window (minutes)");
      this->defaults_.setMinFloat(    "ms1_feature_selection_options:start_elution_window", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:end_elution_window", 180.0, "end elution window (minutes)");
      this->defaults_.setMinFloat(    "ms1_feature_selection_options:end_elution_window", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:mz_range_min", 0.0, "MS1 feature mz range min");
      this->defaults_.setMinFloat(    "ms1_feature_selection_options:mz_range_min", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:mz_range_max", 2000.0, "MS1 feature mz range max");
      this->defaults_.setMinFloat(    "ms1_feature_selection_options:mz_range_max", 0.0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:chrg_range_min", 1, "MS1 feature CHRG range min");
      this->defaults_.setMinInt(      "ms1_feature_selection_options:chrg_range_min", 0);
      // ----------------------------------------------------------------------------------------------------
      this->defaults_.setValue(       "ms1_feature_selection_options:chrg_range_max", 5, "MS1 feature CHRG range max");
      this->defaults_.setMinInt(      "ms1_feature_selection_options:chrg_range_max", 0);
      
      this->check_defaults_ =  false;
    }
    
    unsigned int getNativeScanId(String native_id)
    {
      
      Size start_idx=0;
      while(!isdigit(native_id[start_idx]) && start_idx<native_id.length())
      {
        ++start_idx;
      }
      if(start_idx==native_id.length())
      {
        std::cout << "Native id could not be determined: " << native_id;
        throw Exception::InvalidParameter(__FILE__,__LINE__,__PRETTY_FUNCTION__, "Cannot convert native id to unsigned integer");
      }
      
      Size end_idx = start_idx;
      while (isdigit(native_id[end_idx]))
      {
        ++end_idx;
      }
      
      return native_id.substr(start_idx, end_idx-start_idx).toInt();
    }
    
    virtual void run()
    {
      std::cout << "SuperHirn feature extraction...\n";
      
      map_ = *(FeatureFinderAlgorithm<PeakType, FeatureType>::map_);
      
      MyMap dummyMap;
      Vec* datavec = new Vec(map_.size(), dummyMap);
      unsigned int scanId = 0;
      
      // Ordering by native IDs order by scan numbers
      // To achieve the exact same results as the original
      // superhirn does, this is necessary.
      // However, its is very experimental and will work
      // for all data since its based on string comparison.
      bool orderByNativeIds = false;
      
      for (unsigned int s = 0; s < map_.size(); s++)
      {
        const SpectrumType& spectrum = map_[s];
        double rt = spectrum.getRT();
        
        if (orderByNativeIds) 
        {
          scanId = getNativeScanId(spectrum.getNativeID());
          if (scanId == 0) {
            std::cout << "Order by native ids not working, turning it off.\n";
            orderByNativeIds = false;
            scanId = 1;
          }
        }
        else {
          scanId++;
        }
        
        std::vector<double>* vmzvals = new std::vector<double>();
        std::vector<double>* vintvals = new std::vector<double>();
        
        for (Size p=0; p<spectrum.size(); ++p)
        {
          vmzvals->push_back(spectrum[p].getMZ());
          vintvals->push_back(spectrum[p].getIntensity());
        }
        
        RawData* data = new RawData(*vmzvals, *vintvals);
        
        MyMap m;
        m[rt/60.0] = data;
        unsigned int scanIndex = scanId - 1;
        datavec->at(scanIndex) = m;
      }
      
      FeatureFinderAlgorithmSHCtrl ctrl;
      ctrl.initParams(this->param_);
      std::vector<Feature> thefeatures = ctrl.extractPeaks(*datavec);
      
      for (unsigned int i=0; i<thefeatures.size(); ++i)
        features_->push_back(thefeatures[i]);
    }
    
    static FeatureFinderAlgorithm<Peak1D,Feature>* create()
    {
      return new FeatureFinderAlgorithmSH();
    }
    
    static const String getProductName()
    {
      return "superhirn";
    }
    
    
  protected:
    MapType map_;
    
  };
  
}

#endif
