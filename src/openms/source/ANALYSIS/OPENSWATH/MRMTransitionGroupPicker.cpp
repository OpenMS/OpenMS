// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

namespace OpenMS
{

  // Simple linear interpolation at point x between x0 and x1
  double lin_interpolate(double x, double x0, double x1, double y0, double y1)
  {
    double slope = (y1 - y0) / (x1 - x0);
    double delta_y = (x - x0) * slope;
    return y0 + delta_y;
  }

  MRMTransitionGroupPicker::MRMTransitionGroupPicker() :
    DefaultParamHandler("MRMTransitionGroupPicker")
  {
    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("min_peak_width", 0.001, "Minimal peak width (s), discard all peaks below this value (-1 means no action).", {"advanced"});

    defaults_.setValue("peak_integration", "original", "Calculate the peak area and height either the smoothed or the raw chromatogram data.", {"advanced"});
    defaults_.setValidStrings("peak_integration", {"original","smoothed"});

    defaults_.setValue("background_subtraction", "none", "Remove background from peak signal using estimated noise levels. The 'original' method is only provided for historical purposes, please use the 'exact' method and set parameters using the PeakIntegrator: settings. The same original or smoothed chromatogram specified by peak_integration will be used for background estimation.", {"advanced"});
    defaults_.setValidStrings("background_subtraction", {"none","original","exact"});

    defaults_.setValue("recalculate_peaks", "false", "Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if the variation within the picked peaks is too large.", {"advanced"});
    defaults_.setValidStrings("recalculate_peaks", {"true","false"});

    defaults_.setValue("use_precursors", "false", "Use precursor chromatogram for peak picking (note that this may lead to precursor signal driving the peak picking)", {"advanced"});
    defaults_.setValidStrings("use_precursors", {"true","false"});

    defaults_.setValue("use_consensus", "true", "Use consensus peak boundaries when computing transition group picking (if false, compute independent peak boundaries for each transition)", {"advanced"});
    defaults_.setValidStrings("use_consensus", {"true","false"});

    defaults_.setValue("recalculate_peaks_max_z", 1.0, "Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).", {"advanced"});

    defaults_.setValue("minimal_quality", -10000.0, "Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold", {"advanced"});

    defaults_.setValue("resample_boundary", 15.0, "For computing peak quality, how many extra seconds should be sample left and right of the actual peak", {"advanced"});

    defaults_.setValue("compute_peak_quality", "false", "Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.", {"advanced"});
    defaults_.setValidStrings("compute_peak_quality", {"true","false"});
    
    defaults_.setValue("compute_peak_shape_metrics", "false", "Calculates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.", {"advanced"});
    defaults_.setValidStrings("compute_peak_shape_metrics", {"true","false"});

    defaults_.setValue("compute_total_mi", "false", "Compute mutual information metrics for individual transitions that can be used for OpenSWATH/IPF scoring.", {"advanced"});
    defaults_.setValidStrings("compute_total_mi", {"true","false"});

    defaults_.setValue("boundary_selection_method", "largest", "Method to use when selecting the best boundaries for peaks.", {"advanced"});
    defaults_.setValidStrings("boundary_selection_method", {"largest","widest"});

    defaults_.insert("PeakPickerChromatogram:", PeakPickerChromatogram().getDefaults());
    defaults_.insert("PeakIntegrator:", PeakIntegrator().getDefaults());

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMTransitionGroupPicker::~MRMTransitionGroupPicker() = default;

  MRMTransitionGroupPicker& MRMTransitionGroupPicker::operator=(const MRMTransitionGroupPicker& rhs)
  {
    if (&rhs == this)
      return *this;

    // don't copy parameters

    return *this;
  }

  void MRMTransitionGroupPicker::updateMembers_()
  {
    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
    stop_after_intensity_ratio_ = (double)param_.getValue("stop_after_intensity_ratio");
    peak_integration_ = param_.getValue("peak_integration").toString();
    background_subtraction_ = param_.getValue("background_subtraction").toString();
    recalculate_peaks_ = (bool)param_.getValue("recalculate_peaks").toBool();
    use_precursors_ = (bool)param_.getValue("use_precursors").toBool();
    use_consensus_ = (bool)param_.getValue("use_consensus").toBool();
    recalculate_peaks_max_z_ = (double)param_.getValue("recalculate_peaks_max_z");
    compute_peak_quality_ = (bool)param_.getValue("compute_peak_quality").toBool();
    compute_peak_shape_metrics_ = (bool)param_.getValue("compute_peak_shape_metrics").toBool();
    compute_total_mi_ = (bool)param_.getValue("compute_total_mi").toBool();
    min_qual_ = (double)param_.getValue("minimal_quality");
    min_peak_width_ = (double)param_.getValue("min_peak_width");
    resample_boundary_ = (double)param_.getValue("resample_boundary");
    boundary_selection_method_ = param_.getValue("boundary_selection_method").toString();

    picker_.setParameters(param_.copy("PeakPickerChromatogram:", true));
    pi_.setParameters(param_.copy("PeakIntegrator:", true));
  }

  void MRMTransitionGroupPicker::findLargestPeak(const std::vector<MSChromatogram >& picked_chroms, int& chr_idx, int& peak_idx)
  {
    double largest = 0.0;
    for (Size k = 0; k < picked_chroms.size(); k++)
    {
      for (Size i = 0; i < picked_chroms[k].size(); i++)
      {
        if (picked_chroms[k][i].getIntensity() > largest)
        {
          largest = picked_chroms[k][i].getIntensity();
          chr_idx = (int)k;
          peak_idx = (int)i;
        }
      }
    }
  }

  void MRMTransitionGroupPicker::findWidestPeakIndices(const std::vector<MSChromatogram>& picked_chroms, Int& chrom_idx, Int& point_idx) const
  {
    double max_width{0};
    for (Size i = 0; i < picked_chroms.size(); ++i)
    {
      for (Size k = 0; k < picked_chroms[i].size(); ++k)
      {
        const double left_rt = picked_chroms[i].getFloatDataArrays()[PeakPickerChromatogram::IDX_LEFTBORDER][k];
        const double right_rt = picked_chroms[i].getFloatDataArrays()[PeakPickerChromatogram::IDX_RIGHTBORDER][k];
        const double local_peak_width = right_rt - left_rt;
        OPENMS_LOG_DEBUG << "findWidestPeakIndices(): local_peak_width=" << local_peak_width << std::endl;
        if (local_peak_width > max_width)
        {
          max_width = local_peak_width;
          chrom_idx = static_cast<Int>(i);
          point_idx = static_cast<Int>(k);
          OPENMS_LOG_DEBUG << "findWidestPeakIndices(): max_width=" << max_width << "; chrom_idx=" << chrom_idx << "; point_idx=" << point_idx << std::endl;
        }
      }
    }
  }
}
