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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMTransitionGroupPicker.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

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
    defaults_.setValue("min_peak_width", -1.0, "Minimal peak width (s), discard all peaks below this value (-1 means no action).", ListUtils::create<String>("advanced"));

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated at the peak boundaries, either the smoothed or the raw chromatogram data can be used for that.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("background_subtraction", ListUtils::create<String>("none,smoothed_average,smoothed_exact,original_average,original_exact"));

    defaults_.setValue("recalculate_peaks", "false", "Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if theof variation within the picked peaks is too large.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("recalculate_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("use_precursors", "false", "Use precursor chromatogram for peak picking", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("use_precursors", ListUtils::create<String>("true,false"));

    defaults_.setValue("recalculate_peaks_max_z", 1.0, "Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).", ListUtils::create<String>("advanced"));

    defaults_.setValue("minimal_quality", -10000.0, "Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold", ListUtils::create<String>("advanced"));

    defaults_.setValue("resample_boundary", 15.0, "For computing peak quality, how many extra seconds should be sample left and right of the actual peak", ListUtils::create<String>("advanced"));

    defaults_.setValue("compute_peak_quality", "false", "Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("compute_peak_quality", ListUtils::create<String>("true,false"));
    
    defaults_.setValue("compute_peak_shape_metrics", "false", "Calulates various peak shape metrics (e.g., tailing) that can be used for downstream QC/QA.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("compute_peak_shape_metrics", ListUtils::create<String>("true,false"));

    defaults_.insert("PeakPickerMRM:", PeakPickerMRM().getDefaults());

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  MRMTransitionGroupPicker::~MRMTransitionGroupPicker()
  {
  }

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
    background_subtraction_ = param_.getValue("background_subtraction");
    recalculate_peaks_ = (bool)param_.getValue("recalculate_peaks").toBool();
    use_precursors_ = (bool)param_.getValue("use_precursors").toBool();
    recalculate_peaks_max_z_ = (double)param_.getValue("recalculate_peaks_max_z");
    compute_peak_quality_ = (bool)param_.getValue("compute_peak_quality").toBool();
    compute_peak_shape_metrics_ = (bool)param_.getValue("compute_peak_shape_metrics").toBool();
    min_qual_ = (double)param_.getValue("minimal_quality");
    min_peak_width_ = (double)param_.getValue("min_peak_width");
    resample_boundary_ = (double)param_.getValue("resample_boundary");

    picker_.setParameters(param_.copy("PeakPickerMRM:", true));
  }
  
  void MRMTransitionGroupPicker::calculateBgEstimationAverage_(const MSChromatogram& chromatogram,
      double best_left, double best_right, double & background, double & avg_noise_level)
  {
    // determine (in the chromatogram) the intensity at the left / right border
    MSChromatogram::const_iterator it = chromatogram.begin();
    int nr_points = 0;
    for (; it != chromatogram.end(); ++it)
    {
      if (it->getMZ() > best_left)
      {
        nr_points++;
        break;
      }
    }
    double intensity_left = it->getIntensity();
    for (; it != chromatogram.end(); ++it)
    {
      if (it->getMZ() > best_right)
      {
        break;
      }
      nr_points++;
    }
    if (it == chromatogram.begin() || nr_points < 1)
    {
      // something is fishy, the endpoint of the peak is the beginning of the chromatogram
      std::cerr << "Tried to calculate background but no points were found " << std::endl;
      return;
    }

    // decrease the iterator and the nr_points by one (because we went one too far)
    double intensity_right = (it--)->getIntensity();
    nr_points--;

    avg_noise_level = (intensity_right + intensity_left) / 2;
    background = avg_noise_level * nr_points;
  }

  void MRMTransitionGroupPicker::calculateBgEstimationExact_(const MSChromatogram& chromatogram,
      double best_left, double best_right, double peak_height, double & background, double & avg_noise_level)
  {
    // determine (in the chromatogram) the intensity at the left / right border
    double intensity_left = 0.0;
    double rt_apex = 0.0;
    double intensity_right = 0.0;
    
    // calculate the average noise level
    for (MSChromatogram::const_iterator it = std::next(chromatogram.begin()); it != chromatogram.end(); ++it)
    {
      MSChromatogram::const_iterator it_prev = it;
      --it_prev; //previous point

      if (it->getMZ() >= best_left && it_prev->getMZ() < best_left)
      {
        intensity_left = it->getIntensity();
      }
      else if (it->getIntensity() >= peak_height && it_prev->getIntensity() < peak_height)
      {
        rt_apex = it->getMZ();
      }
      else if (it->getMZ() >= best_right && it_prev->getMZ() < best_right)
      {
        intensity_right = it->getIntensity();
      }
    }

    double intensity_max, intensity_min, rt_min;
    if (intensity_left >= intensity_right)
    {
      intensity_max = intensity_left;
      intensity_min = intensity_right;
      rt_min = best_right;
    }
    else 
    {
      intensity_max = intensity_right;
      intensity_min = intensity_left;
      rt_min = best_left;
    }
    // calculate the average noise level using the sin/cos rule
    double delta_int = intensity_max - intensity_min;
    double delta_rt = best_right - best_left;
    double delta_rt_apex = std::fabs(rt_min-rt_apex);
    double delta_int_apex = delta_int*delta_rt_apex/delta_rt;

    avg_noise_level = intensity_min + delta_int_apex;
    //NOTE: formula for calculating the background using the trapezoidal rule (future PR)
    // background = intensity_min*delta_rt + 0.5*delta_int*delta_rt;

    // calculate the background
    background = 0.0;
    for (MSChromatogram::const_iterator it = std::next(chromatogram.begin()); it != chromatogram.end(); ++it)
    {
      MSChromatogram::const_iterator it_prev = it;
      --it_prev; //previous point

      if (it->getMZ() >= best_left && it_prev->getMZ() < best_right)
      {
        // calculate the background using the formula
        // y = mx + b where x = retention time, m = slope, b = left intensity
        double delta_int = intensity_right - intensity_left; // sign will determine line direction
        double delta_rt_current = (it->getMZ() - best_left);
        double background_int_current = delta_int/delta_rt*delta_rt_current + intensity_left;
        background = background + background_int_current;
      }
    }
  }

  void MRMTransitionGroupPicker::findLargestPeak(std::vector<MSChromatogram >& picked_chroms, int& chr_idx, int& peak_idx)
  {
    double largest = 0.0;
    ChromatogramPeak largest_pos;
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
  
  void MRMTransitionGroupPicker::calculatePeakApexInt_(const MSChromatogram& chromatogram,
    double best_left, double best_right, 
    ConvexHull2D::PointArrayType & hull_points,
    double & intensity_sum, 
    double & intensity_integral,
    double & rt_sum,
    double & peak_apex_int,
    double & peak_apex_rt)
  {
    intensity_sum = 0.0;
    rt_sum = 0.0;
    peak_apex_int = -1;
    double peak_apex_dist = std::fabs(chromatogram.begin()->getMZ() - peak_apex_rt);
    int peak_nr (0);
    bool added_right(false);
    // FEATURE : use RTBegin / MZBegin -> for this we need to know whether the template param is a real chromatogram or a spectrum!
    MSChromatogram::const_iterator prev_it = chromatogram.begin();
    for (MSChromatogram::const_iterator it = chromatogram.begin(); it != chromatogram.end(); ++it)
    {
      if (it->getMZ() > best_left && it->getMZ() < best_right)
      {
        if (peak_nr == 0 && it != chromatogram.begin())
        {
          // add area between first measured peak inside the peak boundaries and start of peak to the left
          double delta_rt = it->getRT() - best_left;
          double interpol_intensity = lin_interpolate(best_left, prev_it->getRT(), it->getRT(), prev_it->getIntensity(), it->getIntensity());
          intensity_integral += (interpol_intensity + it->getIntensity())/2.0 * delta_rt;
        }

        if (peak_nr > 0)
        {
          // add area between last peak (inside boundaries) and current peak (inside boundaries)
          double delta_rt = it->getRT() - prev_it->getRT();
          intensity_integral += (prev_it->getIntensity() + it->getIntensity())/2.0 * delta_rt;
        }

        DPosition<2> p;
        p[0] = it->getMZ();
        p[1] = it->getIntensity();
        hull_points.push_back(p);
        if (std::fabs(it->getMZ() - peak_apex_rt) <= peak_apex_dist)
        {
          peak_apex_int = p[1];
          peak_apex_dist = std::fabs(it->getMZ() - peak_apex_rt);
        }
        rt_sum += it->getMZ();
        intensity_sum += it->getIntensity();
        peak_nr++;
      }
      else if (peak_nr > 0 && !added_right)
      {
        // add area between last measured peak inside the peak boundaries and end of peak to the right
        double delta_rt = best_right - prev_it->getRT();
        double interpol_intensity = lin_interpolate(best_right, prev_it->getRT(), it->getRT(), prev_it->getIntensity(), it->getIntensity());
        intensity_integral += (prev_it->getIntensity() + interpol_intensity)/2.0 * delta_rt;
        added_right = true;
      }
      prev_it = it;
    }
  }

  void MRMTransitionGroupPicker::calculatePeakShapeMetrics_(const MSChromatogram& chromatogram, 
    double best_left, double best_right, 
    double peak_height, double peak_apex_rt, double avg_noise_level,
    PeakShapeMetrics_ & peakShapeMetrics)
  {
    peakShapeMetrics.points_across_baseline = 0;
    double start_intensity(0), end_intensity(0);
    double delta_rt, delta_int, height_5, height_10, height_50;
    
    for (MSChromatogram::const_iterator it = chromatogram.begin() + 1; it != chromatogram.end(); ++it)
    {
      MSChromatogram::const_iterator it_prev = it;
      --it_prev; //previous point
      double intensity = std::max(it->getIntensity()-avg_noise_level, 0.0); //background-subtracted intensity
      double intensity_prev = std::max(it_prev->getIntensity()-avg_noise_level, 0.0); //background-subtracted intensity of the previous point
      double retention_time = it->getMZ();
      double retention_time_prev = it_prev->getMZ();

      // start and end intensities
      if (retention_time_prev < best_left && retention_time >= best_left)
      {
        start_intensity = intensity_prev;
      }
      else if (retention_time_prev < best_right && retention_time >= best_right)
      {
        end_intensity = intensity;
      }

      if (retention_time >= best_left && retention_time <= best_right)
      {
        //start and end retention times
        if (retention_time < peak_apex_rt)
        {
          // start_time_at_5
          if (intensity >= 0.05*peak_height && 
            intensity_prev < 0.05*peak_height && 
            peakShapeMetrics.points_across_baseline > 1)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity - intensity_prev;
            height_5 = intensity - 0.05*peak_height;
            peakShapeMetrics.start_time_at_5 = retention_time - delta_int*delta_rt/height_5;
          }
          // start_time_at_10
          if (intensity >= 0.1*peak_height && 
            intensity_prev < 0.1*peak_height && 
            peakShapeMetrics.points_across_baseline > 1)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity - intensity_prev;
            height_10 = intensity - 0.1*peak_height;
            peakShapeMetrics.start_time_at_10 = retention_time - delta_int*delta_rt/height_10;
          }
          // start_time_at_50
          if (intensity >= 0.5*peak_height && 
            intensity_prev < 0.5*peak_height && 
            peakShapeMetrics.points_across_baseline > 1)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity - intensity_prev;
            height_50 = intensity - 0.5*peak_height;
            peakShapeMetrics.start_time_at_50 = retention_time - delta_int*delta_rt/height_50;
          }
        } 
        else if (retention_time > peak_apex_rt)
        {
          // end_time_at_5
          if (intensity <= 0.05*peak_height && 
            intensity_prev > 0.05*peak_height)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity_prev - intensity;
            height_5 = 0.05*peak_height - intensity;
            peakShapeMetrics.end_time_at_5 = retention_time - delta_int*delta_rt/height_5;
          }
          // start_time_at_10
          if (intensity <= 0.1*peak_height && 
            intensity_prev > 0.1*peak_height)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity_prev - intensity;
            height_10 = 0.1*peak_height - intensity;
            peakShapeMetrics.end_time_at_10 = retention_time - delta_int*delta_rt/height_10;
          }
          // end_time_at_50
          if (intensity <= 0.5*peak_height && 
          intensity_prev > 0.5*peak_height)
          {
            delta_rt = retention_time - retention_time_prev;
            delta_int = intensity_prev - intensity;
            height_50 = 0.5*peak_height - intensity;
            peakShapeMetrics.end_time_at_50 = retention_time - delta_int*delta_rt/height_50;
          }
        } 

        // points across the peak
        peakShapeMetrics.points_across_baseline ++;
        if (intensity >= 0.5*peak_height)
        {
          peakShapeMetrics.points_across_half_height ++;
        }
      }
    }

    // peak widths
    peakShapeMetrics.width_at_5 = peakShapeMetrics.end_time_at_5 - peakShapeMetrics.start_time_at_5;
    peakShapeMetrics.width_at_10 = peakShapeMetrics.end_time_at_10 - peakShapeMetrics.start_time_at_10;
    peakShapeMetrics.width_at_50 = peakShapeMetrics.end_time_at_50 - peakShapeMetrics.start_time_at_50;
    peakShapeMetrics.total_width = best_right - best_left;
    peakShapeMetrics.slope_of_baseline = end_intensity - start_intensity;
    peakShapeMetrics.baseline_delta_2_height = peakShapeMetrics.slope_of_baseline / peak_height;

    // other
    peakShapeMetrics.tailing_factor = peakShapeMetrics.width_at_5 / std::min(peak_apex_rt - peakShapeMetrics.start_time_at_5, peakShapeMetrics.end_time_at_5 - peak_apex_rt);
    peakShapeMetrics.asymmetry_factor = std::min(peak_apex_rt - peakShapeMetrics.start_time_at_10, peakShapeMetrics.end_time_at_10 - peak_apex_rt) / std::max(peak_apex_rt - peakShapeMetrics.start_time_at_10, peakShapeMetrics.end_time_at_10 - peak_apex_rt);
  }

}

