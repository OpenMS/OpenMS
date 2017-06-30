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

namespace OpenMS
{

  MRMTransitionGroupPicker::MRMTransitionGroupPicker() :
    DefaultParamHandler("MRMTransitionGroupPicker")
  {
    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("min_peak_width", -1.0, "Minimal peak width (s), discard all peaks below this value (-1 means no action).", ListUtils::create<String>("advanced"));

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated at the peak boundaries, either the smoothed or the raw chromatogram data can be used for that.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("background_subtraction", ListUtils::create<String>("none,smoothed,original"));

    defaults_.setValue("recalculate_peaks", "false", "Tries to get better peak picking by looking at peak consistency of all picked peaks. Tries to use the consensus (median) peak border if theof variation within the picked peaks is too large.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("recalculate_peaks", ListUtils::create<String>("true,false"));

    defaults_.setValue("use_precursors", "false", "Use precursor chromatogram for peak picking", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("use_precursors", ListUtils::create<String>("true,false"));

    defaults_.setValue("recalculate_peaks_max_z", 1.0, "Determines the maximal Z-Score (difference measured in standard deviations) that is considered too large for peak boundaries. If the Z-Score is above this value, the median is used for peak boundaries (default value 1.0).", ListUtils::create<String>("advanced"));

    defaults_.setValue("minimal_quality", -10000.0, "Only if compute_peak_quality is set, this parameter will not consider peaks below this quality threshold", ListUtils::create<String>("advanced"));

    defaults_.setValue("resample_boundary", 15.0, "For computing peak quality, how many extra seconds should be sample left and right of the actual peak", ListUtils::create<String>("advanced"));

    defaults_.setValue("compute_peak_quality", "false", "Tries to compute a quality value for each peakgroup and detect outlier transitions. The resulting score is centered around zero and values above 0 are generally good and below -1 or -2 are usually bad.", ListUtils::create<String>("advanced"));
    defaults_.setValidStrings("compute_peak_quality", ListUtils::create<String>("true,false"));

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
    min_qual_ = (double)param_.getValue("minimal_quality");
    min_peak_width_ = (double)param_.getValue("min_peak_width");
    resample_boundary_ = (double)param_.getValue("resample_boundary");
  }

  void MRMTransitionGroupPicker::calculateBgEstimation_(const MSChromatogram<>& chromatogram,
      double best_left, double best_right, double & background, double & avg_noise_level)
  {
    // determine (in the chromatogram) the intensity at the left / right border
    MSChromatogram<>::const_iterator it = chromatogram.begin();
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

  void MRMTransitionGroupPicker::findLargestPeak(std::vector<MSChromatogram<> >& picked_chroms, int& chr_idx, int& peak_idx)
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

}

