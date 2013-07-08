// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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

  MRMTransitionGroupPicker::MRMTransitionGroupPicker() :
    DefaultParamHandler("MRMTransitionGroupPicker")
  {
    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("min_peak_width", -1.0, "Minimal peak width (s)");

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated at the peak boundaries, either the smoothed or the raw chromatogram data can be used for that."); //, StringList::create("advanced"));
    defaults_.setValidStrings("background_subtraction", StringList::create("none,smoothed,original"));

    defaults_.setValue("recalculate_peaks", "false", "Tries to get better peak picking by looking at peak consistency");

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

    // dont copy parameters

    return *this;
  }

  void MRMTransitionGroupPicker::updateMembers_()
  {
    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
    stop_after_intensity_ratio_ = (DoubleReal)param_.getValue("stop_after_intensity_ratio");
    background_subtraction_ = param_.getValue("background_subtraction");
    recalculate_peaks_ = (bool)param_.getValue("recalculate_peaks").toBool();
    min_peak_width_ = (DoubleReal)param_.getValue("min_peak_width");
  }

  double MRMTransitionGroupPicker::calculateBgEstimation_(const RichPeakChromatogram& chromatogram, double best_left, double best_right)
  {
    // determine (in the chromatogram) the intensity at the left / right border
    RichPeakChromatogram::const_iterator it = chromatogram.begin();
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
      return 0;
    }

    // decrease the iterator and the nr_points by one (because we went one too far)
    double intensity_right = (it--)->getIntensity();
    nr_points--;

    double avg_noise_level = (intensity_right + intensity_left) / 2;
    return avg_noise_level * nr_points;
  }

  void MRMTransitionGroupPicker::findLargestPeak(std::vector<RichPeakChromatogram>& picked_chroms, int& chr_idx, int& peak_idx)
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
