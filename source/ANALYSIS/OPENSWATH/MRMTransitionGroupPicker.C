// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

namespace OpenMS
{

  MRMTransitionGroupPicker::MRMTransitionGroupPicker() :
    DefaultParamHandler("MRMTransitionGroupPicker")
  {
    defaults_.setValue("sgolay_frame_length", 15, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("sgolay_polynomial_order", 3, "Order or the polynomial that is fitted.");
    defaults_.setValue("gauss_width", 50.0, "Gaussian width in seconds, estimated peak size.");
    defaults_.setValue("use_gauss", "true", "Use gauss for smoothing (other option is sgolay)");

    defaults_.setValue("peak_width", 40.0, "Estimated peak width in seconds.");
    defaults_.setMinFloat("peak_width", 0.0);
    defaults_.setValue("signal_to_noise", 1.0, "Signal to noise.");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("sn_win_len", 1000.0, "Signal to noise window length.");
    defaults_.setValue("sn_bin_count", 30, "Signal to noise bin count.");

    defaults_.setValue("stop_after_feature", -1, "Stop finding after feature (ordered by intensity; -1 means do not stop).");
    defaults_.setValue("stop_after_intensity_ratio", 0.0001, "Stop after reaching intensity ratio");
    defaults_.setValue("stop_report_after_feature", -1, "Stop reporting after feature (ordered by quality; 1 means do not stop).");

    defaults_.setValue("background_subtraction", "none", "Try to apply a background subtraction to the peak (experimental). The background is estimated at the peak boundaries, either the smoothed or the raw chromatogram data can be used for that."); //, StringList::create("advanced"));
    defaults_.setValidStrings("background_subtraction", StringList::create("none,smoothed,original"));

    // write defaults into Param object param_
    defaultsToParam_();
    handle_params();
  }

  MRMTransitionGroupPicker::~MRMTransitionGroupPicker()
  {
  }

  MRMTransitionGroupPicker & MRMTransitionGroupPicker::operator = (const MRMTransitionGroupPicker &rhs)
  {
    if (&rhs == this)
      return *this;

    // dont copy parameters

    return *this;
  }

  void MRMTransitionGroupPicker::updateMembers_()
  {
    handle_params();
  }

  void MRMTransitionGroupPicker::handle_params()
  {
    sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
    sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
    gauss_width_ = (DoubleReal)param_.getValue("gauss_width");
    peak_width_ = (DoubleReal)param_.getValue("peak_width");
    signal_to_noise_ = (DoubleReal)param_.getValue("signal_to_noise");
    sn_win_len_ = (DoubleReal)param_.getValue("sn_win_len");
    sn_bin_count_ = (UInt)param_.getValue("sn_bin_count");
    use_gauss_ = (bool)param_.getValue("use_gauss").toBool();

    stop_after_feature_ = (int)param_.getValue("stop_after_feature");
    stop_after_intensity_ratio_ = (DoubleReal)param_.getValue("stop_after_intensity_ratio");

    background_subtraction_ = param_.getValue("background_subtraction");
  }

  void MRMTransitionGroupPicker::pickChromatogram(const RichPeakChromatogram & chromatogram, RichPeakChromatogram & smoothed_chrom, RichPeakChromatogram & picked_chrom)
  {
    // Smooth the chromatogram
    smoothed_chrom = chromatogram;
    if (!use_gauss_)
    {
      SavitzkyGolayFilter sgolay;
      Param filter_parameters = sgolay.getParameters();
      filter_parameters.setValue("frame_length", sgolay_frame_length_);
      filter_parameters.setValue("polynomial_order", sgolay_polynomial_order_);
      sgolay.setParameters(filter_parameters);
      sgolay.filter(smoothed_chrom);
    }
    else
    {
      GaussFilter gauss;
      Param filter_parameters = gauss.getParameters();
      filter_parameters.setValue("gaussian_width", gauss_width_);
      gauss.setParameters(filter_parameters);
      gauss.filter(smoothed_chrom);
    }

    // Find initial seeds (peak picking)
    PeakPickerHiRes pp;
    Param pepi_param = PeakPickerHiRes().getDefaults();
    pp.setParameters(pepi_param);
    pp.pick(smoothed_chrom, picked_chrom);

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");

    SignalToNoiseEstimatorMedian<RichPeakChromatogram> snt;
    Param snt_parameters = snt.getParameters();
    snt_parameters.setValue("win_len", sn_win_len_);
    snt_parameters.setValue("bin_count", sn_bin_count_);
    snt.setParameters(snt_parameters);

    snt.init(chromatogram);
    for (Size i = 0; i < picked_chrom.size(); i++)
    {

      // Find the peak width and best RT
      // FEATURE : we could use #include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h> 
      // --> this sounds costly, fitting a spline and then linear regression ...
      double central_peak_mz = picked_chrom[i].getMZ();
      double min_d = std::fabs(central_peak_mz - chromatogram[0].getMZ());
      Size min_i = 0;
      for (Size j = 0; j < chromatogram.size(); j++)
      {
        if (std::fabs(central_peak_mz - chromatogram[j].getMZ()) < min_d)
        {
          min_d = std::fabs(central_peak_mz - chromatogram[j].getMZ());
          min_i = j;
        }
      }

      std::map<double, double> peak_raw_data;
      peak_raw_data[chromatogram[min_i].getMZ()] = chromatogram[min_i].getIntensity();

      // peak core found, now extend it to the left
      Size k = 1;
      while ((min_i - k + 1) > 0
             //&& std::fabs(chromatogram[min_i-k].getMZ() - peak_raw_data.begin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i - k].getIntensity() < peak_raw_data.begin()->second
               || std::fabs(chromatogram[min_i - k].getMZ() - central_peak_mz) < peak_width_)
            && snt.getSignalToNoise(chromatogram[min_i - k]) >= signal_to_noise_)
      {
        peak_raw_data[chromatogram[min_i - k].getMZ()] = chromatogram[min_i - k].getIntensity();
        ++k;
      }
      Real leftborder = chromatogram[min_i - k + 1].getMZ();

      // to the right
      k = 1;
      while ((min_i + k) < chromatogram.size()
             //&& std::fabs(chromatogram[min_i+k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i + k].getIntensity() < peak_raw_data.rbegin()->second
               || std::fabs(chromatogram[min_i + k].getMZ() - central_peak_mz) < peak_width_)
            && snt.getSignalToNoise(chromatogram[min_i + k]) >= signal_to_noise_)
      {
        peak_raw_data[chromatogram[min_i + k].getMZ()] = chromatogram[min_i + k].getIntensity();
        ++k;
      }

      Real rightborder = chromatogram[min_i + k - 1].getMZ();

      double weighted_mz = 0;
      Real integrated_intensity = 0;
      for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it)
      {
        weighted_mz += map_it->first * map_it->second;
        integrated_intensity += map_it->second;
      }
      // For peaks in chromatograms, it does not make much sense to recenter
      // them since we are not concerned with accuracy as with m/z peaks but
      // rather want to see coelution. For this relative peak shape is much
      // more interesting.
      weighted_mz /= integrated_intensity;

#ifdef DEBUG_MRMPEAKPICKER
      double central_peak_int = picked_chrom[i].getIntensity();
      std::cout << "Found peak at " << central_peak_mz << " and "  << central_peak_int << " with borders " << leftborder << " " << rightborder <<  " (" << rightborder - leftborder << ") " << integrated_intensity << " weighted RT " << weighted_mz << std::endl;
#endif


      picked_chrom.getFloatDataArrays()[0].push_back(integrated_intensity);
      picked_chrom.getFloatDataArrays()[1].push_back(leftborder);
      picked_chrom.getFloatDataArrays()[2].push_back(rightborder);

      /*
      picked_chrom.getFloatDataArrays()[0].push_back( integrated_intensity );
      picked_chrom.getFloatDataArrays()[1].push_back( picked_chrom[i].getMZ() - peak_width / 2.0);
      picked_chrom.getFloatDataArrays()[2].push_back( picked_chrom[i].getMZ() + peak_width / 2.0);
      */
    }

  }

  void MRMTransitionGroupPicker::findLargestPeak(std::vector<RichPeakChromatogram> & picked_chroms, int & chr_idx, int & peak_idx)
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
          chr_idx = k;
          peak_idx = i;
        }
      }
    }
  }

}
