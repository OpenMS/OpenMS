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

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMRM.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/FILTERING/SMOOTHING/GaussFilter.h>

namespace OpenMS
{
  PeakPickerMRM::PeakPickerMRM() :
      DefaultParamHandler("PeakPickerMRM")
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

    defaults_.setValue("remove_overlapping_peaks", "false", "Try to remove overlappign peaks during peak picking");

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void PeakPickerMRM::pickChromatogram(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& smoothed_chrom, RichPeakChromatogram& picked_chrom)
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
    LOG_DEBUG << " ====  Picking chromatogram " << chromatogram.getNativeID() << std::endl;
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

      // double weighted_mz = 0; // not used
      Real integrated_intensity = 0;
      for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it)
      {
        // weighted_mz += map_it->first * map_it->second;
        integrated_intensity += map_it->second;
      }
      // For peaks in chromatograms, it does not make much sense to recenter
      // them since we are not concerned with accuracy as with m/z peaks but
      // rather want to see coelution. For this relative peak shape is much
      // more interesting.
      // weighted_mz /= integrated_intensity;

      LOG_DEBUG << "Found peak at " << central_peak_mz << " and "  << picked_chrom[i].getIntensity() 
        << " with borders " << leftborder << " " << rightborder <<  " (" << rightborder - leftborder << ") " 
        << integrated_intensity << " weighted RT " << /* weighted_mz << */ std::endl;

      picked_chrom.getFloatDataArrays()[0].push_back(integrated_intensity);
      picked_chrom.getFloatDataArrays()[1].push_back(leftborder);
      picked_chrom.getFloatDataArrays()[2].push_back(rightborder);
    }

    if (remove_overlapping_)
      removeOverlappingPeaks_(chromatogram, picked_chrom);
  }

  void PeakPickerMRM::removeOverlappingPeaks_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
  {
    LOG_DEBUG  << "Remove overlapping peaks now" << std::endl;
    Size current_peak = 0;
    // Find overlapping peaks
    for (Size i = 0; i < picked_chrom.size() - 1; i++)
    {
      const double current_left   = picked_chrom.getFloatDataArrays()[1][i];
      const double current_right  = picked_chrom.getFloatDataArrays()[2][i];
      const double next_left   = picked_chrom.getFloatDataArrays()[1][i+1];
      const double next_right  = picked_chrom.getFloatDataArrays()[2][i+1];
      if ( current_right > next_left)
      {
        LOG_DEBUG << " Found overlapping " << i << " : " << current_left << " " << current_right << std::endl;
        LOG_DEBUG << "                   -- with  " << i +1 << " : " << next_left << " " << next_right << std::endl;
        // See whether we can correct this and find some border between the two features ... 

        // Find the peak width and best RT
        double central_peak_mz = picked_chrom[i].getMZ();
        double next_peak_mz = picked_chrom[i+1].getMZ();

        current_peak = findClosestPeak_(chromatogram, central_peak_mz, current_peak);
        Size next_peak = findClosestPeak_(chromatogram, next_peak_mz, current_peak);

        // adjust the right border of the current and left border of next
        
        // extend current peak to the right
        Size k = 1;
        while ((current_peak + k) < chromatogram.size()
              && (chromatogram[current_peak + k].getIntensity() < chromatogram[current_peak + k - 1].getIntensity() ))
        {
          ++k;
        }
        Size new_right_border = current_peak + k - 1;

        // extend next peak to the left
        k = 1;
        while ((next_peak - k + 1) > 0
              && (chromatogram[next_peak - k].getIntensity() < chromatogram[next_peak - k + 1].getIntensity() ))
        {
          ++k;
        }
        Size new_left_border = next_peak - k + 1;

        // check that the peaks are now not overlapping any more ... 
        if ( new_left_border < new_right_border) 
          {throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
              "Something went wrong, peaks are still overlapping!");}

        // set new right/left borders
        picked_chrom.getFloatDataArrays()[2][i] = chromatogram[new_right_border].getMZ();
        picked_chrom.getFloatDataArrays()[1][i+1] = chromatogram[new_left_border].getMZ();
      }
    }
  }

  Size PeakPickerMRM::findClosestPeak_(const RichPeakChromatogram& chromatogram, double central_peak_mz, Size current_peak)
  {
    while( current_peak < chromatogram.size() )
    {
      // check if we have walked past the RT of the peak
      if (central_peak_mz - chromatogram[current_peak].getMZ() < 0.0)
      {
        // see which one is closer, the current one or the one before
        if (current_peak > 0 && 
            std::fabs(central_peak_mz - chromatogram[current_peak-1].getMZ()) < 
            std::fabs(central_peak_mz - chromatogram[current_peak].getMZ()) )
        {
          current_peak--;
        }

        return current_peak;
      }
      current_peak++;
    }
    return current_peak;
  }

  void PeakPickerMRM::updateMembers_()
  {
    sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
    sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
    gauss_width_ = (DoubleReal)param_.getValue("gauss_width");
    peak_width_ = (DoubleReal)param_.getValue("peak_width");
    signal_to_noise_ = (DoubleReal)param_.getValue("signal_to_noise");
    sn_win_len_ = (DoubleReal)param_.getValue("sn_win_len");
    sn_bin_count_ = (UInt)param_.getValue("sn_bin_count");
    use_gauss_ = (bool)param_.getValue("use_gauss").toBool();
    remove_overlapping_ = (bool)param_.getValue("remove_overlapping_peaks").toBool();
  }
}
