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
    // NEW default settings: recommeded:
    //
    // sgolay_frame_length = 9  (29.7s on our data)
    // gauss_width = 30  (if even gauss is used)
    // use_gauss = false
    //
    //
    // THIS is the most important change !! this caused a lot of trouble ...
    // peak_width = -1 (do not force a certain width)
    // method = corrected
    //
    defaults_.setValue("sgolay_frame_length", 15, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("sgolay_polynomial_order", 3, "Order or the polynomial that is fitted.");
    defaults_.setValue("gauss_width", 50.0, "Gaussian width in seconds, estimated peak size.");
    defaults_.setValue("use_gauss", "true", "Use gauss for smoothing (other option is sgolay)");

    defaults_.setValue("peak_width", 40.0, "Estimated peak width in seconds.");
    defaults_.setValue("signal_to_noise", 1.0, "Signal to noise.");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("sn_win_len", 1000.0, "Signal to noise window length.");
    defaults_.setValue("sn_bin_count", 30, "Signal to noise bin count.");

    defaults_.setValue("remove_overlapping_peaks", "false", "Try to remove overlappign peaks during peak picking");

    defaults_.setValue("method", "legacy", "Which method to choose for chromatographic peak-picking (OpenSWATH legacy, corrected picking or Crawdad)");
    defaults_.setValidStrings("method", StringList::create("legacy,corrected,crawdad"));

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();
  }

  void PeakPickerMRM::pickChromatogram(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
  {
    if (!chromatogram.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Chromatogram must be sorted by position");
    }

    LOG_DEBUG << " ====  Picking chromatogram " << chromatogram.getNativeID() << std::endl;
    picked_chrom.clear(true);

    // Crowdad has its own methods, so we can call the wrapper directly
    if (method_ == "crawdad")
    {
      pickChromatogramCrowdad(chromatogram, picked_chrom);
      return;
    }

    // Smooth the chromatogram
    RichPeakChromatogram smoothed_chrom = chromatogram;
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

    if (method_ == "legacy")
    {
      // Legacy is to use the original chromatogram for peak-detection
      pickChromatogram_(chromatogram, picked_chrom);
      if (remove_overlapping_)
        removeOverlappingPeaks_(chromatogram, picked_chrom);

      // for peak integration, we want to use the raw data
      integratePeaks_(chromatogram);
    }
    else if (method_ == "corrected")
    {
      // use the smoothed chromatogram to derive the peak boundaries
      pickChromatogram_(smoothed_chrom, picked_chrom);
      if (remove_overlapping_)
        removeOverlappingPeaks_(smoothed_chrom, picked_chrom);

      // for peak integration, we want to use the raw data
      integratePeaks_(chromatogram);
    }

    // Store the result in the picked_chromatogram
    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
    for (Size i = 0; i < picked_chrom.size(); i++)
    {
      Real leftborder = chromatogram[left_width[i]].getMZ();
      Real rightborder = chromatogram[right_width[i]].getMZ();
      picked_chrom.getFloatDataArrays()[0].push_back(integrated_intensities[i]);
      picked_chrom.getFloatDataArrays()[1].push_back(leftborder);
      picked_chrom.getFloatDataArrays()[2].push_back(rightborder);
    }
  }

  void PeakPickerMRM::pickChromatogram_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
  {
    SignalToNoiseEstimatorMedian<RichPeakChromatogram> snt;
    Param snt_parameters = snt.getParameters();
    snt_parameters.setValue("win_len", sn_win_len_);
    snt_parameters.setValue("bin_count", sn_bin_count_);
    snt.setParameters(snt_parameters);

    integrated_intensities.clear();
    left_width.clear();
    right_width.clear();
    integrated_intensities.reserve(picked_chrom.size());
    left_width.reserve(picked_chrom.size());
    right_width.reserve(picked_chrom.size());

    snt.init(chromatogram);
    Size current_peak = 0;
    for (Size i = 0; i < picked_chrom.size(); i++)
    {
      const double central_peak_mz = picked_chrom[i].getMZ();
      current_peak = findClosestPeak_(chromatogram, central_peak_mz, current_peak);
      const Size min_i = current_peak;

      // peak core found, now extend it to the left
      Size k = 2;
      while ((min_i - k + 1) > 0
             //&& std::fabs(chromatogram[min_i-k].getMZ() - peak_raw_data.begin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i - k].getIntensity() < chromatogram[min_i - k + 1].getIntensity()
               || (peak_width_ > 0.0 && std::fabs(chromatogram[min_i - k].getMZ() - central_peak_mz) < peak_width_)
                )
            && snt.getSignalToNoise(chromatogram[min_i - k]) >= signal_to_noise_)
      {
        ++k;
      }
      int left_idx = min_i - k + 1;

      // to the right
      k = 2;
      while ((min_i + k) < chromatogram.size()
             //&& std::fabs(chromatogram[min_i+k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i + k].getIntensity() < chromatogram[min_i + k - 1].getIntensity()
               || (peak_width_ > 0.0 && std::fabs(chromatogram[min_i + k].getMZ() - central_peak_mz) < peak_width_)
                )
            && snt.getSignalToNoise(chromatogram[min_i + k]) >= signal_to_noise_)
      {
        ++k;
      }
      int right_idx = min_i + k - 1;

      left_width.push_back(left_idx);
      right_width.push_back(right_idx);
      integrated_intensities.push_back(0);

      LOG_DEBUG << "Found peak at " << central_peak_mz << " and "  << picked_chrom[i].getIntensity()
                << " with borders " << chromatogram[left_width[i]].getMZ() << " " << chromatogram[right_width[i]].getMZ() <<
        " (" << chromatogram[right_width[i]].getMZ() - chromatogram[left_width[i]].getMZ() << ") "
                << 0 << " weighted RT " << /* weighted_mz << */ std::endl;
    }
  }

#ifdef WITH_CRAWDAD
  void PeakPickerMRM::pickChromatogramCrowdad(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
  {
    LOG_DEBUG << "Picking chromatogram using crawdad " << std::endl;

    std::vector<double> time;
    std::vector<double> intensity;
    for (Size i = 0; i < chromatogram.size(); i++)
    {
      time.push_back(chromatogram[i].getRT());
      intensity.push_back(chromatogram[i].getIntensity());
    }

    CrawdadWrapper crawdad_pp;
    crawdad_pp.SetChromatogram(time, intensity);
    std::vector<crawpeaks::SlimCrawPeak> result = crawdad_pp.CalcPeaks();

    picked_chrom.getFloatDataArrays().clear();
    picked_chrom.getFloatDataArrays().resize(3);
    picked_chrom.getFloatDataArrays()[0].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[1].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[2].setName("rightWidth");
    for (std::vector<crawpeaks::SlimCrawPeak>::iterator it = result.begin(); it != result.end(); it++)
    {
      ChromatogramPeak p;
      p.setRT(chromatogram[it->peak_rt_idx].getRT());
      p.setIntensity(it->peak_area); //chromatogram[it->peak_rt_idx].getIntensity() );

      picked_chrom.getFloatDataArrays()[0].push_back(it->peak_area);
      picked_chrom.getFloatDataArrays()[1].push_back(chromatogram[it->start_rt_idx].getRT());
      picked_chrom.getFloatDataArrays()[2].push_back(chromatogram[it->stop_rt_idx].getRT());
      /*
      int peak_rt_idx, start_rt_idx, stop_rt_idx, max_rt_idx;
      int mz_idx;
      int len;
      float fwhm;
      bool fwhm_calculated_ok;
      float bg_area;
      float raw_area; // total area under the curve, including background
      float peak_area;
      float bgslope;
      ///cutoff level for extending past the peak

      ///maximum height, calculated above background
      float peak_height;
      float raw_height;

      */

      LOG_DEBUG << "Found peak at " << p.getRT() << " and "  << chromatogram[it->peak_rt_idx].getIntensity()
                << " with borders " << chromatogram[it->start_rt_idx].getRT() << " " << chromatogram[it->stop_rt_idx].getRT()  <<  " (" << chromatogram[it->start_rt_idx].getRT() - chromatogram[it->stop_rt_idx].getRT() << ") "
                << it->peak_area << " weighted RT " << /* weighted_mz << */ std::endl;

      picked_chrom.push_back(p);

    }

  }

#else
  void PeakPickerMRM::pickChromatogramCrowdad(const RichPeakChromatogram& /* chromatogram */, RichPeakChromatogram& /* picked_chrom */)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                     "PeakPickerMRM was not compiled with crawdad, please choose a different algorithm!");
  }

#endif


  void PeakPickerMRM::removeOverlappingPeaks_(const RichPeakChromatogram& chromatogram, RichPeakChromatogram& picked_chrom)
  {
    if (picked_chrom.empty()) {return; }
    LOG_DEBUG << "Remove overlapping peaks now (size " << picked_chrom.size() << ")" << std::endl;
    Size current_peak = 0;
    // Find overlapping peaks
    for (Size i = 0; i < picked_chrom.size() - 1; i++)
    {
      // Check whether the current right overlaps with the next left
      // See whether we can correct this and find some border between the two
      // features ...
      if (right_width[i] > left_width[i + 1])
      {
        const int current_left_idx = left_width[i];
        const int current_right_idx = right_width[i];
        const int next_left_idx = left_width[i + 1];
        const int next_right_idx = right_width[i + 1];
        LOG_DEBUG << " Found overlapping " << i << " : " << current_left_idx << " " << current_right_idx << std::endl;
        LOG_DEBUG << "                   -- with  " << i + 1 << " : " << next_left_idx << " " << next_right_idx << std::endl;

        // Find the peak width and best RT
        double central_peak_mz = picked_chrom[i].getMZ();
        double next_peak_mz = picked_chrom[i + 1].getMZ();
        current_peak = findClosestPeak_(chromatogram, central_peak_mz, current_peak);
        Size next_peak = findClosestPeak_(chromatogram, next_peak_mz, current_peak);

        // adjust the right border of the current and left border of next
        Size k = 1;
        while ((current_peak + k) < chromatogram.size()
              && (chromatogram[current_peak + k].getIntensity() < chromatogram[current_peak + k - 1].getIntensity()))
        {
          ++k;
        }
        Size new_right_border = current_peak + k - 1;
        k = 1;
        while ((next_peak - k + 1) > 0
              && (chromatogram[next_peak - k].getIntensity() < chromatogram[next_peak - k + 1].getIntensity()))
        {
          ++k;
        }
        Size new_left_border = next_peak - k + 1;

        // assert that the peaks are now not overlapping any more ...
        if (new_left_border < new_right_border)
        {
          std::cerr << "Something went wrong, peaks are still overlapping!" << " - new left border " << new_left_border << " vs " << new_right_border << " -- will take the mean" << std::endl;
          new_left_border = (new_left_border + new_right_border) / 2;
          new_right_border = (new_left_border + new_right_border) / 2;

        }

        LOG_DEBUG << "New peak l: " << chromatogram[current_left_idx].getMZ() << " " << chromatogram[new_right_border].getMZ() << " int " << integrated_intensities[i] << std::endl;
        LOG_DEBUG << "New peak r: " << chromatogram[new_left_border].getMZ() << " " << chromatogram[next_right_idx].getMZ() << " int " << integrated_intensities[i + 1] << std::endl;


        right_width[i] = new_right_border;
        left_width[i + 1] = new_left_border;

      }
    }
  }

  Size PeakPickerMRM::findClosestPeak_(const RichPeakChromatogram& chromatogram, double central_peak_mz, Size current_peak)
  {
    while (current_peak < chromatogram.size())
    {
      // check if we have walked past the RT of the peak
      if (central_peak_mz - chromatogram[current_peak].getMZ() < 0.0)
      {
        // see which one is closer, the current one or the one before
        if (current_peak > 0 &&
            std::fabs(central_peak_mz - chromatogram[current_peak - 1].getMZ()) <
            std::fabs(central_peak_mz - chromatogram[current_peak].getMZ()))
        {
          current_peak--;
        }

        return current_peak;
      }
      current_peak++;
    }
    return current_peak;
  }

  void PeakPickerMRM::integratePeaks_(const RichPeakChromatogram& chromatogram)
  {
    for (Size i = 0; i < left_width.size(); i++)
    {
      const int current_left_idx = left_width[i];
      const int current_right_idx = right_width[i];

      // Also integrate the intensities
      integrated_intensities[i] = 0;
      for (Size k = current_left_idx; k <= current_right_idx; k++)
      {
        integrated_intensities[i] += chromatogram[k].getIntensity();
      }
    }
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
    // TODO make list, not boolean
    use_gauss_ = (bool)param_.getValue("use_gauss").toBool();
    remove_overlapping_ = (bool)param_.getValue("remove_overlapping_peaks").toBool();
    method_ = (String)param_.getValue("method");

    if (method_ != "crawdad" && method_ != "corrected" && method_ != "legacy")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "Method needs to be one of: crawdad, corrected, legacy");
    }

#ifndef WITH_CRAWDAD
    if (method_ == "crawdad")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __PRETTY_FUNCTION__,
                                       "PeakPickerMRM was not compiled with crawdad, please choose a different algorithm!");
    }
#endif
  }

}
