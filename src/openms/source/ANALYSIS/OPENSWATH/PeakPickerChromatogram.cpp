// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerChromatogram.h>

namespace OpenMS
{
  PeakPickerChromatogram::PeakPickerChromatogram() :
    DefaultParamHandler("PeakPickerChromatogram")
  {
    // For SWATH-MS data from 5600 TripleTOF, these settings are recommended: 
    //
    // sgolay_frame_length = 9  (29.7s on our data)
    // gauss_width = 30  (if even gauss is used)
    // use_gauss = false
    //
    defaults_.setValue("sgolay_frame_length", 15, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
    defaults_.setValue("sgolay_polynomial_order", 3, "Order of the polynomial that is fitted.");
    defaults_.setValue("gauss_width", 50.0, "Gaussian width in seconds, estimated peak size.");
    defaults_.setValue("use_gauss", "true", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
    defaults_.setValidStrings("use_gauss", {"false","true"});

    defaults_.setValue("peak_width", -1.0, "Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.");
    defaults_.setValue("signal_to_noise", 1.0, "Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.");
    defaults_.setMinFloat("signal_to_noise", 0.0);

    defaults_.setValue("sn_win_len", 1000.0, "Signal to noise window length.");
    defaults_.setValue("sn_bin_count", 30, "Signal to noise bin count.");
    defaults_.setValue("write_sn_log_messages", "false", "Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin");
    defaults_.setValidStrings("write_sn_log_messages", {"true","false"});

    defaults_.setValue("remove_overlapping_peaks", "false", "Try to remove overlapping peaks during peak picking");
    defaults_.setValidStrings("remove_overlapping_peaks", {"false","true"});

    defaults_.setValue("method", "corrected", "Which method to choose for chromatographic peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed chromatogram or Crawdad on smoothed chromatogram).");
    defaults_.setValidStrings("method", {"legacy","corrected","crawdad"});

    // write defaults into Param object param_
    defaultsToParam_();
    updateMembers_();

    // PeakPickerHiRes pp_;
    Param pepi_param = pp_.getDefaults();
    pepi_param.setValue("signal_to_noise", signal_to_noise_);
    // disable spacing constraints, since we're dealing with chromatograms
    pepi_param.setValue("spacing_difference", 0.0);
    pepi_param.setValue("spacing_difference_gap", 0.0);
    pepi_param.setValue("report_FWHM", "true");
    pepi_param.setValue("report_FWHM_unit", "absolute");
    pp_.setParameters(pepi_param);
  }

  void PeakPickerChromatogram::pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom)
  {
    MSChromatogram s;
    pickChromatogram(chromatogram, picked_chrom, s);
  }
  
  void PeakPickerChromatogram::pickChromatogram(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom, MSChromatogram& smoothed_chrom)
  {
    if (!chromatogram.isSorted())
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Chromatogram must be sorted by position");
    }

    if (chromatogram.empty())
    {
      OPENMS_LOG_DEBUG << " ====  Chromatogram " << chromatogram.getNativeID() << "empty. Skip picking.";
      return;
    }
    else
    {
        OPENMS_LOG_DEBUG << " ====  Picking chromatogram " << chromatogram.getNativeID() << 
        " with " << chromatogram.size() << " peaks (start at RT " << chromatogram[0].getRT() << " to RT " << chromatogram.back().getRT() << ") "
        "using method \'" << method_ << "\'" << std::endl;
    }
    picked_chrom.clear(true);
    // Crawdad has its own methods, so we can call the wrapper directly
    if (method_ == "crawdad")
    {
      pickChromatogramCrawdad_(chromatogram, picked_chrom);
      return;
    }

    // Smooth the chromatogram
    smoothed_chrom = chromatogram;
    if (!use_gauss_)
    {
      sgolay_.filter(smoothed_chrom);
    }
    else
    {
      gauss_.filter(smoothed_chrom);
    }

    // Find initial seeds (peak picking)
    pp_.pick(smoothed_chrom, picked_chrom);
    OPENMS_LOG_DEBUG << "Picked " << picked_chrom.size() << " chromatographic peaks." << std::endl;

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
    OPENMS_POSTCONDITION(picked_chrom.getFloatDataArrays().size() == 1 &&
                         picked_chrom.getFloatDataArrays()[IDX_FWHM].getName() == "FWHM", "Swath: PeakPicking did not deliver FWHM attributes.")

    picked_chrom.getFloatDataArrays().resize(SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[IDX_ABUNDANCE].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[IDX_RIGHTBORDER].setName("rightWidth");
    // just copy FWHM from initial peak picking

    picked_chrom.getFloatDataArrays()[IDX_ABUNDANCE].reserve(picked_chrom.size());
    picked_chrom.getFloatDataArrays()[IDX_LEFTBORDER].reserve(picked_chrom.size());
    picked_chrom.getFloatDataArrays()[IDX_RIGHTBORDER].reserve(picked_chrom.size());
    for (Size i = 0; i < picked_chrom.size(); i++)
    {
      picked_chrom.getFloatDataArrays()[IDX_ABUNDANCE].push_back(integrated_intensities_[i]);
      picked_chrom.getFloatDataArrays()[IDX_LEFTBORDER].push_back((float)chromatogram[left_width_[i]].getRT());
      picked_chrom.getFloatDataArrays()[IDX_RIGHTBORDER].push_back((float)chromatogram[right_width_[i]].getRT());
    }
  }

  void PeakPickerChromatogram::pickChromatogram_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom)
  {

    integrated_intensities_.clear();
    left_width_.clear();
    right_width_.clear();
    integrated_intensities_.reserve(picked_chrom.size());
    left_width_.reserve(picked_chrom.size());
    right_width_.reserve(picked_chrom.size());

    if (signal_to_noise_ > 0.0)
    {
      snt_.init(chromatogram);
    }
    Size current_peak = 0;
    for (Size i = 0; i < picked_chrom.size(); i++)
    {
      const double central_peak_rt = picked_chrom[i].getRT();
      current_peak = findClosestPeak_(chromatogram, central_peak_rt, current_peak);
      const Size min_i = current_peak;

      // peak core found, now extend it to the left
      Size k = 2;
      while ((min_i - k + 1) > 0
             //&& std::fabs(chromatogram[min_i-k].getMZ() - peak_raw_data.begin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i - k].getIntensity() < chromatogram[min_i - k + 1].getIntensity()
               || (peak_width_ > 0.0 && std::fabs(chromatogram[min_i - k].getRT() - central_peak_rt) < peak_width_))
            && (signal_to_noise_ <= 0.0 || snt_.getSignalToNoise(min_i - k) >= signal_to_noise_))
      {
        ++k;
      }
      int left_idx = min_i - k + 1;

      // to the right
      k = 2;
      while ((min_i + k) < chromatogram.size()
             //&& std::fabs(chromatogram[min_i+k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference*min_spacing
            && (chromatogram[min_i + k].getIntensity() < chromatogram[min_i + k - 1].getIntensity()
               || (peak_width_ > 0.0 && std::fabs(chromatogram[min_i + k].getRT() - central_peak_rt) < peak_width_))
            && (signal_to_noise_ <= 0.0 || snt_.getSignalToNoise(min_i + k) >= signal_to_noise_) )
      {
        ++k;
      }
      int right_idx = min_i + k - 1;

      left_width_.push_back(left_idx);
      right_width_.push_back(right_idx);
      integrated_intensities_.push_back(0);

      OPENMS_LOG_DEBUG << "Found peak at " << central_peak_rt << " with intensity "  << picked_chrom[i].getIntensity()
                << " and borders " << chromatogram[left_width_[i]].getRT() << " " << chromatogram[right_width_[i]].getRT() <<
        " (" << chromatogram[right_width_[i]].getRT() - chromatogram[left_width_[i]].getRT() << ") "
                << 0 << " weighted RT " << /* weighted_mz << */ std::endl;
    }
  }

#ifdef WITH_CRAWDAD
  void PeakPickerChromatogram::pickChromatogramCrawdad_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom)
  {
    OPENMS_LOG_DEBUG << "Picking chromatogram using crawdad " << std::endl;

    // copy meta data of the input chromatogram
    picked_chrom.clear(true);
    picked_chrom.ChromatogramSettings::operator=(chromatogram);
    picked_chrom.MetaInfoInterface::operator=(chromatogram);
    picked_chrom.setName(chromatogram.getName());

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
    picked_chrom.getFloatDataArrays().resize(SIZE_OF_FLOATINDICES);
    picked_chrom.getFloatDataArrays()[IDX_ABUNDANCE].setName("IntegratedIntensity");
    picked_chrom.getFloatDataArrays()[IDX_LEFTBORDER].setName("leftWidth");
    picked_chrom.getFloatDataArrays()[IDX_RIGHTBORDER].setName("rightWidth");
    for (std::vector<crawpeaks::SlimCrawPeak>::iterator it = result.begin(); it != result.end(); ++it)
    {
      ChromatogramPeak p;
      p.setRT(chromatogram[it->peak_rt_idx].getRT());
      p.setIntensity(it->peak_area); //chromatogram[it->peak_rt_idx].getIntensity() );

      picked_chrom.getFloatDataArrays()[IDX_ABUNDANCE].push_back(it->peak_area);
      picked_chrom.getFloatDataArrays()[IDX_LEFTBORDER].push_back(chromatogram[it->start_rt_idx].getRT());
      picked_chrom.getFloatDataArrays()[IDX_RIGHTBORDER].push_back(chromatogram[it->stop_rt_idx].getRT());

      OPENMS_LOG_DEBUG << "Found peak at " << p.getRT() << " and "  << chromatogram[it->peak_rt_idx].getIntensity()
                << " with borders " << chromatogram[it->start_rt_idx].getRT() << " " << chromatogram[it->stop_rt_idx].getRT()  <<  " (" << chromatogram[it->start_rt_idx].getRT() - chromatogram[it->stop_rt_idx].getRT() << ") "
                << it->peak_area << " weighted RT " << /* weighted_mz << */ std::endl;

      picked_chrom.push_back(p);

    }

  }
#else
  void PeakPickerChromatogram::pickChromatogramCrawdad_(const MSChromatogram& /* chromatogram */, MSChromatogram& /* picked_chrom */)
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "PeakPickerChromatogram was not compiled with crawdad, please choose a different algorithm!");
  }
#endif

  void PeakPickerChromatogram::removeOverlappingPeaks_(const MSChromatogram& chromatogram, MSChromatogram& picked_chrom)
  {
    if (picked_chrom.empty()) {return; }
    OPENMS_LOG_DEBUG << "Remove overlapping peaks now (size " << picked_chrom.size() << ")" << std::endl;
    Size current_peak = 0;
    // Find overlapping peaks
    for (Size i = 0; i < picked_chrom.size() - 1; i++)
    {
      // Check whether the current right overlaps with the next left
      // See whether we can correct this and find some border between the two
      // features ...
      if (right_width_[i] > left_width_[i + 1])
      {
        const int current_left_idx = left_width_[i];
        const int current_right_idx = right_width_[i];
        const int next_left_idx = left_width_[i + 1];
        const int next_right_idx = right_width_[i + 1];
        OPENMS_LOG_DEBUG << " Found overlapping " << i << " : " << current_left_idx << " " << current_right_idx << std::endl;
        OPENMS_LOG_DEBUG << "                   -- with  " << i + 1 << " : " << next_left_idx << " " << next_right_idx << std::endl;

        // Find the peak width and best RT
        double central_peak_rt = picked_chrom[i].getPos();
        double next_peak_rt = picked_chrom[i + 1].getPos();
        current_peak = findClosestPeak_(chromatogram, central_peak_rt, current_peak);
        Size next_peak = findClosestPeak_(chromatogram, next_peak_rt, current_peak);

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

        OPENMS_LOG_DEBUG << "New peak l: " << chromatogram[current_left_idx].getPos() 
          << " " << chromatogram[new_right_border].getPos() 
          << " int " << integrated_intensities_[i] << std::endl;
        OPENMS_LOG_DEBUG << "New peak r: " << chromatogram[new_left_border].getPos() 
          << " " << chromatogram[next_right_idx].getPos() 
          << " int " << integrated_intensities_[i + 1] << std::endl;


        right_width_[i] = new_right_border;
        left_width_[i + 1] = new_left_border;

      }
    }
  }

  Size PeakPickerChromatogram::findClosestPeak_(const MSChromatogram& chromatogram, double target_rt, Size current_peak)
  {
    while (current_peak < chromatogram.size())
    {
      // check if we have walked past the RT of the peak
      if (target_rt < chromatogram[current_peak].getRT())
      {
        // see which one is closer, the current one or the one before
        if (current_peak > 0 &&
            std::fabs(target_rt - chromatogram[current_peak - 1].getRT()) <
            std::fabs(target_rt - chromatogram[current_peak].getRT()))
        {
          current_peak--;
        }

        return current_peak;
      }
      current_peak++;
    }
    return current_peak;
  }

  void PeakPickerChromatogram::integratePeaks_(const MSChromatogram& chromatogram)
  {
    for (Size i = 0; i < left_width_.size(); i++)
    {
      const int current_left_idx = left_width_[i];
      const int current_right_idx = right_width_[i];

      // Also integrate the intensities
      integrated_intensities_[i] = 0;
      for (int k = current_left_idx; k <= current_right_idx; k++)
      {
        integrated_intensities_[i] += chromatogram[k].getIntensity();
      }
    }
  }

  void PeakPickerChromatogram::updateMembers_()
  {
    sgolay_frame_length_ = (UInt)param_.getValue("sgolay_frame_length");
    sgolay_polynomial_order_ = (UInt)param_.getValue("sgolay_polynomial_order");
    gauss_width_ = (double)param_.getValue("gauss_width");
    peak_width_ = (double)param_.getValue("peak_width");
    signal_to_noise_ = (double)param_.getValue("signal_to_noise");
    sn_win_len_ = (double)param_.getValue("sn_win_len");
    sn_bin_count_ = (UInt)param_.getValue("sn_bin_count");
    // TODO make list, not boolean
    use_gauss_ = (bool)param_.getValue("use_gauss").toBool();
    remove_overlapping_ = (bool)param_.getValue("remove_overlapping_peaks").toBool();
    write_sn_log_messages_ = (bool)param_.getValue("write_sn_log_messages").toBool();
    method_ = (String)param_.getValue("method").toString();

    if (method_ != "crawdad" && method_ != "corrected" && method_ != "legacy")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "Method needs to be one of: crawdad, corrected, legacy");
    }

    Param sg_filter_parameters = sgolay_.getParameters();
    sg_filter_parameters.setValue("frame_length", sgolay_frame_length_);
    sg_filter_parameters.setValue("polynomial_order", sgolay_polynomial_order_);
    sgolay_.setParameters(sg_filter_parameters);

    Param gfilter_parameters = gauss_.getParameters();
    gfilter_parameters.setValue("gaussian_width", gauss_width_);
    gauss_.setParameters(gfilter_parameters);

    Param snt_parameters = snt_.getParameters();
    snt_parameters.setValue("win_len", sn_win_len_);
    snt_parameters.setValue("bin_count", sn_bin_count_);
    snt_parameters.setValue("write_log_messages", param_.getValue("write_sn_log_messages"));
    snt_.setParameters(snt_parameters);

#ifndef WITH_CRAWDAD
    if (method_ == "crawdad")
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                       "PeakPickerChromatogram was not compiled with crawdad, please choose a different algorithm!");
    }
#endif
  }

}
