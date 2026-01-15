// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/PeakPickerMobilogram.h>
#include <iostream>
#include <iomanip>    // For std::setw
#include <vector>
#include <algorithm>  // For std::min_element and std::max_element
#include <fstream>
#include <iostream>
#include <random>
#include <sstream>

namespace OpenMS
{
    PeakPickerMobilogram::PeakPickerMobilogram() :
            DefaultParamHandler("PeakPickerMobilogram")
//            ProgressLogger()
    {
        defaults_.setValue("sgolay_frame_length", 9, "The number of subsequent data points used for smoothing.\nThis number has to be uneven. If it is not, 1 will be added.");
        defaults_.setValue("sgolay_polynomial_order", 3, "Order of the polynomial that is fitted.");
        defaults_.setValue("gauss_width", 0.002, "Gaussian width in seconds, estimated peak size.");
        defaults_.setValue("use_gauss", "false", "Use Gaussian filter for smoothing (alternative is Savitzky-Golay filter)");
        defaults_.setValidStrings("use_gauss", {"false","true"});

        defaults_.setValue("peak_width", -1.0, "Force a certain minimal peak_width on the data (e.g. extend the peak at least by this amount on both sides) in seconds. -1 turns this feature off.");
        defaults_.setValue("signal_to_noise", 1.0, "Signal-to-noise threshold at which a peak will not be extended any more. Note that setting this too high (e.g. 1.0) can lead to peaks whose flanks are not fully captured.");
        defaults_.setMinFloat("signal_to_noise", 0.0);

        defaults_.setValue("sn_win_len", 1, "Signal to noise window length.");
        defaults_.setValue("sn_bin_count", 4, "Signal to noise bin count.");
        defaults_.setValue("write_sn_log_messages", "false", "Write out log messages of the signal-to-noise estimator in case of sparse windows or median in rightmost histogram bin");
        defaults_.setValidStrings("write_sn_log_messages", {"true","false"});

        defaults_.setValue("remove_overlapping_peaks", "false", "Try to remove overlapping peaks during peak picking");
        defaults_.setValidStrings("remove_overlapping_peaks", {"false","true"});

        defaults_.setValue("method", "corrected", "Which method to choose for mobilogram peak-picking (OpenSWATH legacy on raw data, corrected picking on smoothed mobilogram).");
        defaults_.setValidStrings("method", {"legacy","corrected","crawdad"});

        // write defaults into Param object param_
        defaultsToParam_();
        updateMembers_();

        // PeakPickerHiRes pp_;
        Param pepi_param = pp_.getDefaults();
        pepi_param.setValue("signal_to_noise", signal_to_noise_);
        // disable spacing constraints, since we're dealing with mobilograms
        pepi_param.setValue("spacing_difference", 0.0);
        pepi_param.setValue("spacing_difference_gap", 0.0);
        pepi_param.setValue("report_FWHM", "true");
        pepi_param.setValue("report_FWHM_unit", "absolute");
        pp_.setParameters(pepi_param);
    }

    void PeakPickerMobilogram::pickMobilogram(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {
      Mobilogram s;
      pickMobilogram(mobilogram, picked_mobilogram, s);
    }

    void PeakPickerMobilogram::pickMobilogram(Mobilogram mobilogram, Mobilogram& picked_mobilogram, Mobilogram& smoothed_mobilogram)
    {
      if (!mobilogram.isSorted())
      {
        throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                         "Mobilogram needs to be sorted by position.");
      }

      if (mobilogram.empty())
      {
          OPENMS_LOG_DEBUG << " ====  Mobilogram empty. Skip picking.";
          return;
      }
      else
      {
        OPENMS_LOG_DEBUG << " ====  Picking mobilogram with " << mobilogram.size() << " peaks (start at IM " << mobilogram[0].getMobility() << " to IM " << mobilogram.back().getMobility() << ") "
        "using method \'" << method_ << "\'" << std::endl;
      }
      picked_mobilogram.clear();

      // Smooth the mobilogram
      smoothed_mobilogram = mobilogram;
      if (!use_gauss_)
      {
          sgolay_.filter(smoothed_mobilogram);
      }
      else
      {
          gauss_.filter(smoothed_mobilogram);
      }

      // Find initial seeds (peak picking)
      pp_.pick(smoothed_mobilogram, picked_mobilogram);

      if (method_ == "legacy")
      {
        // Legacy is to use the original chromatogram for peak-detection
        pickMobilogram_(mobilogram, picked_mobilogram);
        if (remove_overlapping_)
          removeOverlappingPeaks_(mobilogram, picked_mobilogram);

        // for peak integration, we want to use the raw data
        integratePeaks_(mobilogram);
      }
      else if (method_ == "corrected")
      {
        // use the smoothed chromatogram to derive the peak boundaries
        pickMobilogram_(smoothed_mobilogram, picked_mobilogram);
        if (remove_overlapping_)
          removeOverlappingPeaks_(smoothed_mobilogram, picked_mobilogram);

        // for peak integration, we want to use the raw data
        integratePeaks_(mobilogram);
      }

      picked_mobilogram.getFloatDataArrays().resize(SIZE_OF_FLOATINDICES);
      picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].setName("IntegratedIntensity");
      picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].setName("leftWidth");
      picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].setName("rightWidth");
      picked_mobilogram.getFloatDataArrays()[IDX_OF_LEFTBORDER_IDX].setName("leftWidthIndex");
      picked_mobilogram.getFloatDataArrays()[IDX_OF_RIGHTBORDER_IDX].setName("rightWidthIndex");
      // just copy FWHM from initial peak picking
      picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_OF_LEFTBORDER_IDX].reserve(picked_mobilogram.size());
      picked_mobilogram.getFloatDataArrays()[IDX_OF_RIGHTBORDER_IDX].reserve(picked_mobilogram.size());
      for (Size i = 0; i < picked_mobilogram.size(); i++)
      {
        picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE].push_back(integrated_intensities_[i]);
        picked_mobilogram.getFloatDataArrays()[IDX_LEFTBORDER].push_back((float)mobilogram[left_width_[i]].getMobility());
        picked_mobilogram.getFloatDataArrays()[IDX_RIGHTBORDER].push_back((float)mobilogram[right_width_[i]].getMobility());
        picked_mobilogram.getFloatDataArrays()[IDX_OF_LEFTBORDER_IDX].push_back(left_width_[i]);
        picked_mobilogram.getFloatDataArrays()[IDX_OF_RIGHTBORDER_IDX].push_back(right_width_[i]);
      }

    }

    void PeakPickerMobilogram::filterTopPeak(Mobilogram& picked_mobilogram, std::vector<Mobilogram>& mobilograms, PeakPickerMobilogram::PeakPositions& peak_pos)
    {
      peak_pos = filterTopPeak_(picked_mobilogram, mobilograms);
    }

    void PeakPickerMobilogram::filterTopPeak(Mobilogram& picked_mobilogram, Mobilogram& mobilogram, PeakPickerMobilogram::PeakPositions& peak_pos)
    {
      peak_pos = filterTopPeak_(picked_mobilogram, mobilogram);
    }

    void PeakPickerMobilogram::pickMobilogram_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {

        integrated_intensities_.clear();
        left_width_.clear();
        right_width_.clear();
        integrated_intensities_.reserve(picked_mobilogram.size());
        left_width_.reserve(picked_mobilogram.size());
        right_width_.reserve(picked_mobilogram.size());

        Size current_peak = 0;
        for (Size i = 0; i < picked_mobilogram.size(); i++)
        {
          const double central_peak_im = picked_mobilogram[i].getMobility();
          current_peak = findClosestPeak_(mobilogram, central_peak_im, current_peak);
          const Size min_i = current_peak;

          // peak core found, now extend it to the left
          Size k = 2;
          while ((min_i - k + 1) > 0
                 && (mobilogram[min_i - k].getIntensity() < mobilogram[min_i - k + 1].getIntensity()
                     || (peak_width_ > 0.0 && std::fabs(mobilogram[min_i - k].getMobility() - central_peak_im) < peak_width_)))
          {
            ++k;
          }
          int left_idx = min_i - k + 1;

          // to the right
          k = 2;
          while ((min_i + k) < mobilogram.size()
                 && (mobilogram[min_i + k].getIntensity() < mobilogram[min_i + k - 1].getIntensity()
                     || (peak_width_ > 0.0 && std::fabs(mobilogram[min_i + k].getMobility() - central_peak_im) < peak_width_)))
          {
            ++k;
          }
          int right_idx = min_i + k - 1;

          left_width_.push_back(left_idx);
          right_width_.push_back(right_idx);
          integrated_intensities_.push_back(0);

          OPENMS_LOG_DEBUG << "Found peak at " << central_peak_im << " with intensity "  << picked_mobilogram[i].getIntensity()
                           << " and borders " << mobilogram[left_width_[i]].getMobility() << " " << mobilogram[right_width_[i]].getMobility() <<
            " (" << mobilogram[right_width_[i]].getMobility() - mobilogram[left_width_[i]].getMobility() << ") "
                           << 0 << " weighted IM " << /* weighted_mz << */ std::endl;
        }
    }

    void PeakPickerMobilogram::integratePeaks_(const Mobilogram& mobilogram)
    {
      for (Size i = 0; i < left_width_.size(); i++)
      {
        const int current_left_idx = left_width_[i];
        const int current_right_idx = right_width_[i];

        // Also integrate the intensities
        integrated_intensities_[i] = 0;
        for (int k = current_left_idx; k <= current_right_idx; k++)
        {
          integrated_intensities_[i] += mobilogram[k].getIntensity();
        }
      }
    }

    Size PeakPickerMobilogram::findClosestPeak_(const Mobilogram& mobilogram, double target_im, Size current_peak)
    {
      while (current_peak < mobilogram.size())
      {
        // check if we have walked past the IM of the peak
        if (target_im < mobilogram[current_peak].getMobility())
        {
          // see which one is closer, the current one or the one before
          if (current_peak > 0 &&
              std::fabs(target_im - mobilogram[current_peak - 1].getMobility()) <
                std::fabs(target_im - mobilogram[current_peak].getMobility()))
          {
            current_peak--;
          }

          return current_peak;
        }
        current_peak++;
      }
      return current_peak;
    }

    PeakPickerMobilogram::PeakPositions PeakPickerMobilogram::findHighestPeak_(const std::vector<double> intensities,
                                                                              const std::vector<Size> left_widths,
                                                                              const std::vector<Size> right_widths,
                                                                              const size_t im_size)
    {
      // If no peaks were found, return a peak at the center of the mobilogram
      if (intensities.empty())
      {
        OPENMS_LOG_DEBUG << "No peaks found in mobilogram. Returning peak at center of original mobilogram." << std::endl;
        return PeakPickerMobilogram::PeakPositions{0, im_size / 2, im_size-1};
      }

      // Find the iterator pointing to the maximum element
      auto max_it = std::max_element(intensities.begin(), intensities.end());

      // Get the index of the maximum element
      size_t max_index = std::distance(intensities.begin(), max_it);

      // Return the tuple
      return PeakPickerMobilogram::PeakPositions{left_widths[max_index], max_index, right_widths[max_index]};
    }

    void PeakPickerMobilogram::filterPeakIntensities_(Mobilogram& mobilogram,
                               size_t left_index,
                               size_t right_index) 
    {
      // Create a temporary vector to hold the filtered peaks
      std::vector<MobilityPeak1D> filtered_peaks;

      for (size_t i = left_index; i <= right_index; ++i) {
        const auto& peak = mobilogram[i];
        // Collect the peaks within the range
        filtered_peaks.push_back(peak); 
      }

      // Clear existing data and replace with filtered peaks
      mobilogram.clear();
      for (const auto& peak : filtered_peaks) {
        mobilogram.push_back(peak);
      }
    }

    void PeakPickerMobilogram::filterPeakIntensities_(std::vector<Mobilogram>& mobilograms,
                                size_t left_index,
                                size_t right_index) 
    {
      for (auto& mobilogram : mobilograms) {
        // Create a temporary vector to hold the filtered peaks
        std::vector<MobilityPeak1D> filtered_peaks;

        for (size_t i = left_index; i <= right_index; ++i) {
          const auto& peak = mobilogram[i];
          // Collect the peaks within the range
          filtered_peaks.push_back(peak); 
        }

        // Clear existing data and replace with filtered peaks
        mobilogram.clear(); 
        for (const auto& peak : filtered_peaks) {
          mobilogram.push_back(peak);
        }
      }
    }

    std::vector<double> PeakPickerMobilogram::extractFloatValues_(const OpenMS::DataArrays::FloatDataArray& floatDataArray)
    {
      std::vector<double> result;

      if (floatDataArray.empty()) {
        return result;
      }

      result.reserve(floatDataArray.size()); 

      for (size_t i = 0; i < floatDataArray.size(); ++i) {
        result.push_back(static_cast<double>(floatDataArray[i])); 
      }

      return result;
    }

    std::vector<std::size_t> PeakPickerMobilogram::extractIntValues_(const OpenMS::DataArrays::FloatDataArray& floatDataArray)
    {
      std::vector<std::size_t> result;

      if (floatDataArray.empty()) {
        return result;
      }

      result.reserve(floatDataArray.size());

      for (size_t i = 0; i < floatDataArray.size(); ++i) {
        result.push_back(static_cast<std::size_t>(floatDataArray[i])); // Convert and add to vector
      }

      return result;
    }

    PeakPickerMobilogram::PeakPositions PeakPickerMobilogram::filterTopPeak_(Mobilogram& picked_mobilogram, std::vector<Mobilogram>& mobilograms)
    {
      const auto& apex_abundance_data = extractFloatValues_(picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE]);
      const auto& leftwidth_data = extractIntValues_(picked_mobilogram.getFloatDataArrays()[IDX_OF_LEFTBORDER_IDX]);
      const auto& rightwidth_data = extractIntValues_(picked_mobilogram.getFloatDataArrays()[IDX_OF_RIGHTBORDER_IDX]);
      PeakPositions peak_pos = findHighestPeak_(apex_abundance_data, leftwidth_data, rightwidth_data, mobilograms[0].size());

      OPENMS_LOG_DEBUG << "  -- filtering mobilograms for highest peak at positions " << "(" << peak_pos.left << " - " << peak_pos.right << ")" << std::endl;

       filterPeakIntensities_(mobilograms, peak_pos.left, peak_pos.right);

       return peak_pos;
    }

    PeakPickerMobilogram::PeakPositions PeakPickerMobilogram::filterTopPeak_(Mobilogram& picked_mobilogram, Mobilogram& mobilogram)
    {
      const auto& apex_abundance_data = extractFloatValues_(picked_mobilogram.getFloatDataArrays()[IDX_ABUNDANCE]);
      const auto& leftwidth_data = extractIntValues_(picked_mobilogram.getFloatDataArrays()[IDX_OF_LEFTBORDER_IDX]);
      const auto& rightwidth_data = extractIntValues_(picked_mobilogram.getFloatDataArrays()[IDX_OF_RIGHTBORDER_IDX]);
      PeakPositions peak_pos = findHighestPeak_(apex_abundance_data, leftwidth_data, rightwidth_data, mobilogram.size());

      OPENMS_LOG_DEBUG << "  -- filtering mobilogram for highest peak at positions " << "(" << peak_pos.left << " - " << peak_pos.right << ")" << std::endl;

      filterPeakIntensities_(mobilogram, peak_pos.left, peak_pos.right);

      return peak_pos;
    }

    void PeakPickerMobilogram::removeOverlappingPeaks_(const Mobilogram& mobilogram, Mobilogram& picked_mobilogram)
    {
      if (picked_mobilogram.empty()) {return; }
      OPENMS_LOG_DEBUG << "Remove overlapping peaks now (size " << picked_mobilogram.size() << ")" << std::endl;
      Size current_peak = 0;
      // Find overlapping peaks
      for (Size i = 0; i < picked_mobilogram.size() - 1; i++)
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

          // Find the peak width and best IM
          double central_peak_mz = picked_mobilogram[i].getMobility();
          double next_peak_mz = picked_mobilogram[i + 1].getMobility();
          current_peak = findClosestPeak_(mobilogram, central_peak_mz, current_peak);
          Size next_peak = findClosestPeak_(mobilogram, next_peak_mz, current_peak);

          // adjust the right border of the current and left border of next
          Size k = 1;
          while ((current_peak + k) < mobilogram.size()
                 && (mobilogram[current_peak + k].getIntensity() < mobilogram[current_peak + k - 1].getIntensity()))
          {
            ++k;
          }
          Size new_right_border = current_peak + k - 1;
          k = 1;
          while ((next_peak - k + 1) > 0
                 && (mobilogram[next_peak - k].getIntensity() < mobilogram[next_peak - k + 1].getIntensity()))
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

          OPENMS_LOG_DEBUG << "New peak l: " << mobilogram[current_left_idx].getMobility() << " " << mobilogram[new_right_border].getMobility() << " int " << integrated_intensities_[i] << std::endl;
          OPENMS_LOG_DEBUG << "New peak r: " << mobilogram[new_left_border].getMobility() << " " << mobilogram[next_right_idx].getMobility() << " int " << integrated_intensities_[i + 1] << std::endl;


          right_width_[i] = new_right_border;
          left_width_[i + 1] = new_left_border;

        }
      }
    }

    void PeakPickerMobilogram::updateMembers_()
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
        write_sn_log_messages_ = (bool)param_.getValue("write_sn_log_messages").toBool();
        method_ = (String)param_.getValue("method").toString();

        Param sg_filter_parameters = sgolay_.getParameters();
        sg_filter_parameters.setValue("frame_length", sgolay_frame_length_);
        sg_filter_parameters.setValue("polynomial_order", sgolay_polynomial_order_);
        sgolay_.setParameters(sg_filter_parameters);

        Param gfilter_parameters = gauss_.getParameters();
        gfilter_parameters.setValue("gaussian_width", gauss_width_);
        gauss_.setParameters(gfilter_parameters);
    }
    }