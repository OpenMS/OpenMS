// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/NOISEESTIMATION/SignalToNoiseEstimatorMedian.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>

// #define DEBUG_PEAK_PICKING

namespace OpenMS
{

  /**
   * @brief A small structure to hold peak candidates
   *
  */
  struct PeakCandidate
  {
    int index;
    double peak_apex_intensity;

    double integrated_intensity;
    double leftWidth;
    double rightWidth;
    float mz;
  };

  bool sort_peaks_by_intensity(const PeakCandidate& a, const PeakCandidate& b); // prototype
  bool sort_peaks_by_intensity(const PeakCandidate& a, const PeakCandidate& b)
  {
    return a.peak_apex_intensity > b.peak_apex_intensity;
  }

  /**
  @brief This class implements a peak-picking algorithm for high-resolution MS
  data (specifically designed for TOF-MS data).

  This peak-picking algorithm detects ion signals in profile data and
  reconstructs the corresponding peak shape by identifying the left and right
  borders of the peak. It reports the area under the peak as intensity and the
  weighted m/z values as the m/z value as well as left/right border.
  Furthermore, it next tries to improve the peak positioning iteratively using
  the m/z center computed in the last iteration. This allows for refinement in
  the peak boundaries and more accurate determination of peak center and borders.

  Its approach is similar to the PeakPickerHiRes but additionally uses an
  iterative approach to find and re-center peaks.

  - First, it uses the PeakPickerHiRes to find seeds or candidate peaks.
  - Next it uses n iterations to re-center those peaks and compute left/right
    borders for each peak.
  - Finally it removes peaks that are within the borders of other peaks.

  So far, this peak picker was mainly tested on high resolution TOF-MS data.

  @htmlinclude OpenMS_PeakPickerIterative.parameters

  @note The peaks must be sorted according to ascending m/z!

  @ingroup PeakPicking

  */
  class OPENMS_DLLAPI PeakPickerIterative :
    public DefaultParamHandler,
    public ProgressLogger
  {

private:
    double signal_to_noise_;
    double peak_width_;
    double spacing_difference_;
    int sn_bin_count_;
    int nr_iterations_;
    double sn_win_len_;
    bool check_width_internally_;

public:

    /// Constructor
    PeakPickerIterative() :
      DefaultParamHandler("PeakPickerIterative"),
      ProgressLogger()
    {
      defaults_.setValue("signal_to_noise_", 1.0, "Signal to noise value, each peak is required to be above this value (turn off by setting it to 0.0)");
      defaults_.setValue("peak_width", 0.0, "Expected peak width half width in Dalton - peaks will be extended until this half width is reached (even if the intensitity is increasing). In conjunction with check_width_internally it will also be used to remove peaks whose spacing is larger than this value.");


      defaults_.setValue("spacing_difference", 1.5, "Difference between peaks in multiples of the minimal difference to continue. The higher this value is set, the further apart peaks are allowed to be to still extend a peak. E.g. if the value is set to 1.5 and in a current peak the minimal spacing between peaks is 10 mDa, then only peaks at most 15 mDa apart will be added to the peak.", {"advanced"});
      defaults_.setValue("sn_bin_count_", 30, "Bin count for the Signal to Noise estimation.", {"advanced"});
      defaults_.setValue("nr_iterations_", 5, "Nr of iterations to perform (how many times the peaks are re-centered).", {"advanced"});
      defaults_.setMinInt("nr_iterations_", 1);
      defaults_.setValue("sn_win_len_", 20.0, "Window length for the Signal to Noise estimation.", {"advanced"});

      defaults_.setValue("check_width_internally", "false", "Delete peaks where the spacing is larger than the peak width (should be set to true to avoid artefacts)", {"advanced"});
      defaults_.setValidStrings("check_width_internally", {"true","false"});

      defaults_.setValue("ms1_only", "false", "Only do MS1");
      defaults_.setValidStrings("ms1_only", {"true","false"});
      defaults_.setValue("clear_meta_data", "false", "Delete meta data about peak width");
      defaults_.setValidStrings("clear_meta_data", {"true","false"});

      // write defaults into Param object param_
      defaultsToParam_();
    }

    void updateMembers_() override
    {
      signal_to_noise_ = (double)param_.getValue("signal_to_noise_");
      peak_width_ = (double)param_.getValue("peak_width");
      spacing_difference_ = (double)param_.getValue("spacing_difference");
      sn_bin_count_ = (double)param_.getValue("sn_bin_count_");
      nr_iterations_ = (double)param_.getValue("nr_iterations_");
      sn_win_len_ = (double)param_.getValue("sn_win_len_");

      check_width_internally_ = param_.getValue("check_width_internally").toBool();
    }

    /// Destructor
    ~PeakPickerIterative() override {}

private:

    /*
     * This will re-center the peaks by using the seeds (ordered by intensity) to
     * find raw signals that may belong to this peak. Then the peak is centered
     * using a weighted average.
     * Signals are added to the peak as long as they are still inside the
     * peak_width or as long as the signal intensity keeps falling. Also the
     * distance to the previous signal and the whether the signal is below the
     * noise level is taken into account.
     * This function implements a single iteration of this algorithm.
     *
    */
    void pickRecenterPeaks_(const MSSpectrum& input,
                              std::vector<PeakCandidate>& PeakCandidates,
                              SignalToNoiseEstimatorMedian<MSSpectrum>& snt) const
    {
      for (Size peak_it = 0; peak_it < PeakCandidates.size(); peak_it++)
      {
        int i = PeakCandidates[peak_it].index;
        double central_peak_mz = input[i].getMZ(), central_peak_int = input[i].getIntensity();
        double left_neighbor_mz = input[i - 1].getMZ(), left_neighbor_int = input[i - 1].getIntensity();
        double right_neighbor_mz = input[i + 1].getMZ(), right_neighbor_int = input[i + 1].getIntensity();

        // MZ spacing sanity checks
        double left_to_central = std::fabs(central_peak_mz - left_neighbor_mz);
        double central_to_right = std::fabs(right_neighbor_mz - central_peak_mz);
        double min_spacing = (left_to_central < central_to_right) ? left_to_central : central_to_right;
        double est_peak_width = peak_width_;

        if (check_width_internally_ && (left_to_central > est_peak_width || central_to_right > est_peak_width))
        {
          // something has gone wrong, the points are further away than the peak width -> delete this peak
          PeakCandidates[peak_it].integrated_intensity = -1;
          PeakCandidates[peak_it].leftWidth = -1;
          PeakCandidates[peak_it].rightWidth = -1;
          PeakCandidates[peak_it].mz = -1;
          continue;
        }

        std::map<double, double> peak_raw_data;

        peak_raw_data[central_peak_mz] = central_peak_int;
        peak_raw_data[left_neighbor_mz] = left_neighbor_int;
        peak_raw_data[right_neighbor_mz] = right_neighbor_int;

        // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //
        // COPY - from PeakPickerHiRes
        // peak core found, now extend it to the left
        Size k = 2;
        while ((i - k + 1) > 0
              && std::fabs(input[i - k].getMZ() - peak_raw_data.begin()->first) < spacing_difference_ * min_spacing
              && (input[i - k].getIntensity() < peak_raw_data.begin()->second
                 || std::fabs(input[i - k].getMZ() - central_peak_mz) < est_peak_width)
               )
        {
          if (signal_to_noise_ > 0.0)
          {
            if (snt.getSignalToNoise(i - k) < signal_to_noise_)
            {
              break;
            }
          }
          peak_raw_data[input[i - k].getMZ()] = input[i - k].getIntensity();
          ++k;
        }
        double leftborder = input[i - k + 1].getMZ();

        // to the right
        k = 2;
        while ((i + k) < input.size()
              && std::fabs(input[i + k].getMZ() - peak_raw_data.rbegin()->first) < spacing_difference_ * min_spacing
              && (input[i + k].getIntensity() < peak_raw_data.rbegin()->second
                 || std::fabs(input[i + k].getMZ() - central_peak_mz) < est_peak_width)
               )
        {
          if (signal_to_noise_ > 0.0)
          {
            if (snt.getSignalToNoise(i + k) < signal_to_noise_)
            {
              break;
            }
          }

          peak_raw_data[input[i + k].getMZ()] = input[i + k].getIntensity();
          ++k;
        }
        // END COPY - from PeakPickerHiRes
        // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // //

        double rightborder = input[i + k - 1].getMZ();

        double weighted_mz = 0;
        double integrated_intensity = 0;
        for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it)
        {
          weighted_mz += map_it->first * map_it->second;
          integrated_intensity += map_it->second;
        }
        weighted_mz /= integrated_intensity;

        // store the data
        PeakCandidates[peak_it].integrated_intensity = integrated_intensity;
        PeakCandidates[peak_it].leftWidth = leftborder;
        PeakCandidates[peak_it].rightWidth = rightborder;
        PeakCandidates[peak_it].mz = weighted_mz;

        // find the closest raw signal peak to where we just put our peak and store it
        double min_diff = std::fabs(weighted_mz - input[i].getMZ());
        int min_i = i;

        // Search to the left
        for (int m = 1; i - m > 0 && leftborder < input[i - m].getMZ(); m++)
        {
          if (std::fabs(weighted_mz - input[i - m].getMZ()) < min_diff)
          {
            min_diff = std::fabs(weighted_mz - input[i - m].getMZ());
            min_i = i - m;
          }
        }
        // Search to the right
        for (int m = 1; i - m > 0 && rightborder > input[i + m].getMZ(); m++)
        {
          if (std::fabs(weighted_mz - input[i + m].getMZ()) < min_diff)
          {
            min_diff = std::fabs(weighted_mz - input[i + m].getMZ());
            min_i = i + m;
          }
        }
        PeakCandidates[peak_it].index = min_i;
      }
    }

public:

    /*
      This will pick one single spectrum. The PeakPickerHiRes is used to
      generate seeds, these seeds are then used to re-center the mass and
      compute peak width and integrated intensity of the peak.
     
      Finally, other peaks that would fall within the primary peak are
      discarded
     
      The output are the remaining peaks.
    */
    void pick(const MSSpectrum& input, MSSpectrum& output)
    {
      // don't pick a spectrum with less than 3 data points
      if (input.size() < 3) return;

      // copy the spectrum meta data
      copySpectrumMeta(input, output);
      output.setType(SpectrumSettings::CENTROID);
      output.getFloatDataArrays().clear();

      std::vector<PeakCandidate> PeakCandidates;
      MSSpectrum picked_spectrum;

      // Use the PeakPickerHiRes to find candidates ...
      OpenMS::PeakPickerHiRes pp;
      Param pepi_param = OpenMS::PeakPickerHiRes().getDefaults();
      pepi_param.setValue("signal_to_noise", signal_to_noise_);
      pepi_param.setValue("spacing_difference", spacing_difference_);
      pp.setParameters(pepi_param);
      pp.pick(input, picked_spectrum);

      // after picking peaks, we store the closest index of the raw spectrum and the picked intensity
      std::vector<PeakCandidate> newPeakCandidates_;
      Size j = 0;
      OPENMS_LOG_DEBUG << "Candidates " << picked_spectrum.size() << std::endl;
      for (Size k = 0; k < input.size() && j < picked_spectrum.size(); k++)
      {
        if (input[k].getMZ() > picked_spectrum[j].getMZ())
        {
          OPENMS_LOG_DEBUG << "got a value " << k << " @ " << input[k] << std::endl;
          PeakCandidate pc = { /*.index=*/ static_cast<int>(k), /*.intensity=*/ picked_spectrum[j].getIntensity(), -1, -1, -1, -1};
          newPeakCandidates_.push_back(pc);
          j++;
        }
      }

      PeakCandidates = newPeakCandidates_;
      std::sort(PeakCandidates.begin(), PeakCandidates.end(), sort_peaks_by_intensity);

      // signal-to-noise estimation
      SignalToNoiseEstimatorMedian<MSSpectrum > snt;
      if (signal_to_noise_ > 0.0)
      {
        Param snt_parameters = snt.getParameters();
        snt_parameters.setValue("win_len", sn_win_len_);
        snt_parameters.setValue("bin_count", sn_bin_count_);
        snt.setParameters(snt_parameters);
        snt.init(input);
      }

      // The peak candidates are re-centered and the width is computed for each peak
      for (int i = 0; i < nr_iterations_; i++)
      {
        pickRecenterPeaks_(input, PeakCandidates, snt);
      }

      output.getFloatDataArrays().resize(3);
      output.getFloatDataArrays()[0].setName("IntegratedIntensity");
      output.getFloatDataArrays()[1].setName("leftWidth");
      output.getFloatDataArrays()[2].setName("rightWidth");

      // Go through all candidates and exclude all lower-intensity candidates
      // that are within the borders of another peak
      OPENMS_LOG_DEBUG << "Will now merge candidates" << std::endl;
      for (Size peak_it = 0; peak_it < PeakCandidates.size(); peak_it++)
      {
        if (PeakCandidates[peak_it].leftWidth < 0) continue;

        //Remove all peak candidates that are enclosed by this peak
        for (Size m = peak_it + 1; m < PeakCandidates.size(); m++)
        {
          if (PeakCandidates[m].mz >= PeakCandidates[peak_it].leftWidth && PeakCandidates[m].mz <= PeakCandidates[peak_it].rightWidth)
          {
            OPENMS_LOG_DEBUG << "Remove peak " << m <<  " : " << PeakCandidates[m].mz << " "  <<
              PeakCandidates[m].peak_apex_intensity << " (too close to " << PeakCandidates[peak_it].mz <<
              " " << PeakCandidates[peak_it].peak_apex_intensity <<  ")" << std::endl;
            PeakCandidates[m].leftWidth = PeakCandidates[m].rightWidth = -1;
          }
        }

        Peak1D peak;
        peak.setMZ(PeakCandidates[peak_it].mz);
        peak.setIntensity(PeakCandidates[peak_it].integrated_intensity);
        output.push_back(peak);

        OPENMS_LOG_DEBUG << "Push peak " << peak_it << "  " << peak << std::endl;
        output.getFloatDataArrays()[0].push_back(PeakCandidates[peak_it].integrated_intensity);
        output.getFloatDataArrays()[1].push_back(PeakCandidates[peak_it].leftWidth);
        output.getFloatDataArrays()[2].push_back(PeakCandidates[peak_it].rightWidth);
      }

      OPENMS_LOG_DEBUG << "Found seeds: " << PeakCandidates.size() << " / Found peaks: " << output.size() << std::endl;
      output.sortByPosition();
    }

    void pickExperiment(const PeakMap& input, PeakMap& output)
    {
      // make sure that output is clear
      output.clear(true);

      // copy experimental settings
      static_cast<ExperimentalSettings&>(output) = input;

      // resize output with respect to input
      output.resize(input.size());

      bool ms1_only = param_.getValue("ms1_only").toBool();
      bool clear_meta_data = param_.getValue("clear_meta_data").toBool();

      Size progress = 0;
      startProgress(0, input.size(), "picking peaks");
      for (Size scan_idx = 0; scan_idx != input.size(); ++scan_idx)
      {
        if (ms1_only && (input[scan_idx].getMSLevel() != 1))
        {
          output[scan_idx] = input[scan_idx];
        }
        else
        {
          pick(input[scan_idx], output[scan_idx]);
          if (clear_meta_data) {output[scan_idx].getFloatDataArrays().clear();}
        }
        setProgress(progress++);
      }
      endProgress();
    }

  };

}

