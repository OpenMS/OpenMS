// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Author: Timo Sachsenberg, Mohammed Alhigaylan $
// $Maintainer: Timo Sachsenberg $
// -------------------------------------------------------------------------------------------------------------------------------------------

#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerIM.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/PROCESSING/SMOOTHING/GaussFilter.h>
#include <OpenMS/PROCESSING/SMOOTHING/SavitzkyGolayFilter.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/PROCESSING/RESAMPLING/LinearResamplerAlign.h>
#include <OpenMS/MATH/MISC/CubicSpline2d.h>
#include <OpenMS/MATH/MISC/SplineBisection.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FEATUREFINDER/MassTraceDetection.h>
#include <OpenMS/FEATUREFINDER/ElutionPeakDetection.h>
#include <iostream>
#include <deque>
#include <algorithm>
#include <limits>
#include <numeric>
#include <unordered_map>
#include <utility>


using namespace std;

// #define DEBUG_PICKER
#ifdef DEBUG_PICKER
#include <OpenMS/FORMAT/MzMLFile.h>
#endif 

namespace OpenMS
{
    double PeakPickerIM::computeOptimalSamplingRate(const vector<MSSpectrum>& spectra)
    {
      vector<double> mz_diffs;
      Size upper_peak_limit = 0;
      for (size_t s = 0; s < spectra.size(); ++s)
      {
        upper_peak_limit += spectra[s].size();
      }
      mz_diffs.reserve(upper_peak_limit);

      for (size_t s = 0; s < spectra.size(); ++s)
      {
        const MSSpectrum& spectrum = spectra[s];
        // The spectrum could have multiple ion mobility peaks at the same x position.
        // Sum the peak intensity
        MSSpectrum summed_trace;
        sumFrame_(spectrum, summed_trace, 0.0006, false);

        if (summed_trace.size() < 20)
        {
#ifdef DEBUG_PICKER
          std::cerr << "Skipping trace " << s << " because it has too few points ("
                    << summed_trace.size() << ")." << std::endl;
#endif
          continue; // skip this spectrum
        }

        for (size_t i = 1; i < summed_trace.size(); ++i)
        {
          double diff = summed_trace[i].getMZ() - summed_trace[i - 1].getMZ();
          mz_diffs.push_back(diff);
        }
        if (mz_diffs.size() > 1000) break; // 1000 diffs should be enough to estimate sampling
      }

      // If we found no valid m/z differences (traces too short)
      if (mz_diffs.empty())
      {
#ifdef DEBUG_PICKER
        std::cerr << "Warning: No valid m/z differences found in any spectra. Using default sampling rate of 0.01" << std::endl;
#endif
        return 0.01; // Fallback value
      }

      // Sort the differences to compute the 75th percentile threshold
      // This is needed in case there is a gap in the mobilogram. i+1 peak will skew the computed
      // sampling rate.
      std::sort(mz_diffs.begin(), mz_diffs.end());

      size_t percentile_index = static_cast<size_t>(mz_diffs.size() * 0.75);
      double threshold = mz_diffs[percentile_index];

#ifdef DEBUG_PICKER
      std::cerr << "75th percentile of position differences is: " << threshold << std::endl;
#endif

      // Filter out large differences (keep diffs <= threshold)
      vector<double> small_mz_diffs;
      for (double diff : mz_diffs)
      {
        if (diff <= threshold)
        {
          small_mz_diffs.push_back(diff);
        }
      }

      if (small_mz_diffs.empty())
      {
        std::cerr << "Warning: No valid small m/z differences found after filtering. Using default sampling rate of 0.01" << std::endl;
        return 0.01;
      }

      // Compute the mode
      std::unordered_map<double, int> freq_map;
      for (double diff : small_mz_diffs)
      {
        freq_map[diff]++;
      }

      double mode_sampling_rate = small_mz_diffs.front();
      int max_count = 0;

      for (const auto& [diff, count] : freq_map)
      {
        if (count > max_count)
        {
          mode_sampling_rate = diff;
          max_count = count;
        }
      }

#ifdef DEBUG_PICKER
      std::cout << "Computed optimal sampling rate: " << mode_sampling_rate << std::endl;
#endif

      return mode_sampling_rate;
    }

    // Function to compute the lower and upper m/z bounds based on ppm tolerance
    std::pair<double, double> PeakPickerIM::ppmBounds(double mz, double ppm)
    {
      ppm = ppm / 1e6;
      double delta_mz = (ppm * mz) / 2.0;

      double low = mz - delta_mz;
      double high = mz + delta_mz;

      return std::make_pair(low, high);
    }

    // PRECONDITION: input_spectrum is sorted by position (m/z peak or ion mobility)
    // This function sums peaks if they are nearly identical
    // OpenMS represents TimsTOF data in MSSpectrum() objects as one-array.
    // Example: There could be multiple 500.0 m/z peaks with different ion mobility values.
    // Example2: extracted mobilogram could have multiple 0.88 1/k values from different m/z peaks.
    // Peak picking (such as HiRes) will not work properly if there are multiple y measurements at a given x position.
    // Note: does not clear the output_spectrum but add peaks to it (required for fast padding)
    void PeakPickerIM::sumFrame_(const MSSpectrum& input_spectrum,
                                 MSSpectrum& output_spectrum,
                                 double tolerance,
                                 bool use_ppm)
    {
      if (input_spectrum.empty()) return;

      OPENMS_PRECONDITION(input_spectrum.isSorted(), "Spectrum must be sorted by m/z before summing peaks.");

      double current_mz = input_spectrum[0].getMZ();
      double current_intensity = input_spectrum[0].getIntensity();

      for (Size i = 1; i < input_spectrum.size(); ++i)
      {
        double next_mz = input_spectrum[i].getMZ();
        double next_intensity = input_spectrum[i].getIntensity();

        double delta_mz = std::abs(next_mz - current_mz);
        bool within_tolerance = use_ppm
                                  ? ((delta_mz / current_mz) * 1e6 <= tolerance)
                                  : (delta_mz <= tolerance);

        if (within_tolerance)
        {
          current_intensity += next_intensity;
        }
        else
        {
          output_spectrum.emplace_back(current_mz, current_intensity);
          current_mz = next_mz;
          current_intensity = next_intensity;
        }
      }
      output_spectrum.emplace_back(current_mz, current_intensity);
    }

    // We use peak FWHM (from PeakPickerHiRes) to extract ion mobility traces.
    // Given a picked m/z peak, we write a temporary MSSpectrum() object with ion mobility measurements
    // in place of m/z in Peak1D object. This facilitates peak picking in the ion mobility dimension.
    // To enable recomputing of m/z center after ion mobility peak picking, we tack raw m/z peak values
    // in FloatDataArrays().

    vector<MSSpectrum> PeakPickerIM::extractIonMobilityTraces(
      const MSSpectrum& picked_spectrum,
      const MSSpectrum& raw_spectrum)
    {
      const auto& float_data_arrays = picked_spectrum.getFloatDataArrays();

      // Find FWHM array in picked_spectrum
      const MSSpectrum::FloatDataArray* fwhm_array = nullptr;

      for (const auto& array : float_data_arrays)
      {
        if (array.getName() == "FWHM_ppm")
        {
          fwhm_array = &array;
          break;
        }
      }

      if (!fwhm_array)
      {
        std::cerr << "FWHM data array not found!" << std::endl;
        return {};
      }

      if (fwhm_array->size() != picked_spectrum.size())
      {
        std::cerr << "Size mismatch between FWHM array and picked peaks!" << std::endl;
        return {};
      }
      // Get the Ion Mobility array index from raw_spectrum
      const auto [im_data_index, im_unit] = raw_spectrum.getIMData();
      const auto& ion_mobility_array = raw_spectrum.getFloatDataArrays()[im_data_index];

      // Vector of MSSpectra for each picked m/z peak (each spectrum is a mobilogram trace)
      vector<MSSpectrum> mobility_traces;

      for (size_t i = 0; i < picked_spectrum.size(); ++i)
      {
        double picked_mz = picked_spectrum[i].getMZ();
        double fwhm_ppm = (*fwhm_array)[i];

        auto bounds = ppmBounds(picked_mz, fwhm_ppm);
        double lower_bound = bounds.first;
        double upper_bound = bounds.second;

        SignedSize center_idx = raw_spectrum.findNearest(picked_mz);

        if (center_idx == -1)
        {
          std::cerr << "No raw peaks found near picked m/z: " << picked_mz << std::endl;
          mobility_traces.emplace_back();
          continue;
        }

        MSSpectrum trace_spectrum; // A single mobilogram trace
        // Prepare FloatDataArray to store raw m/z values
        MSSpectrum::FloatDataArray raw_mz_array;
        raw_mz_array.setName("raw_mz");

        // Expand left
        SignedSize left_idx = center_idx;
        while (left_idx >= 0 && raw_spectrum[left_idx].getMZ() >= lower_bound)
        {
          trace_spectrum.emplace_back(ion_mobility_array[left_idx], raw_spectrum[left_idx].getIntensity()); // IM stored as m/z temporarily

          // Store the raw m/z
          raw_mz_array.push_back(raw_spectrum[left_idx].getMZ());

          --left_idx;
        }

        // Expand right
        SignedSize right_idx = center_idx + 1;
        while (right_idx < static_cast<SignedSize>(raw_spectrum.size()) &&
               raw_spectrum[right_idx].getMZ() <= upper_bound)
        {
          trace_spectrum.emplace_back(ion_mobility_array[right_idx], raw_spectrum[right_idx].getIntensity());

          // Store the raw m/z data in floatDataArrays()
          raw_mz_array.push_back(raw_spectrum[right_idx].getMZ());

          ++right_idx;
        }

        // Attach the raw m/z array to trace_spectrum
        auto& trace_float_arrays = trace_spectrum.getFloatDataArrays();
        trace_float_arrays.push_back(std::move(raw_mz_array));

        // Sort the trace_spectrum by ion mobility (m/z), while keeping raw m/z aligned
        trace_spectrum.sortByPosition(); // Note: having the float arrays attached ensures that sorting is performed on everything

        mobility_traces.push_back(std::move(trace_spectrum));
      }

      return mobility_traces;
    }

    // Function to compute m/z centers from mobilogram_traces and picked_traces
    MSSpectrum PeakPickerIM::computeCentroids_(const vector<MSSpectrum>& mobilogram_traces,
                              const vector<MSSpectrum>& picked_traces)
    {
      MSSpectrum centroided_frame;

      // Create float data arrays to house ion mobility data and peaks FWHM
      MSSpectrum::FloatDataArray ion_mobility_array;
      ion_mobility_array.setName("Ion Mobility");

      MSSpectrum::FloatDataArray ion_mobility_fwhm;
      ion_mobility_fwhm.setName("IM FWHM");

      MSSpectrum::FloatDataArray mz_fwhm_array;
      mz_fwhm_array.setName("MZ FWHM");

#ifdef DEBUG_PICKER
      std::cout << "picked_traces.size(): " << picked_traces.size() << std::endl;
#endif
      // Loop over picked traces and their corresponding raw mobilogram traces
      for (size_t i = 0; i < picked_traces.size(); ++i)
      {
        // std::cout << "Looping through picked_trace that has .. " << picked_traces[i].size() << std::endl;
        const MSSpectrum& picked_trace = picked_traces[i];
        const MSSpectrum& raw_trace = mobilogram_traces[i];

        const auto& picked_float_arrays = picked_trace.getFloatDataArrays();

        if (picked_float_arrays.empty())
        {
          std::cerr << "No IM FWHM array found for picked_trace " << i << "!" << std::endl;
          continue;
        }

        // Assuming the first FloatDataArray holds the ion mobility peak FWHM values
        const auto& fwhm_array = picked_float_arrays[0];

        if (fwhm_array.size() != picked_trace.size())
        {
          std::cerr << "FWHM array size mismatch with picked_trace size!" << std::endl;
          continue;
        }

        // Get the FloatDataArrays from raw_trace (assumed to hold the raw m/z values)
        const auto& raw_float_arrays = raw_trace.getFloatDataArrays();

        if (raw_float_arrays.empty())
        {
          std::cerr << "No raw m/z peaks found for raw_trace " << i << "!" << std::endl;
          continue;
        }

        // Assume the first array holds the raw m/z values
        const auto& raw_mz_values = raw_float_arrays[0];

        if (raw_mz_values.size() != raw_trace.size())
        {
          std::cerr << "raw_mz_values size mismatch with raw_trace size!" << std::endl;
          continue;
        }

#ifdef DEBUG_PICKER
        std::cout << "\n--- Processing picked_trace " << i << " ---\n";
#endif

        // Create reusable objects outside the loop to reduce memory allocations
        MSSpectrum raw_peaks_within_bounds;
        MSSpectrum raw_mz_peaks;
        vector<double> mz_values;
        vector<double> intensity_values;
        vector<size_t> indices;
        vector<double> sorted_mz;
        vector<double> sorted_intensity;
        
        // Iterate through picked peaks in this trace
        for (Size j = 0; j < picked_trace.size(); ++j)
        {
          double centroid_im = picked_trace[j].getMZ();   // Ion mobility centroid (stored as m/z)
          double fwhm = fwhm_array[j];

          double im_lower = centroid_im - (fwhm / 2.0);
          double im_upper = centroid_im + (fwhm / 2.0);

#ifdef DEBUG_PICKER
          std::cout << "Picked peak " << j << " IM centroid: " << centroid_im
                    << " ion mobility FWHM: " << fwhm
                    << " --> IM bounds: [" << im_lower << ", " << im_upper << "]" << std::endl;
#endif
          // Use findNearest() to get the index of the closest peak in the raw mobilogram trace
          SignedSize center_idx = raw_trace.findNearest(centroid_im);

          if (center_idx == -1)
          {
            std::cerr << "Could not find nearest peak to centroid_im in raw_trace!" << std::endl;
            continue;
          }

          // Clear the spectrum for reuse
          raw_peaks_within_bounds.clear(true);

          // --- Expand Left ---
          SignedSize left_idx = center_idx;
          while (left_idx >= 0 && raw_trace[left_idx].getMZ() >= im_lower)
          {
            Peak1D new_peak;
            new_peak.setMZ(raw_mz_values[left_idx]);                      // m/z from FloatDataArray
            new_peak.setIntensity(raw_trace[left_idx].getIntensity());    // intensity from raw_trace
            raw_peaks_within_bounds.push_back(new_peak);

            --left_idx;
          }

          // --- Expand Right ---
          SignedSize right_idx = center_idx + 1;
          while (right_idx < static_cast<SignedSize>(raw_trace.size()) &&
                 raw_trace[right_idx].getMZ() <= im_upper)
          {
            Peak1D new_peak;
            new_peak.setMZ(raw_mz_values[right_idx]);
            new_peak.setIntensity(raw_trace[right_idx].getIntensity());
            raw_peaks_within_bounds.push_back(new_peak);

            ++right_idx;
          }

#ifdef DEBUG_PICKER
          std::cout << "Picked IM peak " << j << ": collected " << raw_peaks_within_bounds.size()
                    << " raw m/z points between IM [" << im_lower << ", " << im_upper << "]" << std::endl;
#endif

          // If we only retrieved one raw peak, pass it over to centroided_frame as is
          // Resampling and smoothing the raw data distorts the intensity values.
          // We recompute the m/z peak maxima and intensity using spline
          if (raw_peaks_within_bounds.size() == 1)
          {
            const Peak1D& single_peak = raw_peaks_within_bounds.front();

            // Add it directly to centroided_frame
            centroided_frame.push_back(single_peak);

            // Push corresponding ion mobility and FWHM arrays
            ion_mobility_array.push_back(centroid_im);
            ion_mobility_fwhm.push_back(fwhm);
            mz_fwhm_array.push_back(0.0);

#ifdef DEBUG_PICKER
            std::cout << "[INFO] Only one raw peak found. Added directly to centroided_frame. m/z: "
                      << single_peak.getMZ() << " intensity: " << single_peak.getIntensity() << std::endl;
#endif
            // Skip the rest of the loop and move on to the next picked_trace peak
            continue;
          }
        
          // If not sorted, sort both the peaks array and the corresponding raw m/z values array
          if (!raw_peaks_within_bounds.isSorted())
          {
            raw_peaks_within_bounds.sortByPosition();
          }
           
          // Clear the spectrum for reuse
          raw_mz_peaks.clear(true);
          sumFrame_(raw_peaks_within_bounds, raw_mz_peaks, 0.1, true);
          if (raw_mz_peaks.empty())
          {
            OPENMS_LOG_DEBUG << "No data in raw_mz_peaks for picked IM peak " << j << "!" << std::endl;
            continue;
          }

          // Summing peaks could result in spectrum.size() == 1 which causes error.
          // in this case, simply sum the intensity values
          if (raw_mz_peaks.size() == 1)
          {
            centroided_frame.push_back(raw_mz_peaks[0]);
            ion_mobility_array.push_back(centroid_im);
            ion_mobility_fwhm.push_back(fwhm);
            mz_fwhm_array.push_back(0.0);

#ifdef DEBUG_PICKER
            const Peak1D& single_peak = raw_mz_peaks[0];
            std::cout << "[INFO] sumFrame_ reduced peaks to a single entry. Added directly to centroided_frame. m/z: " << single_peak.getMZ()
                      << " intensity: " << single_peak.getIntensity() << std::endl;
#endif
            continue;
          }

          // Clear and reuse vectors for spline data
          mz_values.clear();
          intensity_values.clear();
          
          // Reserve space for efficiency
          mz_values.reserve(raw_mz_peaks.size());
          intensity_values.reserve(raw_mz_peaks.size());

          // Initialize sorting flag
          bool is_sorted = true;

          // Populate the vectors and check if sorted at the same time
          for (size_t i = 0; i < raw_mz_peaks.size(); ++i)
          {
            double current_mz = raw_mz_peaks[i].getMZ();
            double current_intensity = raw_mz_peaks[i].getIntensity();
            
            // Check if still sorted (compare with previous value if not the first element)
            if (i > 0 && current_mz < mz_values.back())
            {
              is_sorted = false;
            }
            
            mz_values.push_back(current_mz);
            intensity_values.push_back(current_intensity);
          }

          // Sort the vectors if needed TODO: check with Mohammed if this is needed
          if (!is_sorted)
          {
            // Reuse indices vector
            indices.resize(mz_values.size());
            for (size_t i = 0; i < indices.size(); ++i)
            {
              indices[i] = i;
            }
            
            // Sort indices based on m/z values
            std::sort(indices.begin(), indices.end(),
                    [&mz_values](size_t i1, size_t i2) {
                      return mz_values[i1] < mz_values[i2];
                    });
            
            // Reuse sorted vectors
            sorted_mz.resize(mz_values.size());
            sorted_intensity.resize(intensity_values.size());
            
            // Reorder both vectors using the sorted indices
            for (size_t i = 0; i < indices.size(); ++i)
            {
              sorted_mz[i] = mz_values[indices[i]];
              sorted_intensity[i] = intensity_values[indices[i]];
            }
            
            // Replace the original vectors with the sorted ones
            mz_values = std::move(sorted_mz);
            intensity_values = std::move(sorted_intensity);
          }

          // Initialize spline with the two vectors
          CubicSpline2d spline(mz_values, intensity_values);

          // Define boundaries
          const double left_bound = mz_values.front();
          const double right_bound = mz_values.back();

          // Find maximum via spline bisection
          double apex_mz = (left_bound + right_bound) / 2.0;
          double apex_intensity = 0.0;

          const double max_search_threshold = 1e-6;

          Math::spline_bisection(spline, left_bound, right_bound, apex_mz, apex_intensity, max_search_threshold);

#ifdef DEBUG_PICKER
          std::cout << "Apex m/z: " << apex_mz << std::endl;
          std::cout << "Apex intensity: " << apex_intensity << std::endl;
#endif

          // FWHM calculation (same binary search as before)
          double half_height = apex_intensity / 2.0;
          const double fwhm_search_threshold = 0.01 * half_height;

          // ---- Left side search ----
          double mz_left = left_bound;
          double mz_center = apex_mz;
          double int_mid = 0.0;
          double mz_mid = mz_left;

          if (spline.eval(mz_left) > half_height)
          {
            mz_mid = mz_left;
          }
          else
          {
            do
            {
              mz_mid = (mz_left + mz_center) / 2.0;
              int_mid = spline.eval(mz_mid);

              if (int_mid < half_height)
              {
                mz_left = mz_mid;
              }
              else
              {
                mz_center = mz_mid;
              }
            } while (std::fabs(int_mid - half_height) > fwhm_search_threshold);
          }
          double fwhm_left_mz = mz_mid;

          // ---- Right side search ----
          double mz_right = right_bound;
          mz_center = apex_mz;

          if (spline.eval(mz_right) > half_height)
          {
            mz_mid = mz_right;
          }
          else
          {
            do
            {
              mz_mid = (mz_right + mz_center) / 2.0;
              int_mid = spline.eval(mz_mid);

              if (int_mid < half_height)
              {
                mz_right = mz_mid;
              }
              else
              {
                mz_center = mz_mid;
              }

            } while (std::fabs(int_mid - half_height) > fwhm_search_threshold);
          }
          double fwhm_right_mz = mz_mid;

          // ---- FWHM result ----
          double mz_fwhm = fwhm_right_mz - fwhm_left_mz;

#ifdef DEBUG_PICKER
          std::cout << "Left m/z at half height: " << fwhm_left_mz << std::endl;
          std::cout << "Right m/z at half height: " << fwhm_right_mz << std::endl;
          std::cout << "m/z FWHM: " << mz_fwhm << std::endl;
#endif

          centroided_frame.emplace_back(apex_mz, apex_intensity);
          ion_mobility_array.push_back(centroid_im);
          ion_mobility_fwhm.push_back(fwhm);
          mz_fwhm_array.push_back(mz_fwhm);
        }

#ifdef DEBUG_PICKER
        std::cout << "--- Finished processing picked_trace " << i << " ---\n\n";
#endif
      }

      auto& centroided_frame_fda = centroided_frame.getFloatDataArrays();
      centroided_frame_fda.push_back(std::move(ion_mobility_array));
      centroided_frame_fda.push_back(std::move(ion_mobility_fwhm));
      centroided_frame_fda.push_back(std::move(mz_fwhm_array));
      centroided_frame.sortByPosition();

#ifdef DEBUG_PICKER
      std::cout << "Peaks in centroided frame: " << centroided_frame.size() << std::endl;
      std::cout << "Printing centroided_frame inside ComputerCenters function " << std::endl;
      for (const auto& peak : centroided_frame)
      {
        std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
      }
#endif
      return centroided_frame;
    }

    void removeAllFloatDataArraysExcept(OpenMS::MSSpectrum& spectrum, const String& keep_name)
    {
      auto& float_arrays = spectrum.getFloatDataArrays();
  
      // Use remove_if to move all elements to remove to the end
      auto new_end = std::remove_if(float_arrays.begin(), float_arrays.end(),
                               [&keep_name](const MSSpectrum::FloatDataArray& array) {
                                 return array.getName() != keep_name;  // Remove if NOT the one we want to keep
                               });
  
      // Erase the removed elements
      float_arrays.erase(new_end, float_arrays.end());
    }

    /*
    Param PeakPickerIM::getDefaultParameters() const
    {
      Param p;
      // internal PickIMTraces parameters that users should not change.

      p.setValue("signal_to_noise", 0.0, "Signal to noise threshold for peak picking");
      p.setValue("spacing_difference_gap", 0.0, "The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.");
      p.setValue("spacing_difference", 0.0, "Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.");
      p.setValue("missing", 0, "Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.");
      p.setValue("report_FWHM", "true");
      p.setValue("report_FWHM_unit", "absolute");

      return p;
    }
     */
    PeakPickerIM::PeakPickerIM()
        : DefaultParamHandler("PeakPickerIM")
    {
      // --- PickIMTraces parameters ---
      defaults_.setValue("pickIMTraces:sum_tolerance_mz",        1.0,   "Tolerance for summing adjacent m/z peaks (ppm)");
      defaults_.setValue("pickIMTraces:sum_tolerance_im",        0.0006,"Tolerance for summing adjacent ion mobility peaks (1/k0)");
      defaults_.setValue("pickIMTraces:gauss_ppm_tolerance",     5.0,   "Gaussian smoothing m/z tolerance in ppm");
      defaults_.setValue("pickIMTraces:sgolay_frame_length",     5,     "Savitzky-Golay smoothing frame length");
      defaults_.setValue("pickIMTraces:sgolay_polynomial_order", 3,     "Savitzky-Golay smoothing polynomial order");
      // --- PickIMCluster parameters ---
      defaults_.setValue("pickIMCluster:ppm_tolerance_cluster", 50.0, "m/z tolerance in ppm for clustering");
      defaults_.setValue("pickIMCluster:im_tolerance_cluster", 0.1, "Ion mobility tolerance in 1/k for clustering");
      // --- PickIMElutionProfiles parameters ---
      defaults_.setValue("pickIMElutionProfiles:ppm_tolerance_elution", 50.0, "Mass trace m/z tolerance in ppm");

      defaultsToParam_();  // copies defaults_ into param_
      updateMembers_();    // caches into member variables
    }
    void PeakPickerIM::updateMembers_()
    {
      sum_tolerance_mz_      = (double)param_.getValue("pickIMTraces:sum_tolerance_mz");
      sum_tolerance_im_      = (double)param_.getValue("pickIMTraces:sum_tolerance_im");
      gauss_ppm_tolerance_   = (double)param_.getValue("pickIMTraces:gauss_ppm_tolerance");
      sgolay_frame_length_   = (int)param_.getValue("pickIMTraces:sgolay_frame_length");
      sgolay_polynomial_order_= (int)param_.getValue("pickIMTraces:sgolay_polynomial_order");

      ppm_tolerance_cluster_ = (double)param_.getValue("pickIMCluster:ppm_tolerance_cluster");
      im_tolerance_cluster_ = (double)param_.getValue("pickIMCluster:im_tolerance_cluster");

      ppm_tolerance_elution_ = (double)param_.getValue("pickIMElutionProfiles:ppm_tolerance_elution");

    }



    void PeakPickerIM::pickIMTraces(MSSpectrum& spectrum)
    {
      /*
      // IM format determination (Temporarily commented out)
      IMFormat format = IMTypes::determineIMFormat(spectrum);
      switch (format)
      {
          case IMFormat::NONE:
              return; // no IM data
          case IMFormat::CENTROIDED:
              return; // already centroided
          case IMFormat::UNKNOWN:
              throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                  "IMFormat set to UNKNOWN after determineIMFormat. This should never happen.",
                  String(NamesOfIMFormat[(size_t)format]));
      }
      */

      // ************************************************* PART I *****************************************************
      // ------------------------------------------ mass-to-charge peak picking -------------------------------------
      // ------------------------------------------ Step a1: Sum m/z peaks ------------------------------------------
      // project all timsTOF peaks into the m/z axis using sumFrame_
      // The ppm tolerance is a dynamic way of testing m/z floats being almost identical. The raw intensity is summed.
      MSSpectrum summed_spectrum;
      sumFrame_(spectrum, summed_spectrum, sum_tolerance_mz_, true);
#ifdef DEBUG_PICKER
      std::cout << "Spectrum after sumFrame_ has " << summed_spectrum.size() << " peaks." << std::endl;
#endif
      // ------------------------------------------ step 2a: smooth ------------------------------------------
      // Apply gaussian smoothing to the peaks projected into the m/z axis. This facilitates peak picking
      // in the m/z dimension and subseqent mobilogram extraction for each picked m/z peak.
#ifdef DEBUG_PICKER
      std::cout << "Applying Gaussian smoothing..." << std::endl;
#endif
      GaussFilter gauss_filter;
      double gauss_tol = gauss_ppm_tolerance_;
      Param gauss_params;
      gauss_params.setValue("ppm_tolerance", gauss_tol);
      gauss_params.setValue("use_ppm_tolerance", "true");
      gauss_filter.setParameters(gauss_params);
      gauss_filter.filter(summed_spectrum);
#ifdef DEBUG_PICKER
      std::cout << "Spectrum after Gaussian smoothing has " << summed_spectrum.size() << " peaks." << std::endl;
      for (const auto& peak : summed_spectrum)
      {
        std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
      }
#endif

      // ------------------------------------------ step 3a: m/z Peak Picking ------------------------------------------
      // Pick peaks in the m/z axis and toggle reporting peak width at half max (FWHM)
      // we will use the FWHM of each picked m/z peak to extract mobilograms.
      PeakPickerHiRes picker_mz;
      Param picker_mz_p;
      picker_mz_p.setValue("signal_to_noise", 0.0);
      picker_mz_p.setValue("report_FWHM", "true");
      picker_mz_p.setValue("report_FWHM_unit", "relative");
      picker_mz.setParameters(picker_mz_p);
      MSSpectrum picked_spectrum;
      picker_mz.pick(summed_spectrum, picked_spectrum);
#ifdef DEBUG_PICKER
      std::cout << "Size of picked spectrum: " << picked_spectrum.size() << std::endl;
#endif

      // ------------------------------------------ step 4a: Extraction mobilograms ------------------------------------------
      // Using m/z peaks FWHM, we iteratively extract ion mobility traces (mobilograms) from the raw spectrum.
      // To rescue weak signal in the extracted mobilograms, we use linear resampling.
      // For linear resampling, it is recommended to use a sampling rate equal or higher than the raw sampling rate.
      // We dynamically determine the raw sampling rate from well-populated extracted mobilograms
      // (currently we have this hard-coded as +20 raw peaks in a mobilogram to be considered well-populated).

      auto mobilogram_traces = extractIonMobilityTraces(picked_spectrum, spectrum);

      // TODO: Should we add a multipler to the linear sampling rate? We may boost algorithm speed and/or result if we try >1 multiplier
      double sampling_rate = computeOptimalSamplingRate(mobilogram_traces) * 1;
      Param resampler_param;
      resampler_param.setValue("spacing", sampling_rate);
      resampler_param.setValue("ppm", "false");

#ifdef DEBUG_PICKER
      std::cout << "Using sampling rate... : " << sampling_rate << std::endl;
#endif

#ifdef DEBUG_PICKER
      for (size_t i = 0; i < mobilogram_traces.size(); ++i)
      {
        std::cout << "Trace " << i << " contains " << mobilogram_traces[i].size() << " points in ion mobility space." << std::endl;
      }
#endif
      // ************************************************* PART II *****************************************************
      // ------------------------------------------ Ion mobility peak picking ------------------------------------------
      // ------------------------------------------ part 1b: sum ion mobility peaks ------------------------------------
      // An extract ion mobilogram can have two peaks with identicial 1/k value and cuase issues in the peak picking steps.
      // Example: if raw sampling rate is 0.0012 1/k -- then ion mobility peak 0.8800 1/k and 0.8806 1/k should be combined.
      // Use 0.0006 1/k as default. This parameter may need to change depending on ion mobility ramp tamp
      // (it is currently optimized for 100 ms ramp time)


      // prepare picked ion mobility objects (we are internally using MSSpectrum object for downstram peak picking inputs).
      vector<MSSpectrum> picked_traces;
      // TODO: check why there are so many mobilogram_traces that are empty. leads to segfault later
      mobilogram_traces.erase(
        std::remove_if(mobilogram_traces.begin(), mobilogram_traces.end(),
                   [](const auto& trace) { return trace.empty(); }),
        mobilogram_traces.end());

      for (size_t i = 0; i < mobilogram_traces.size(); ++i)
      {
        MSSpectrum& trace = mobilogram_traces[i];

#ifdef DEBUG_PICKER
        std::cout << "\n--- Processing Trace " << i << " ---\n";
        std::cout << "Original trace has " << trace.size() << " peaks." << std::endl;
#endif
        MSSpectrum summed_trace;
        summed_trace.reserve(trace.size() + 1);
        summed_trace.emplace_back(-1.0, -1.0);
        sumFrame_(trace, summed_trace, sum_tolerance_im_, false);
#ifdef DEBUG_PICKER
        std::cout << "Trace after sumFrame_ has " << summed_trace.size() << " peaks." << std::endl;
#endif
        // ------------------------------------------ part 2b: smooth and resample --------------------------------
        // Prepare mobilograms for SGolay smoothing.
        // To avoid edge effect, we will pad the edges with (sgolay_frame_length_ -1 / 2.0) points.
        double im_start = summed_trace[1].getMZ();
        double im_end = summed_trace.back().getMZ();

#ifdef DEBUG_PICKER
        std::cout << "Original summed trace ion mobility range: [" << im_start << ", " << im_end << "]" << std::endl;
#endif
        int padding_points = static_cast<int>(std::ceil((sgolay_frame_length_ - 1) / 2.0));

        Peak1D front_padding;
        front_padding.setMZ(im_start - padding_points * sampling_rate);
        front_padding.setIntensity(0.0);
        summed_trace[0] = front_padding;

        Peak1D back_padding;
        back_padding.setMZ(im_end + padding_points * sampling_rate);
        back_padding.setIntensity(0.0);
        summed_trace.push_back(back_padding);

#ifdef DEBUG_PICKER
        std::cout << "Padded summed trace im range: [" << summed_trace.front().getMZ() << ", " << summed_trace.back().getMZ() << "]" << std::endl;
#endif

        // linear resample to rescue weak signal
        LinearResamplerAlign lin_resampler;
        lin_resampler.setParameters(resampler_param);
        lin_resampler.raster(summed_trace);
#ifdef DEBUG_PICKER
        std::cout << "Size of resampled trace: " << summed_trace.size() << " peaks." << std::endl;
        for (const auto& peak : summed_trace)
        {
          std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
        }
#endif
        // SGolay smooth prior to peak picking
        SavitzkyGolayFilter sgolay_filter;
        Param sgolay_params;
        sgolay_params.setValue("frame_length", sgolay_frame_length_);
        sgolay_params.setValue("polynomial_order", sgolay_polynomial_order_);
        sgolay_filter.setParameters(sgolay_params);
        sgolay_filter.filter(summed_trace);

#ifdef DEBUG_PICKER
        std::cout << "Trace after Savitzky-Golay smoothing has " << summed_trace.size() << " peaks." << std::endl;
        for (const auto& peak : summed_trace)
        {
          std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
        }
#endif

        // ------------------------------------------ part 3b: im peak picking --------------------------------
        // apply PeakPickerHiRes to pick ion mobility peaks.
        // PeakPickerHiRes can be applied to chromatograms. We reasoned the same set of parameters ideal for
        // chromatograms is also applicable for mobilograms.
        // Each raw mobilogram contains a float data array with raw m/z values.
        // We will use the ion mobility peak FWHM to define min/max ion mobility boundary
        // and recompute the m/z center based on the ion mobility peak.
        PeakPickerHiRes picker_im;
        Param picker_im_p;
        picker_im_p.setValue("signal_to_noise", 0.0);
        picker_im_p.setValue("spacing_difference_gap", 0.0);
        picker_im_p.setValue("spacing_difference", 0.0);
        picker_im_p.setValue("missing", 0);
        picker_im_p.setValue("report_FWHM", "true");
        picker_im_p.setValue("report_FWHM_unit", "absolute");
        picker_im.setParameters(picker_im_p);

        MSSpectrum picked_trace;
        picker_im.pick(summed_trace, picked_trace);
        picked_traces.push_back(std::move(picked_trace));

#ifdef DEBUG_PICKER
        std::cout << "--- Finished Processing Trace " << i << " ---\n\n";
#endif
      }

      // Recompute m/z centers and output centroided frame
      MSSpectrum centroided_frame = computeCentroids_(mobilogram_traces, picked_traces);

#ifdef DEBUG_PICKER
      std::cout << "--- Centroided frame has  " << centroided_frame.size() << " --- peaks.\n";
#endif
      // Copy only SpectrumSettings from the input into the centroided result
      static_cast<SpectrumSettings&>(centroided_frame) = static_cast<const SpectrumSettings&>(spectrum);
      centroided_frame.setMSLevel(spectrum.getMSLevel());
      centroided_frame.setName(spectrum.getName());
      centroided_frame.setRT(spectrum.getRT());
      removeAllFloatDataArraysExcept(centroided_frame, "Ion Mobility");
      centroided_frame.setIMFormat(IMFormat::CENTROIDED);
      spectrum = std::move(centroided_frame);
      

#ifdef DEBUG_PICKER
      // Print peaks for debugging
      std::cout << "--- Spectrum final output object has ..  " << spectrum.size() << " --- peaks.\n";
      for (const auto& peak : spectrum)
      {
        std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
      }
#endif
    }

  // static function to pick IM clusters
/*
  void PeakPickerIM::pickIMCluster(MSSpectrum& input, double ppm_tolerance, double im_tolerance)
  {
        if (input.empty()) return;

        // Get IM data array
        const auto [im_data_index, im_unit] = input.getIMData();
        auto& im_data = input.getFloatDataArrays()[im_data_index];

        
         // @brief Internal structure to hold peak information during processing        
        struct Point {
            double mz;        ///< m/z value
            double im;        ///< ion mobility value
            float intensity; ///< intensity value
            Point(double mz_val, double im_val, double int_val) :
                mz(mz_val), im(im_val), intensity(int_val) {}
        };

        // Convert peaks to Points for processing
        vector<Point> points;
        points.reserve(input.size());
        
        for (size_t i = 0; i < input.size(); ++i)
        {
            const auto& peak = input[i];
            points.emplace_back(peak.getMZ(), im_data[i], peak.getIntensity());
        }

        // Sort points by m/z if not already sorted
        if (!input.isSorted())
        {
            sort(points.begin(), points.end(), 
                 [](const Point& a, const Point& b) { return a.mz < b.mz; });
        }

        // Deques for maintaining min/max of IM values in current window
        deque<double> minDeque;  ///< Maintains minimum IM values
        deque<double> maxDeque;  ///< Maintains maximum IM values

        double ppm_factor = 1.0 + (ppm_tolerance * 1e-6);
        size_t start = 0;
        vector<Point> averaged_points;

         // @brief Helper function to compute intensity-weighted averages for a cluster
         // @param s Start index of the cluster
         // @param e End index of the cluster
        auto finalize_cluster = [&](int s, int e) {
            if (s <= e) {
                double sum_intensity = 0.0;
                double sum_mv = 0.0;
                double sum_iv = 0.0;
                
                for (int idx = s; idx <= e; ++idx) {
                    sum_intensity += points[idx].intensity;
                    sum_mv += points[idx].mz * points[idx].intensity;
                    sum_iv += points[idx].im * points[idx].intensity;
                }
                
                averaged_points.emplace_back(
                    sum_mv / sum_intensity,
                    sum_iv / sum_intensity,
                    sum_intensity
                );
            }
        };

        // Main clustering loop
        for (size_t end = 0; end < points.size(); ++end)
        {
            double current_i = points[end].im;

            // Update minDeque
            while (!minDeque.empty() && minDeque.back() > current_i)
                minDeque.pop_back();
            minDeque.push_back(current_i);

            // Update maxDeque
            while (!maxDeque.empty() && maxDeque.back() < current_i)
                maxDeque.pop_back();
            maxDeque.push_back(current_i);

            // Check and update window start
            while (start <= end)
            {
                double m_min = points[start].mz;
                double m_max = points[end].mz;
                bool m_ok = (m_max <= m_min * ppm_factor);

                double i_min = minDeque.front();
                double i_max = maxDeque.front();
                bool i_ok = ((i_max - i_min) <= im_tolerance);

                if (m_ok && i_ok)
                    break;

                // Remove start point from deques
                double i_to_remove = points[start].im;
                if (i_to_remove == minDeque.front())
                    minDeque.pop_front();
                if (i_to_remove == maxDeque.front())
                    maxDeque.pop_front();
                start++;
            }

            // Check if next point would violate constraints
            if (end + 1 < points.size())
            {
                double next_m = points[end + 1].mz;
                double next_i = points[end + 1].im;

                double m_min = points[start].mz;
                bool m_ok = (next_m <= m_min * ppm_factor);

                double i_min_next = min(minDeque.front(), next_i);
                double i_max_next = max(maxDeque.front(), next_i);
                bool i_ok = ((i_max_next - i_min_next) <= im_tolerance);

                if (!(m_ok && i_ok))
                {
                    finalize_cluster(start, end);
                    start = end + 1;
                    minDeque.clear();
                    maxDeque.clear();
                }
            }
        }

        // Finalize last cluster
        if (start < points.size())
        {
            finalize_cluster(start, points.size() - 1);
        }

        // Convert averaged points back to MSSpectrum
        input.resize(averaged_points.size());
        input.shrink_to_fit();        
        
        im_data.resize(averaged_points.size());
        im_data.shrink_to_fit();

        for (size_t i = 0; i != averaged_points.size(); ++i)
        {
          const auto& p = averaged_points[i];
          input[i].setMZ(p.mz);
          input[i].setIntensity(p.intensity);
          im_data[i] = p.im;
        }

        input.sortByPosition();
        input.updateRanges(); // TODO: maybe needed
  }
*/
    void PeakPickerIM::pickIMCluster(OpenMS::MSSpectrum& spectrum, double ppm_tolerance_cluster_, double im_tolerance_cluster_)
    {
      if (spectrum.empty()) return;

      // Get IM data array
      if (!spectrum.containsIMData())
      {
        OPENMS_LOG_WARN << "No ion mobility data found in spectrum." << std::endl;
        return;
      }
      const auto [im_data_index, im_unit] = spectrum.getIMData();
      auto& im_data = spectrum.getFloatDataArrays()[im_data_index];


      struct Point {
        double mz;
        double im;
        float intensity;
        OpenMS::Size original_index;

        Point(double mz_val, double im_val, float int_val, OpenMS::Size idx) :
            mz(mz_val), im(im_val), intensity(int_val), original_index(idx) {}
      };

      // Convert peaks to Points (same as before)
      std::vector<Point> points;
      points.reserve(spectrum.size());
      for (OpenMS::Size i = 0; i < spectrum.size(); ++i) {
        const auto& peak = spectrum[i];
        points.emplace_back(peak.getMZ(), im_data[i], peak.getIntensity(), i);
      }

      // --- Setup for Clustering ---

      // 1. Create m/z-sorted indices (needed for cluster expansion)
      std::vector<OpenMS::Size> mz_sorted_indices(points.size());
      std::iota(mz_sorted_indices.begin(), mz_sorted_indices.end(), 0);
      std::sort(mz_sorted_indices.begin(), mz_sorted_indices.end(),
                [&](OpenMS::Size a, OpenMS::Size b) {
                  return points[a].mz < points[b].mz;
                });

      // Create reverse lookup: original_index -> mz_sorted_position (needed for cluster expansion)
      std::vector<OpenMS::Size> original_to_sorted_pos(points.size());
      for (OpenMS::Size i = 0; i < points.size(); ++i) {
        original_to_sorted_pos[mz_sorted_indices[i]] = i;
      }

      // 2. ***OPTIMIZATION: Create intensity-sorted indices for seed picking***
      std::vector<OpenMS::Size> intensity_sorted_indices(points.size());
      std::iota(intensity_sorted_indices.begin(), intensity_sorted_indices.end(), 0);
      std::sort(intensity_sorted_indices.begin(), intensity_sorted_indices.end(),
                [&](OpenMS::Size a, OpenMS::Size b) {
                  // Sort DESCENDING by intensity
                  return points[a].intensity > points[b].intensity;
                });

      // 3. Keep track of used points
      std::vector<bool> used(points.size(), false);
      OpenMS::Size num_used = 0;

      // Store results temporarily
      std::vector<Point> averaged_points;
      averaged_points.reserve(points.size());

      double ppm_factor = ppm_tolerance_cluster_ * 1e-6;

      // --- Main Clustering Loop ---
      // Iterate through peaks in descending order of intensity to find seeds
      for (OpenMS::Size current_intensity_rank = 0; current_intensity_rank < points.size(); ++current_intensity_rank)
      {
        // 4. Get the original index of the next potential seed (highest intensity first)
        OpenMS::Size seed_original_idx = intensity_sorted_indices[current_intensity_rank];

        // 5. Check if this peak has already been used (part of a previous cluster)
        if (used[seed_original_idx]) {
          continue; // Skip to the next highest intensity peak
        }

        // --- Found an unused seed, start clustering ---

        // 6. Initialize the cluster with the seed
        std::vector<OpenMS::Size> current_cluster_indices;
        current_cluster_indices.push_back(seed_original_idx);
        used[seed_original_idx] = true;
        num_used++;

        double cluster_mz_min = points[seed_original_idx].mz;
        double cluster_mz_max = points[seed_original_idx].mz;
        double cluster_im_min = points[seed_original_idx].im;
        double cluster_im_max = points[seed_original_idx].im;

        // 7. Expand the cluster using m/z sorted neighbors (same logic as before)
        OpenMS::Size seed_sorted_pos = original_to_sorted_pos[seed_original_idx];
        OpenMS::SignedSize left_idx = static_cast<OpenMS::SignedSize>(seed_sorted_pos) - 1;
        OpenMS::SignedSize right_idx = static_cast<OpenMS::SignedSize>(seed_sorted_pos) + 1;

        bool changed = true;
        while (changed) {
          changed = false;

          // Check left neighbor
          while (left_idx >= 0) {
            OpenMS::Size candidate_original_idx = mz_sorted_indices[left_idx];
            if (!used[candidate_original_idx]) {
              const auto& candidate_point = points[candidate_original_idx];
              double potential_mz_min = std::min(cluster_mz_min, candidate_point.mz);
              double potential_mz_max = std::max(cluster_mz_max, candidate_point.mz); // Fixed: Use max for the max value
              double potential_im_min = std::min(cluster_im_min, candidate_point.im);
              double potential_im_max = std::max(cluster_im_max, candidate_point.im);

              // Fixed m/z tolerance calculation
              bool mz_ok = (potential_mz_max - potential_mz_min) <= (potential_mz_min * ppm_factor);
              bool im_ok = (potential_im_max - potential_im_min) <= im_tolerance_cluster_;

              if (mz_ok && im_ok) {
                current_cluster_indices.push_back(candidate_original_idx);
                used[candidate_original_idx] = true;
                num_used++;
                cluster_mz_min = potential_mz_min;
                cluster_mz_max = potential_mz_max; // Fixed: Update max value
                cluster_im_min = potential_im_min;
                cluster_im_max = potential_im_max;
                left_idx--;
                changed = true;
                break; // Added point, break inner while to re-evaluate from outer changed loop
              } else {
                left_idx = -1; // Stop checking left this round
                break;
              }
            } else {
              left_idx--; // Skip used point
            }
          }

          // Check right neighbor
          while (right_idx < static_cast<OpenMS::SignedSize>(points.size())) {
            OpenMS::Size candidate_original_idx = mz_sorted_indices[right_idx];
            if (!used[candidate_original_idx]) {
              const auto& candidate_point = points[candidate_original_idx];
              double potential_mz_min = std::min(cluster_mz_min, candidate_point.mz); // Fixed: Use min
              double potential_mz_max = std::max(cluster_mz_max, candidate_point.mz);
              double potential_im_min = std::min(cluster_im_min, candidate_point.im);
              double potential_im_max = std::max(cluster_im_max, candidate_point.im);

              // Fixed m/z tolerance calculation
              bool mz_ok = (potential_mz_max - potential_mz_min) <= (potential_mz_min * ppm_factor);
              bool im_ok = (potential_im_max - potential_im_min) <= im_tolerance_cluster_;

              if (mz_ok && im_ok) {
                current_cluster_indices.push_back(candidate_original_idx);
                used[candidate_original_idx] = true;
                num_used++;
                cluster_mz_min = potential_mz_min; // Fixed: Update min value
                cluster_mz_max = potential_mz_max;
                cluster_im_min = potential_im_min;
                cluster_im_max = potential_im_max;
                right_idx++;
                changed = true;
                break; // Added point, break inner while to re-evaluate from outer changed loop
              } else {
                right_idx = points.size(); // Stop checking right this round
                break;
              }
            } else {
              right_idx++; // Skip used point
            }
          }
        } // End cluster expansion (while changed)

        // 8. Finalize the current cluster (same logic as before)
        if (!current_cluster_indices.empty()) {
          double sum_intensity = 0.0;
          double sum_mz_intensity = 0.0;
          double sum_im_intensity = 0.0;

          for (OpenMS::Size original_idx : current_cluster_indices) {
            const auto& p = points[original_idx];
            sum_intensity += p.intensity;
            sum_mz_intensity += p.mz * p.intensity;
            sum_im_intensity += p.im * p.intensity;
          }

          if (sum_intensity > std::numeric_limits<double>::epsilon()) {
            averaged_points.emplace_back(
              sum_mz_intensity / sum_intensity,
              sum_im_intensity / sum_intensity,
              static_cast<float>(sum_intensity),
              0 // Original index is meaningless for averaged points
            );
          }
        }

        // Optimization: Check if all points are processed
        if (num_used == points.size()) {
          break; // Exit the main loop early
        }

      } // End main loop (for intensity_sorted_indices)

      // 9. Update spectrum (same logic as before)
      spectrum.resize(averaged_points.size());
      spectrum.shrink_to_fit();
      im_data.resize(averaged_points.size());
      im_data.shrink_to_fit();

      for (size_t i = 0; i != averaged_points.size(); ++i)
      {
        const auto& p = averaged_points[i];
        spectrum[i].setMZ(p.mz);
        spectrum[i].setIntensity(p.intensity);
        im_data[i] = p.im;
      }

      spectrum.sortByPosition();
      spectrum.updateRanges();
    } // End of pickIMCluster function
    void PeakPickerIM::pickIMElutionProfiles(MSSpectrum& input, double ppm_tolerance_elution_)
    {
      if (input.empty()) return;
      // Get IM data array
      if (!input.containsIMData())
      {
        OPENMS_LOG_WARN << "No ion mobility data found in spectrum." << std::endl;
        return;
      }
      const auto [im_data_index, im_unit] = input.getIMData();
      auto& im_data = input.getFloatDataArrays()[im_data_index];
      // convert to MSExperiment and set drift time as RT
      MSExperiment frame_as_spectra = IMDataConverter::reshapeIMFrameToMany(input);
      for (auto& s : frame_as_spectra)
      {
        s.setRT(s.getDriftTime());
        s.setDriftTime(-1);
        s.setMSLevel(1);
      }

#ifdef DEBUG_IM_PICKER
  // write out IM frame as RT/MZ for debugging purposes to test algorithm that yet don't support the IM dimension
  MzMLFile().store("debug" + String(input.getRT()) + ".mzML", frame_as_spectra);
#endif

      if (frame_as_spectra.size() <= 3 ) return;

      // detect mass traces in IM frame
      MassTraceDetection mte;
      Param param = mte.getParameters();
      // disable most filter criteria
      param.setValue("min_trace_length", -1.0);
      param.setValue("max_trace_length", -1.0);
      param.setValue("noise_threshold_int", 0.1); // only ignore 0 peaks
      param.setValue("chrom_peak_snr", 0.0);
      param.setValue("reestimate_mt_sd", "false");
      param.setValue("mass_error_ppm", ppm_tolerance_elution_);
      param.setValue("trace_termination_criterion", "outlier");
      param.setValue("trace_termination_outlier", 1);

      mte.setLogType(ProgressLogger::NONE);
      mte.setParameters(param);
      vector<MassTrace> output_mt;
      mte.run(frame_as_spectra, output_mt);

      ElutionPeakDetection epd;
      param = epd.getParameters();
      param.setValue("chrom_fwhm", 0.01);
      param.setValue("chrom_peak_snr", 0.0);
      param.setValue("width_filtering", "off");
      param.setValue("min_fwhm", -1.0);
      param.setValue("max_fwhm", 1e6);
      param.setValue("masstrace_snr_filtering", "false");
      epd.setParameters(param);

      std::vector<MassTrace> split_mtraces;
      epd.detectPeaks(output_mt, split_mtraces);
      output_mt.clear();

      // copy mass traces centroids back to peaks
      input.resize(split_mtraces.size());
      input.shrink_to_fit();

      im_data.resize(split_mtraces.size());
      im_data.shrink_to_fit();

      for (Size i = 0; i < split_mtraces.size(); ++i)
      {
        const MassTrace& mt = split_mtraces[i];
        input[i].setMZ(mt.getCentroidMZ());
        input[i].setIntensity(mt.computeIntensitySum()); //TODO: somewhat related: why is getIntensity() zero?
        im_data[i] = mt.getCentroidRT(); // IM
      }
      input.sortByPosition(); // TODO: maybe needed
      input.updateRanges(); // TODO: maybe needed
    }

} // namespace OpenMS
