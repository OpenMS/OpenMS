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
#include <OpenMS/CONCEPT/LogStream.h>

#include <iostream>

namespace OpenMS
{
    double PeakPickerIM::computeOptimalSamplingRate(const std::vector<MSSpectrum>& spectra)
    {
      std::vector<double> mz_diffs;
      Size upper_peak_limit = 0;
      for (size_t s = 0; s < spectra.size(); ++s)
      {
        upper_peak_limit += spectra[s].size();
      }
      mz_diffs.reserve(upper_peak_limit);

      for (size_t s = 0; s < spectra.size(); ++s)
      {
        const MSSpectrum& spectrum = spectra[s];
        // The spectrum could have multiple peaks at the same x position.
        // Sum the peak intensity
        MSSpectrum summed_trace;
        SumFrame(spectrum, summed_trace, 7000.0);

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
      std::vector<double> small_mz_diffs;
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

    // PRECONDITION: input_spectrum is sorted by m/z
    // This function sums peaks if they are nearly identical
    // OpenMS represents TimsTOF data in MSSpectrum() objects as one-array.
    // There could be multiple 500.0 m/z peaks with different ion mobility values.
    // Peak picking (such as HiRes) will not work properly if there are multiple y measurements at a given x m/z position.
    // Note: does not clear the output_spectrum but add peaks to it (required for fast padding)
    void PeakPickerIM::SumFrame(const MSSpectrum& input_spectrum, MSSpectrum& output_spectrum, double ppm_tolerance)
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
        double ppm_diff = (delta_mz / current_mz) * 1e6;

        if (ppm_diff <= ppm_tolerance)
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

    std::vector<MSSpectrum> PeakPickerIM::extractIonMobilityTraces(
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

      // Get the Ion Mobility array from raw_spectrum
      const auto& raw_float_data_arrays = raw_spectrum.getFloatDataArrays();
      const MSSpectrum::FloatDataArray* ion_mobility_array = nullptr;

      for (const auto& array : raw_float_data_arrays)
      {
        if (array.getName() == "Ion Mobility")
        {
          ion_mobility_array = &array;
          break;
        }
      }

      if (!ion_mobility_array)
      {
        std::cerr << "Ion Mobility data array not found in raw_spectrum!" << std::endl;
        return {};
      }

      // Vector of MSSpectra for each picked m/z peak (each spectrum is a mobilogram trace)
      std::vector<MSSpectrum> mobility_traces;

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
          trace_spectrum.emplace_back((*ion_mobility_array)[left_idx], raw_spectrum[left_idx].getIntensity()); // currently Ion Mobility as m/z

          // Store the raw m/z
          raw_mz_array.push_back(raw_spectrum[left_idx].getMZ());

          --left_idx;
        }

        // Expand right
        SignedSize right_idx = center_idx + 1;
        while (right_idx < static_cast<SignedSize>(raw_spectrum.size()) &&
               raw_spectrum[right_idx].getMZ() <= upper_bound)
        {
          trace_spectrum.emplace_back((*ion_mobility_array)[right_idx], raw_spectrum[right_idx].getIntensity());

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
    MSSpectrum PeakPickerIM::ComputeCenters(const std::vector<MSSpectrum>& mobilogram_traces,
                              const std::vector<MSSpectrum>& picked_traces)
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
        std::vector<double> mz_values;
        std::vector<double> intensity_values;
        std::vector<size_t> indices;
        std::vector<double> sorted_mz;
        std::vector<double> sorted_intensity;
        
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
          SumFrame(raw_peaks_within_bounds, raw_mz_peaks, 0.1);
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
            std::cout << "[INFO] SumFrame reduced peaks to a single entry. Added directly to centroided_frame. m/z: " << single_peak.getMZ()
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


    PeakPickerIM::PeakPickerIM() :
        parameters_(getDefaultParameters())
    {
    }

    // Destructor
    PeakPickerIM::~PeakPickerIM()
    {
    }

    Param PeakPickerIM::getDefaultParameters() const
    {
      Param p;
      p.setValue("signal_to_noise", 0.0, "Signal to noise threshold for peak picking");
      p.setValue("spacing_difference_gap", 0.0, "The extension of a peak is stopped if the spacing between two subsequent data points exceeds 'spacing_difference_gap * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. '0' to disable the constraint. Not applicable to chromatograms.");
      p.setValue("spacing_difference", 0.0, "Maximum allowed difference between points during peak extension, in multiples of the minimal difference between the peak apex and its two neighboring points. If this difference is exceeded a missing point is assumed (see parameter 'missing'). A higher value implies a less stringent peak definition, since individual signals within the peak are allowed to be further apart. '0' to disable the constraint. Not applicable to chromatograms.");
      p.setValue("missing", 0, "Maximum number of missing points allowed when extending a peak to the left or to the right. A missing data point occurs if the spacing between two subsequent data points exceeds 'spacing_difference * min_spacing'. 'min_spacing' is the smaller of the two spacings from the peak apex to its two neighboring points. Not applicable to chromatograms.");
      p.setValue("report_FWHM", "true");
      p.setValue("report_FWHM_unit", "absolute");
      return p;
    }
    void PeakPickerIM::updateMembers_()
    {
      // This function is a placeholder for potential future use.
    }

    void PeakPickerIM::setParameters(const Param& param)
    {
      parameters_ = param;
      updateMembers_();
    }

    Param PeakPickerIM::getParameters() const
    {
      return parameters_;
    }

    void PeakPickerIM::pickIMTraces(MSSpectrum& spectrum)
    {
      // Debugging: Print the input spectrum size
      std::cout << "Processing spectrum with " << spectrum.size() << " peaks." << std::endl;

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


      // --- Step 1a: Sum m/z peaks
      // First, we project all timsTOF peaks into the m/z axis using SumFrame
      // The ppm tolerance is a dynamic way of testing m/z floats being almost identical. Set it to 0.1 ppm
      MSSpectrum summed_spectrum;
      SumFrame(spectrum, summed_spectrum, 0.1);
#ifdef DEBUG_PICKER
      std::cout << "Spectrum after SumFrame has " << summed_spectrum.size() << " peaks." << std::endl;
#endif
      // --- step 2a: smooth the data projected to the m/z axis.
#ifdef DEBUG_PICKER
      std::cout << "Applying Gaussian smoothing..." << std::endl;
#endif
      GaussFilter gauss_filter;

      // Set Gaussian filter parameters. 5 ppm m/z is good approximation (make this a user parameter!)
      Param gauss_params;
      gauss_params.setValue("ppm_tolerance", 5.0);
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

      // ---step 3a: Apply PeakPickerHiRes and make sure the peak FWHM is reported in ppm.
      // (Maybe this should change to absolute FWHM value...?).
      PeakPickerHiRes picker;
      Param hirez_mz_p;
      hirez_mz_p.setValue("signal_to_noise", 0.0);
      hirez_mz_p.setValue("report_FWHM", "true");
      hirez_mz_p.setValue("report_FWHM_unit", "relative");

      picker.setParameters(hirez_mz_p);
      MSSpectrum picked_spectrum;

      picker.pick(summed_spectrum, picked_spectrum);
#ifdef DEBUG_PICKER
      std::cout << "Size of picked spectrum: " << picked_spectrum.size() << std::endl;
#endif

      // ---step 4a: Extract ion mobility traces for each picked m/z peak
      auto mobilogram_traces = extractIonMobilityTraces(picked_spectrum, spectrum);

      // --- compute optimal sampling rate from well-populated mobilograms in this frame.
      // This is currently set to +20 peaks in a mobilogram. (Should this be a user parameter?)

      // Add a parameter to allow user to control sampling). Here we simply multiply by 4.
      double sampling_rate = computeOptimalSamplingRate(mobilogram_traces) * 1;
      Param resampler_param;
      resampler_param.setValue("spacing", sampling_rate, "Spacing of the resampled output peaks.");
      resampler_param.setValue("ppm", "false", "Whether spacing is in ppm or Th");

#ifdef DEBUG_PICKER
      std::cout << "Using sampling rate... : " << sampling_rate << std::endl;
#endif


      // *************** PART II ***************

      // for each ion mobility trace, we process the raw signal, peak pick
      // and recompute m/z and ion mobility centroid.
#ifdef DEBUG_PICKER
      for (size_t i = 0; i < mobilogram_traces.size(); ++i)
      {
        std::cout << "Trace " << i << " contains " << mobilogram_traces[i].size() << " points in ion mobility space." << std::endl;
      }
#endif
      std::vector<MSSpectrum> picked_traces;


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
        // --- Step 1b: Sum peaks that are too close ---
        MSSpectrum summed_trace;
        summed_trace.reserve(trace.size() + 1);
        summed_trace.emplace_back(-1.0, -1.0); // add now as we need an entry for padding later (faster than insert at front)
        SumFrame(trace, summed_trace, 7000.0); // 7000 ppm tolerance (temporary. Change function to accept absolute number)
#ifdef DEBUG_PICKER
        std::cout << "Trace after SumFrame has " << summed_trace.size() << " peaks." << std::endl;
#endif
        // Determine im boundaries of current mobilogram. Add 10 padding points (should this be a parameter?)
        // If you do not pad the edges, peaks on the edge will have an odd shape and not be picked by PeakPickerHiRes!
        double im_start = summed_trace[1].getMZ(); // first entry after padding
        double im_end = summed_trace.back().getMZ();

#ifdef DEBUG_PICKER
        std::cout << "Original summed trace ion mobility range: [" << im_start << ", " << im_end << "]" << std::endl;
#endif
        // for SGolay smoothing -- pad the edges with the same number as SGolay window size.
        // If the padded region is not matching, it will borrow the left points from the right
        // aka window 15 --> 7 points to the right and 7 points to the left pushed to the right.
        // This can create a fake peak in the padded edges.
        Peak1D front_padding;
        front_padding.setMZ(im_start - 15 * sampling_rate);
        front_padding.setIntensity(0.0);
        summed_trace[0] = front_padding;

        Peak1D back_padding;
        back_padding.setMZ(im_end + 15 * sampling_rate);
        back_padding.setIntensity(0.0);
        summed_trace.push_back(back_padding);

#ifdef DEBUG_PICKER
        std::cout << "Padded summed trace im range: [" << summed_trace.front().getMZ() << ", " << summed_trace.back().getMZ() << "]" << std::endl;
#endif

        // --- Step 2b: Resample the trace ---
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
        /*
        // --- Step 3b: Apply Gaussian Smoothing ---
        std::cout << "Applying Gaussian smoothing..." << std::endl;

        GaussFilter gauss_filter_2;
        Param gauss_params;
        gauss_params.setValue("gaussian_width", 0.005);
        gauss_params.setValue("use_ppm_tolerance", "false"); // or "true" for ppm-based width

        gauss_filter_2.setParameters(gauss_params);
        gauss_filter_2.filter(summed_trace);

#ifdef DEBUG_PICKER
        std::cout << "Trace after Gaussian smoothing has " << summed_trace.size() << " peaks." << std::endl;
        for (const auto& peak : summed_trace)
        {
          std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
        }
#endif
        */


        // --- Step 3b: Apply SGolay Smoothing ---
        SavitzkyGolayFilter sgolay_filter;
        Param sgolay_params;

        sgolay_params.setValue("frame_length", 15);
        sgolay_params.setValue("polynomial_order", 3);
        sgolay_filter.setParameters(sgolay_params);

        sgolay_filter.filter(summed_trace);

#ifdef DEBUG_PICKER
        std::cout << "Trace after Savitzky-Golay smoothing has " << summed_trace.size() << " peaks." << std::endl;
        for (const auto& peak : summed_trace)
        {
          std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
        }
#endif

        // --- Step 4b: Apply PeakPickerHiRes ---
        // We will use ion mobility peak FWHM to define min/max ion mobility boundaries.
        // Revisit the raw traces and compute intensity weighted ion mobility centroids.
        // The ion mobility traces also contains raw m/z peaks in FloatDataArrays.
        // This makes it convenient  to re-compute m/z centroid.
        PeakPickerHiRes picker;
        picker.setParameters(parameters_);

        MSSpectrum picked_trace;
        picker.pick(summed_trace, picked_trace);
        // populated picked_traces vector
        picked_traces.push_back(picked_trace);

#ifdef DEBUG_PICKER
        std::cout << "--- Finished Processing Trace " << i << " ---\n\n";
#endif
      }

      // Recompute m/z centers and output centroided frame
      MSSpectrum centroided_frame = ComputeCenters(mobilogram_traces, picked_traces);

#ifdef DEBUG_PICKER
      std::cout << "--- Centroided frame has  " << centroided_frame.size() << " --- peaks.\n";
#endif

      // Replace the input spectrum with the centroided result
      spectrum = centroided_frame;

#ifdef DEBUG_PICKER
      // Print peaks for debugging
      std::cout << "--- Spectrum final output object has ..  " << spectrum.size() << " --- peaks.\n";
      for (const auto& peak : spectrum)
      {
        std::cout << "m/z: " << peak.getMZ() << ", intensity: " << peak.getIntensity() << std::endl;
      }
#endif
    }
} // namespace OpenMS
