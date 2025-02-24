// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/INTERFACES/DataStructures.h>
#include <OpenMS/INTERFACES/ISpectrumAccess.h>
#include <cassert>
#include <vector>

namespace OpenMS
{

  /**
    @brief Estimates the signal/noise (S/N) ratio of each data point in a scan by using the median (window based)

    For each scan, we define a set of windows of a pre-defined size (param:
    <i>window_length_</i>) in m/z domain for which the intensity median is
    calculated.  The noise for a data point is estimated to be the median of the
    intensities of the current window.

    To get a more robust noise estimate, the nose value is calculated two times
    for two sets of windows that are shifted by 1/2 of the window size and the
    reported noise value is the average of the two.

    A call to estimateNoise will return an object of type NoiseEstimator which
    then provides a function get_noise_value which will return the noise value
    for a given m/z value.

    The idea behind this class is to have an estimator for signal to noise that
    gives similar results to SignalToNoiseEstimatorMedian but performs faster.
    Note that it will not give identical results as SignalToNoiseEstimatorMedian
    but for many application the results from this class will be sufficient.

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI SignalToNoiseEstimatorMedianRapid
  {
    /// Window length parameter
    double window_length_;

public:

    /**
      @brief Class to compute the noise value at a given position

      This class implements a method to obtain the noise value at any given m/z
      position. For a median based noise estimator, the noise at position m/z
      is given by the median intensity in a window around this position. This
      noise estimator has median estimates for a set of precomputed windows and
      retrieves the appropriate noise value from the closest window. To lower
      errors at the bin borders, two noise binning values are provided (for a
      set of windows offset by 1/2 of the window width) and the reported value
      is the average of these two values.
    */
    struct OPENMS_DLLAPI NoiseEstimator
    {
      /// Number of windows in m/z direction for which noise values are stored
      int nr_windows;
      /// Start of m/z domain
      double mz_start;
      /// Length of the window in m/z direction
      double window_length;
      /// Noise values for window starting at mz_start (length = nr_windows)
      std::vector<double> result_windows_even;
      /// Noise values for window starting at mz_start - 0.5 * window_length (length = nr_windows + 1)
      std::vector<double> result_windows_odd;

      /// Constructor
      NoiseEstimator() {}

      /// Constructor
      NoiseEstimator(double nr_windows_, double mz_start_, double win_len_) :
        nr_windows(nr_windows_),
        mz_start(mz_start_),
        window_length(win_len_),
        result_windows_even(nr_windows_),
        result_windows_odd(nr_windows_+1)
      {}

      /**
        @brief Return the noise value at a given m/z position

        Will return the noise value at a given m/z position.

        @note Will return 1.0 if the noise would be lower than 1.0
      */
      double get_noise_value (double mz)
      {
        // Take the average of the two stored values
        // Avoid division by 0 (since most clients will divide by the noise value)
        return std::max(1.0, (get_noise_even(mz)+get_noise_odd(mz))/2.0 );
      }

      double get_noise_even (double mz)
      {
        // PRECONDITION
        int window_nr = (int)((mz - mz_start)/window_length);
        assert(window_nr >= 0);
        assert(window_nr < (int)result_windows_even.size());

        double noise = result_windows_even[window_nr];
        return noise;
      }

      double get_noise_odd (double mz)
      {
        // PRECONDITION
        int window_nr = (int)((mz - mz_start + window_length/2.0)/window_length);
        assert(window_nr >= 0);
        assert(window_nr < (int)result_windows_odd.size());

        double noise = result_windows_odd[window_nr];
        return noise;
      }
    };

    /// default constructor
    SignalToNoiseEstimatorMedianRapid(double window_length) :
      window_length_(window_length)
    {
    }

    /** @brief Compute noise estimator for an m/z and intensity array using windows
     *
     * Will return a noise estimator object.
    */
    inline NoiseEstimator estimateNoise(OpenMS::Interfaces::SpectrumPtr spectrum)
    {
      return estimateNoise(spectrum->getMZArray()->data, spectrum->getIntensityArray()->data);
    }

    /** @brief Compute noise estimator for an m/z and intensity array using windows
     *
     * Will return a noise estimator object.
    */
    inline NoiseEstimator estimateNoise(OpenMS::Interfaces::ChromatogramPtr chrom)
    {
      return estimateNoise(chrom->getTimeArray()->data, chrom->getIntensityArray()->data);
    }

    /** @brief Compute noise estimator for an m/z and intensity array using windows
     *
     * Will return a noise estimator object.
    */
    NoiseEstimator estimateNoise(const std::vector<double>& mz_array, const std::vector<double>& int_array)
    {
      // PRECONDITION
      assert(mz_array.size() == int_array.size());
      assert(mz_array.size() > 2);

      int nr_windows = (int)((mz_array[mz_array.size()-1] - mz_array[0])/window_length_) + 1;
      NoiseEstimator eval(nr_windows, mz_array[0], window_length_);

      // Compute even windows
      computeNoiseInWindows_(mz_array, int_array, eval.result_windows_even, mz_array[0]);
      // Compute odd windows
      computeNoiseInWindows_(mz_array, int_array, eval.result_windows_odd, mz_array[0] - window_length_/2.0);

      return eval;
    }

private:

    /** @brief Computes the noise in windows for two input arrays and stores the median intensity in the result (internal)
     *
     * Note that int_array is copied on purpose, since it is modified while sorting, a copy is needed.
     *
    */
    void computeNoiseInWindows_(const std::vector<double>& mz_array, std::vector<double> int_array, std::vector<double> & result, double mz_start);

    /** @brief Median computation on a part of an array [first,last)
     *
     *  @note Does not guarantee that the elements between [first, last) are in the
     *  same order as before (they most likely will not be).
     *
    */
    double computeMedian_(std::vector<double>::iterator & first, std::vector<double>::iterator & last);

  };

} // namespace OpenMS


