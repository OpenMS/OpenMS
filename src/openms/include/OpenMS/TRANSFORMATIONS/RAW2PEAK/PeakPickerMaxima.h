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

#ifndef OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERMAXIMA_H
#define OPENMS_TRANSFORMATIONS_RAW2PEAK_PEAKPICKERMAXIMA_H

#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <map>

namespace OpenMS
{

  /**
    @brief This class implements a fast peak-picking algorithm best suited for
    high resolution MS data (FT-ICR-MS, Orbitrap). In high resolution data, the
    signals of ions with similar mass-to-charge ratios (m/z) exhibit little or
    no overlapping and therefore allow for a clear separation. Furthermore, ion
    signals tend to show well-defined peak shapes with narrow peak width.

    This peak-picking algorithm detects ion signals in raw data and
    reconstructs the corresponding peak shape by cubic spline interpolation.
    Signal detection depends on the signal-to-noise ratio which is adjustable
    by the user (see parameter signal_to_noise). A picked peak's m/z and
    intensity value is given by the maximum of the underlying peak spline.

    So far, this peak picker was mainly tested on high resolution data. With
    appropriate preprocessing steps (e.g. noise reduction and baseline
    subtraction), it might be also applied to low resolution data.

    @htmlinclude OpenMS_PeakPickerHiRes.parameters

    @note The peaks must be sorted according to ascending m/z!

    @ingroup PeakPicking
  */
  class OPENMS_DLLAPI PeakPickerMaxima 
  {
public:
    /// Constructor
    PeakPickerMaxima(double signal_to_noise, double spacing_difference = 1.5, double sn_window_length = 200) :
      signal_to_noise_(signal_to_noise),
      spacing_difference_(spacing_difference),
      sn_window_length_(sn_window_length)
    {}

    /// Destructor
    virtual ~PeakPickerMaxima() {}

    /**
      @brief The PeakCandidate describes the output of the peak picker

      It contains the m/z and intensity value of the peak candidate.

      It also contains the original index in the m/z axis where the peak was
      found as well as an estimate of its right and left boundary. 
    */
    struct PeakCandidate 
    {
      /// index of the peak apex (relative to the input data) 
      int pos;
      /// index of the left boundary (relative to the input data) 
      int left_boundary;
      /// index of the right boundary (relative to the input data) 
      int right_boundary;
      /// m/z value of the peak apex
      double mz_max;
      /// intensity value of the peak apex
      double int_max;
    };

    /**
      @brief Will find local maxima in raw data

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates

      @note This function will directly report peak apices with right and left
      boundaries but will not use any fitting to estimate the true m/z and
      intensity of the peak. Note that the mz_max and int_max fields will be
      empty in the result (set to -1).

    */
    void findMaxima(std::vector<double>& mz_array, std::vector<double>& int_array, std::vector<PeakCandidate>& pc)
    {
      if (mz_array.size() < 5) return;

      SignalToNoiseEstimatorMedianRapid::NoiseEstimator noise_estimator(0,0,0);
      if (signal_to_noise_ > 0.0)
      {
        SignalToNoiseEstimatorMedianRapid rapid_sne(sn_window_length_);
        noise_estimator = rapid_sne.estimateNoise(mz_array, int_array);
      }

      // find local maxima in raw data
      for (Size i = 2; i < mz_array.size() - 2; ++i)
      {
        double central_peak_mz = mz_array[i], central_peak_int = int_array[i];
        double left_neighbor_mz = mz_array[i - 1], left_neighbor_int = int_array[i - 1];
        double right_neighbor_mz = mz_array[i + 1], right_neighbor_int = int_array[i + 1];

        // MZ spacing sanity checks
        double left_to_central = std::fabs(central_peak_mz - left_neighbor_mz);
        double central_to_right = std::fabs(right_neighbor_mz - central_peak_mz);
        double min_spacing = (left_to_central < central_to_right) ? left_to_central : central_to_right;

        double act_snt = 0.0, act_snt_l1 = 0.0, act_snt_r1 = 0.0;

        if (signal_to_noise_ > 0.0)
        {
          act_snt = central_peak_int / noise_estimator.get_noise_value(central_peak_mz);
          act_snt_l1 = left_neighbor_int / noise_estimator.get_noise_value(left_neighbor_mz);
          act_snt_r1 = right_neighbor_int / noise_estimator.get_noise_value(right_neighbor_mz);
        }

        // look for peak cores meeting MZ and intensity/SNT criteria
        if (act_snt >= signal_to_noise_
           && left_to_central < spacing_difference_ * min_spacing
           && central_peak_int > left_neighbor_int
           && act_snt_l1 >= signal_to_noise_
           && central_to_right < spacing_difference_ * min_spacing
           && central_peak_int > right_neighbor_int
           && act_snt_r1 >= signal_to_noise_)
        {
          // special case: if a peak core is surrounded by more intense
          // satellite peaks (indicates oscillation rather than
          // real peaks) -> remove

          double act_snt_l2 = 0.0, act_snt_r2 = 0.0;
          PeakCandidate candidate;
          candidate.pos = i;
          candidate.mz_max = -1;
          candidate.int_max = -1;

          if (signal_to_noise_ > 0.0)
          {
            act_snt_l2 = int_array[i-2] / noise_estimator.get_noise_value(mz_array[i-2]);
            act_snt_r2 = int_array[i+2] / noise_estimator.get_noise_value(mz_array[i+2]);
          }

          if ((i > 1
              && std::fabs(left_neighbor_mz - mz_array[i - 2]) < spacing_difference_ * min_spacing
              && left_neighbor_int < int_array[i - 2]
              && act_snt_l2 >= signal_to_noise_)
             &&
              ((i + 2) < mz_array.size()
              && std::fabs(mz_array[i + 2] - right_neighbor_mz) < spacing_difference_ * min_spacing
              && right_neighbor_int < int_array[i + 2]
              && act_snt_r2 >= signal_to_noise_)
              )
          {
            ++i;
            continue;
          }

          std::map<double, double> peak_raw_data;

          peak_raw_data[central_peak_mz] = central_peak_int;
          peak_raw_data[left_neighbor_mz] = left_neighbor_int;
          peak_raw_data[right_neighbor_mz] = right_neighbor_int;

          // peak core found, now extend it
          // to the left
          Size k = 2;

          Size missing_left(0);
          Size missing_right(0);

          while ((i - k + 1) > 0
                && (missing_left < 2)
                && int_array[i - k] <= peak_raw_data.begin()->second)
          {

            double act_snt_lk = 0.0;
            if (signal_to_noise_ > 0.0)
            {
              act_snt_lk = int_array[i-k] / noise_estimator.get_noise_value(mz_array[i-k]);
            }

            if (act_snt_lk >= signal_to_noise_ && std::fabs(mz_array[i - k] - peak_raw_data.begin()->first) < spacing_difference_ * min_spacing)
            {
              peak_raw_data[mz_array[i - k]] = int_array[i - k];
            }
            else
            {
              peak_raw_data[mz_array[i - k]] = int_array[i - k];
              ++missing_left;
            }
            ++k;
          }
          candidate.left_boundary = i - k + 1;

          // to the right
          k = 2;
          while ((i + k) < mz_array.size()
                && (missing_right < 2)
                && int_array[i + k] <= peak_raw_data.rbegin()->second)
          {
            double act_snt_rk = 0.0;
            if (signal_to_noise_ > 0.0)
            {
              act_snt_rk = int_array[i+k] / ne.get_noise_value(mz_array[i+k]);
            }

            if (act_snt_rk >= signal_to_noise_ && std::fabs(mz_array[i + k] - peak_raw_data.rbegin()->first) < spacing_difference_ * min_spacing)
            {
              peak_raw_data[mz_array[i + k]] = int_array[i + k];
            }
            else
            {
              peak_raw_data[mz_array[i + k]] = int_array[i + k];
              ++missing_right;
            }
            ++k;
          }
          candidate.right_boundary = i + k - 1;
          // jump over raw data points that have been considered already
          i = i + k - 1;
          pc.push_back(candidate);
        }
      }
    }

    /**
      @brief Will pick peaks in a spectrum

      @param mz_array The array containing m/z values
      @param int_array The array containing intensity values
      @param pc The resulting array containing the peak candidates

      @note This function will first find maxima in the intensity domain and
      then use a spline function to estimate the best m/z and intensity for
      each peak candidate.
    */
    void pick(std::vector<double>& mz_array, std::vector<double>& int_array, std::vector<PeakCandidate>& pc)
    {
      if (mz_array.size() < 5) return;

      findMaxima(mz_array, int_array, pc);

      // Go through all peak candidates and find accurate mz / int values based on the spline interpolation
      for (Size j = 0; j < pc.size(); ++j)
      {
          PeakCandidate candidate = pc[j];

          // output all raw data points selected for one peak
          // TODO: #ifdef DEBUG_ ...
          // for (std::map<double, double>::const_iterator map_it = peak_raw_data.begin(); map_it != peak_raw_data.end(); ++map_it) {
          // PeakType peak;
          // peak.setMZ(map_it->first);
          // peak.setIntensity(map_it->second);
          // output.push_back(peak);
          // std::cout << map_it->first << " " << map_it->second << " snt: " << std::endl;
          // }
          // std::cout << "--------------------" << std::endl;

          double central_peak_mz = mz_array[candidate.pos], central_peak_int = int_array[candidate.pos];
          double left_neighbor_mz = mz_array[candidate.pos - 1]; //, left_neighbor_int = int_array[candidate.pos - 1];
          double right_neighbor_mz = mz_array[candidate.pos + 1]; //, right_neighbor_int = int_array[candidate.pos + 1];


          std::vector<double> raw_mz_values;
          std::vector<double> raw_int_values;

          raw_mz_values.reserve(candidate.right_boundary - candidate.left_boundary);
          raw_int_values.reserve(candidate.right_boundary - candidate.left_boundary);

          raw_mz_values.insert(raw_mz_values.begin(), mz_array.begin() + candidate.left_boundary, mz_array.begin() + candidate.right_boundary + 1);
          raw_int_values.insert(raw_int_values.begin(), int_array.begin() + candidate.left_boundary, int_array.begin() + candidate.right_boundary + 1);

          const Size num_raw_points = raw_mz_values.size();

          // setup gsl splines
          gsl_interp_accel * spline_acc = gsl_interp_accel_alloc();
          gsl_interp_accel * first_deriv_acc = gsl_interp_accel_alloc();
          gsl_spline * peak_spline = gsl_spline_alloc(gsl_interp_cspline, num_raw_points);
          gsl_spline_init(peak_spline, &(*raw_mz_values.begin()), &(*raw_int_values.begin()), num_raw_points);


          // calculate maximum by evaluating the spline's 1st derivative
          // (bisection method)
          double max_peak_mz = central_peak_mz, max_peak_int = central_peak_int;
          double threshold = 0.000001;
          double lefthand = left_neighbor_mz;
          double righthand = right_neighbor_mz;

          bool lefthand_sign = 1;
          double eps = std::numeric_limits<double>::epsilon();

          // bisection
          do
          {
            double mid = (lefthand + righthand) / 2;

            double midpoint_deriv_val = gsl_spline_eval_deriv(peak_spline, mid, first_deriv_acc);

            // if deriv nearly zero then maximum already found
            if (!(std::fabs(midpoint_deriv_val) > eps))
            {
              break;
            }

            bool midpoint_sign = (midpoint_deriv_val < 0.0) ? 0 : 1;

            if (lefthand_sign ^ midpoint_sign)
            {
              righthand = mid;
            }
            else
            {
              lefthand = mid;
            }

            // TODO: #ifdef DEBUG_ ...
            // PeakType peak;
            // peak.setMZ(mid);
            // peak.setIntensity(gsl_spline_eval(peak_spline, mid, spline_acc));
            // output.push_back(peak);

          }
          while (std::fabs(lefthand - righthand) > threshold);

          // sanity check?
          max_peak_mz = (lefthand + righthand) / 2;
          max_peak_int = gsl_spline_eval(peak_spline, max_peak_mz, spline_acc);

          // save picked pick into output spectrum
          pc[j].mz_max = max_peak_mz;
          pc[j].int_max = max_peak_int;

          // free allocated gsl memory
          gsl_spline_free(peak_spline);
          gsl_interp_accel_free(spline_acc);
          gsl_interp_accel_free(first_deriv_acc);
      }
    }

protected:
    // signal-to-noise parameter
    double signal_to_noise_;

    // maximal spacing difference
    double spacing_difference_;

    double sn_window_length_;

  }; // end PeakPickerMaxima

} // namespace OpenMS

#endif
