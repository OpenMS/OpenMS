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
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <deque>
#include <boost/tuple/tuple_comparison.hpp>

#include <OpenMS/MATH/gsl_wrapper.h>

namespace OpenMS
{
  void PeakWidthEstimator::estimateSpectrumFWHM(const MSSpectrum<Peak1D> & input, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > & fwhms)
  {
    PeakPickerHiRes picker;

    MSSpectrum<Peak1D> picked;
    // 1. find peaks
    picker.pick(input, picked);
    // 2. determine fwhm

    if (picked.size() < 3)
      return;

    for (Size i = 1; i < picked.size() - 1; ++i)
    {
      DoubleReal mz = picked[i].getMZ();
      DoubleReal intensity = picked[i].getIntensity();

      DoubleReal half_window_size = 0.25;

      DoubleReal half_dist_left_neighbor = (mz - picked[i - 1].getMZ()) / 2.0;
      DoubleReal half_dist_right_neighbor = (picked[i + 1].getMZ() - mz) / 2.0;

      if (half_dist_right_neighbor < half_window_size)
      {
        half_window_size = half_dist_right_neighbor;
      }

      if (half_dist_left_neighbor < half_window_size)
      {
        half_window_size = half_dist_left_neighbor;
      }

      if (half_window_size < 0.01)
      {
        continue;
      }

      MSSpectrum<Peak1D>::ConstIterator begin_window = input.MZBegin(picked[i].getMZ() - half_window_size);
      MSSpectrum<Peak1D>::ConstIterator end_window = input.MZBegin(picked[i].getMZ() + half_window_size);

      std::map<double, double> values;

      // Add peak maximum
      values.insert(std::make_pair(mz, intensity));

      for (; begin_window != end_window; ++begin_window)
      {
        values.insert(std::make_pair(begin_window->getMZ(), begin_window->getIntensity()));
      }

      // Make sure we have some zeroes
      values.insert(std::make_pair(mz - 0.3, 0));
      values.insert(std::make_pair(mz - 0.25, 0));
      values.insert(std::make_pair(mz + 0.25, 0));
      values.insert(std::make_pair(mz + 0.3, 0));

      if (values.size() < 12)
      {
        continue;
      }

      std::vector<double> raw_mz_values, raw_int_values;
      raw_mz_values.reserve(values.size());
      raw_int_values.reserve(values.size());
      for (std::map<double, double>::const_iterator it = values.begin(); it != values.end(); ++it)
      {
        raw_mz_values.push_back(it->first);
        raw_int_values.push_back(it->second);
      }

      // setup gsl splines
      const Size num_raw_points = raw_mz_values.size();
      deprecated_gsl_interp_accel * spline_acc = deprecated_gsl_interp_accel_alloc();
      deprecated_gsl_interp_accel * first_deriv_acc = deprecated_gsl_interp_accel_alloc();
      deprecated_gsl_spline * peak_spline = deprecated_gsl_spline_alloc(deprecated_wrapper_get_gsl_interp_akima(), num_raw_points);
      deprecated_gsl_spline_init(peak_spline, &(*raw_mz_values.begin()), &(*raw_int_values.begin()), num_raw_points);

      // search for half intensity to the left
      DoubleReal MZ_THRESHOLD = 0.00000001;
      DoubleReal INTENSITY_THRESHOLD = 0.0001;

      DoubleReal half_maximum = intensity / 2.0;
      DoubleReal mid;

      // left bisection
      DoubleReal left = mz - 0.25;
      DoubleReal right = mz;

      do
      {
        mid = (left + right) / 2.0;
        DoubleReal mid_int = deprecated_gsl_spline_eval(peak_spline, mid, spline_acc);

        // found half maximum
        if (std::fabs(mid_int - half_maximum) < INTENSITY_THRESHOLD)
        {
          break;
        }

        if (mid_int < half_maximum)   // mid < half maximum ?
        {
          left = mid;
        }
        else
        {
          right = mid;
        }

        // std::cout << "L~" << mz << " : " << mid << " " << mid_int << " # " << half_maximum << std::endl;
      }
      while (std::fabs(left - right) > MZ_THRESHOLD);

      DoubleReal left_fwhm_mz = mid;

      // right bisection
      left = mz;
      right = mz + 0.25;
      do
      {
        mid = (left + right) / 2.0;
        DoubleReal mid_int = deprecated_gsl_spline_eval(peak_spline, mid, spline_acc);

        // found half maximum
        if (std::fabs(mid_int - half_maximum) < INTENSITY_THRESHOLD)
        {
          break;
        }

        if (mid_int > half_maximum)   // mid < half maximum ?
        {
          left = mid;
        }
        else
        {
          right = mid;
        }

        // std::cout << "R~" << mz << " : " << mid << " " << mid_int << " # " << half_maximum << std::endl;
      }
      while (std::fabs(left - right) > MZ_THRESHOLD);

      DoubleReal right_fwhm_mz = mid;

      // free allocated gsl memory
      deprecated_gsl_spline_free(peak_spline);
      deprecated_gsl_interp_accel_free(spline_acc);
      deprecated_gsl_interp_accel_free(first_deriv_acc);

      // sanity check (left distance and right distance should be more or less equal)
      DoubleReal ratio = std::fabs(left_fwhm_mz - mz) / std::fabs(right_fwhm_mz - mz);

      if (ratio < 0.9 || ratio > 1.1)
      {
        continue;
      }

      //std::cout << std::fabs(left_fwhm_mz - mz) << " " << std::fabs(right_fwhm_mz - mz) << std::endl;

      DoubleReal fwhm = std::fabs(left_fwhm_mz - mz) + std::fabs(right_fwhm_mz - mz);

      fwhms.insert(boost::make_tuple(intensity, mz, fwhm));
    }
  }

  PeakWidthEstimator::Result PeakWidthEstimator::estimateFWHM(const MSExperiment<Peak1D> & input)
  {
    MSExperiment<Peak1D> exp;

    // extract ms1 indices
    std::deque<Size> ms1_indices;
    Size ms1_indices_size = 0;
    for (Size c = 0; c != input.size(); ++c)
    {
      if (input[c].getMSLevel() == 1)
      {
        ms1_indices.push_back(c);
        ms1_indices_size++;
      }
    }

    // reduce to max 100 spectra. possible improvement: use spectra with high intensity in TIC
    if (ms1_indices_size <= 100)
    {
      exp = input;
    }
    else
    {
      while (ms1_indices_size > 100)
      {
        if (ms1_indices_size % 2 == 1)
        {
          ms1_indices.pop_front();
        }
        else
        {
          ms1_indices.pop_back();
        }
        --ms1_indices_size;
      }
    }

    // build reduced experiment
    for (Size c = 0; c != ms1_indices_size; ++c)
    {
      MSSpectrum<Peak1D> spectrum = input[ms1_indices[c]];
      exp.addSpectrum(spectrum);
    }

    // set of (intensity, mz, peak-width)
    std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > fwhms;

    // estimate FWHM on every spectrum
    for (Size scan_idx = 0; scan_idx != exp.size(); ++scan_idx)
    {
      estimateSpectrumFWHM(exp[scan_idx], fwhms);
    }

    if (fwhms.empty())
    {
      throw Exception::InvalidSize(__FILE__, __LINE__, __PRETTY_FUNCTION__, fwhms.size());
    }

    // extract mzs and fwhm for linear regression above the median sorted for the intensity
    std::vector<double> keys, values, weights;
    {
      Size count = fwhms.size() / 2;
      std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> >::reverse_iterator it = fwhms.rbegin();
      for (; count && it != fwhms.rend(); --count, ++it)
      {
        // std::cout << it->get<1>() << ',' << it->get<2>() << ',' << it->get<0>() << std::endl;  // generates nice plots
        keys.push_back(std::log(it->get<1>()));
        values.push_back(std::log(it->get<2>()));
        weights.push_back(it->get<0>());
      }
    }

    double c0, c1, cov00, cov01, cov11, chisq;
    int error = deprecated_gsl_fit_wlinear(&keys[0], 1, &weights[0], 1, &values[0], 1, keys.size(),
                                &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    if (error)
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-PeakWidthEstimator", "Error from GSL");
    }

    return Result(c0, c1);
  }

}
