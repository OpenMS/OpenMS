// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <deque>
#include <boost/tuple/tuple_comparison.hpp>

#include <gsl/gsl_fit.h>
#include <gsl/gsl_interp.h>
#include <gsl/gsl_spline.h>

namespace OpenMS
{
  void PeakWidthEstimator::estimateSpectrumFWHM(const MSSpectrum<Peak1D>& input, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> >& fwhms)
  {
    PeakPickerHiRes picker;

    MSSpectrum<Peak1D> picked;
    // 1. find peaks
    picker.pick(input, picked);
    // 2. determine fwhm

    if (picked.size() < 3) return;

    for (Size i = 1; i < picked.size()-1; ++i)
    {
      DoubleReal mz = picked[i].getMZ();
      DoubleReal intensity = picked[i].getIntensity();

      DoubleReal half_window_size = 0.25;

      DoubleReal half_dist_left_neighbor = (mz - picked[i-1].getMZ()) / 2.0;
      DoubleReal half_dist_right_neighbor = (picked[i+1].getMZ() - mz) / 2.0;

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
      gsl_interp_accel *spline_acc = gsl_interp_accel_alloc();
      gsl_interp_accel *first_deriv_acc = gsl_interp_accel_alloc();
      gsl_spline *peak_spline = gsl_spline_alloc(gsl_interp_akima, num_raw_points);
      gsl_spline_init(peak_spline, &(*raw_mz_values.begin()), &(*raw_int_values.begin()), num_raw_points);

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
        DoubleReal mid_int = gsl_spline_eval(peak_spline, mid, spline_acc);

        // found half maximum
        if (std::fabs(mid_int - half_maximum) < INTENSITY_THRESHOLD)
        {
          break;
        }

        if ( mid_int < half_maximum ) // mid < half maximum ?
        {
          left = mid;
        } else
        {
          right = mid;
        }

        // std::cout << "L~" << mz << " : " << mid << " " << mid_int << " # " << half_maximum << std::endl;
      } while ( std::fabs(left - right) > MZ_THRESHOLD );

      DoubleReal left_fwhm_mz = mid;

      // right bisection
      left = mz;
      right = mz + 0.25;
      do
      {
        mid = (left + right) / 2.0;
        DoubleReal mid_int = gsl_spline_eval(peak_spline, mid, spline_acc);

        // found half maximum
        if (std::fabs(mid_int - half_maximum) < INTENSITY_THRESHOLD)
        {
          break;
        }

        if ( mid_int > half_maximum ) // mid < half maximum ?
        {
          left = mid;
        } else
        {
          right = mid;
        }

        // std::cout << "R~" << mz << " : " << mid << " " << mid_int << " # " << half_maximum << std::endl;
      } while ( std::fabs(left - right) > MZ_THRESHOLD );

      DoubleReal right_fwhm_mz = mid;

      // free allocated gsl memory
      gsl_spline_free(peak_spline);
      gsl_interp_accel_free(spline_acc);
      gsl_interp_accel_free(first_deriv_acc);

      // sanity check (left distance and right distance should be more or less equal)
      DoubleReal ratio = std::fabs(left_fwhm_mz - mz) / std::fabs(right_fwhm_mz - mz);

      if ( ratio < 0.9 || ratio > 1.1)
      {
        continue;
      }

      //std::cout << std::fabs(left_fwhm_mz - mz) << " " << std::fabs(right_fwhm_mz - mz) << std::endl;

      DoubleReal fwhm = std::fabs(left_fwhm_mz - mz) + std::fabs(right_fwhm_mz - mz);

      fwhms.insert(boost::make_tuple(intensity, mz, fwhm));
    }
  }

  PeakWidthEstimator::Result PeakWidthEstimator::estimateFWHM(const MSExperiment<Peak1D>& input)
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
      exp.push_back(spectrum);
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
    int error = gsl_fit_wlinear(&keys[0], 1, &weights[0], 1, &values[0], 1, keys.size(),
                                &c0, &c1, &cov00, &cov01, &cov11, &chisq);

    if (error)
    {
      throw Exception::UnableToFit( __FILE__, __LINE__, __PRETTY_FUNCTION__, "UnableToFit-PeakWidthEstimator", "Error from GSL");
    }

    return Result(c0, c1);
  }
}
