// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_TRANSFORMATIONS_PEAKWIDTHESTIMATOR_H
#define OPENMS_TRANSFORMATIONS_PEAKWIDTHESTIMATOR_H


#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerHiRes.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_interp.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>
#include <deque>

namespace OpenMS
{
/**
   @brief This class implements a peak width estimation algorithm best suited for high resolution MS data (FT-ICR-MS, Orbitrap).
          Peaks are detected and a spline is fitted to the raw data in a window around the peak.
          Then a search for to the half-maximum is performed on the spline to the left and right of the peak maximum.
          The Full Width at the Half Maximum is collected.
          Finally a linear regression is performed to determine FWHM(m/z)

   @note The peaks must be sorted according to ascending m/z!

   @experimental This algorithm has not been tested thoroughly yet.
  */
class OPENMS_DLLAPI PeakWidthEstimator
{
public:
  /// Constructor
  PeakWidthEstimator();

  /// Destructor
  virtual ~PeakWidthEstimator();

  template <typename PeakType>
  void estimateFWHM(const MSSpectrum<PeakType>& input, std::multimap<DoubleReal, DoubleReal>& fwhms) const
  {
    MSSpectrum<PeakType> picked;
    // 1. find peaks
    peak_picker_.pick(input, picked);
    // 2. determine fwhm
    for (Size i = 1; i < picked.size() - 1; ++i)
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

      typename MSSpectrum<PeakType>::ConstIterator begin_window = input.MZBegin(picked[i].getMZ() - half_window_size);
      typename MSSpectrum<PeakType>::ConstIterator end_window = input.MZBegin(picked[i].getMZ() + half_window_size);

      std::vector<double> raw_mz_values;
      std::vector<double> raw_int_values;

      raw_mz_values.push_back(mz - 0.3);
      raw_int_values.push_back(0);

      raw_mz_values.push_back(mz - 0.25);
      raw_int_values.push_back(0);

      bool max_peak_added = false;
      for (; begin_window != end_window; ++begin_window)
      {
        if (!max_peak_added && begin_window->getMZ() > mz)
        {
          raw_mz_values.push_back(mz);
          raw_int_values.push_back(intensity);
          max_peak_added = true;
        }
        raw_mz_values.push_back(begin_window->getMZ());
        raw_int_values.push_back(end_window->getIntensity());
      }

      raw_mz_values.push_back(mz + 0.25);
      raw_int_values.push_back(0);

      raw_mz_values.push_back(mz + 0.3);
      raw_int_values.push_back(0);

      if (raw_mz_values.size() < 12)
      {
        continue;
      }

      // setup gsl splines
      const Size num_raw_points = raw_mz_values.size();
      gsl_interp_accel *spline_acc = gsl_interp_accel_alloc();
      gsl_interp_accel *first_deriv_acc = gsl_interp_accel_alloc();
      gsl_spline *peak_spline = gsl_spline_alloc(gsl_interp_cspline, num_raw_points);
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

      // IMPROVEMENT better outlier filtering using RANSAC on the resulting data
      if ( fwhm > 0.1 )
      {
        continue;
      }

      fwhms.insert( std::pair<DoubleReal, DoubleReal>(mz, fwhm) );
    }
  }

  template <typename PeakType>
  void estimateFWHM(const MSExperiment<PeakType>& input, DoubleReal& intercept, DoubleReal& slope) const
  {
    MSExperiment<PeakType> exp;

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
    } else
    {
      while (ms1_indices_size > 100)
      {
        if (ms1_indices_size % 2 == 1)
        {
          ms1_indices.pop_front();
        } else
        {
          ms1_indices.pop_back();
        }
        ms1_indices_size--;
      }
    }

    // build reduced experiment
    for (Size c = 0; c != ms1_indices_size; ++c)
    {
      MSSpectrum<PeakType> spectrum = input[ms1_indices[c]];
      exp.push_back(spectrum);
    }

    std::multimap<DoubleReal, DoubleReal> fwhms; // map peak mz to fwhm

    // estimate FWHM on every spectrum
    for (Size scan_idx = 0; scan_idx != exp.size(); ++scan_idx)
    {
      estimateFWHM(exp[scan_idx], fwhms);
    }

    // extract mzs and fwhm for linear regression
    std::vector<DoubleReal> keys;
    std::vector<DoubleReal> values;
    for (std::multimap<DoubleReal, DoubleReal>::iterator it = fwhms.begin(); it != fwhms.end(); ++it)
    {
      // std::cout << it->first << " " << it->second << std::endl;
      keys.push_back(it->first);
      values.push_back(it->second);
    }

    Math::LinearRegression linear_reg;
    linear_reg.computeRegression(0.95, keys.begin(), keys.end(), values.begin());

    slope = linear_reg.getSlope();
    intercept = linear_reg.getIntercept();
    // std::cout << "intercept: " << intercept << " slope: " << slope << std::endl;
    // std::cout << linear_reg.getRSD() << std::endl;
    return;
  }

  protected:
    PeakPickerHiRes peak_picker_;

};

}

#endif // OPENMS_TRANSFORMATIONS_PEAKWIDTHESTIMATOR_H
