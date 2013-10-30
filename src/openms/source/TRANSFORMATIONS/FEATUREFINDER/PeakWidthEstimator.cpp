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
#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

#include <boost/tuple/tuple_comparison.hpp>

#include "Wm5IntpAkimaNonuniform1.h"

#include <deque>
#include <algorithm>

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

      typedef std::pair<double, double> MzIntPair; //mz-intensity pair
      std::vector< MzIntPair > values;

      // Add peak maximum
      values.push_back(MzIntPair(mz, intensity));

      for (; begin_window != end_window; ++begin_window)
      {
        values.push_back(MzIntPair(begin_window->getMZ(), begin_window->getIntensity()));
      }
      int true_data_points = 0;
      std::vector< MzIntPair >::const_iterator iter;
      for(iter = values.begin(); iter != values.end(); ++iter)
      {
        double mz = iter->first;
        if(mz > std::numeric_limits<double>::epsilon())//null-point
          ++true_data_points;
      }
      //abort if we have less than 4 data points for spline fitting
      if(true_data_points < 4)
        continue;

      // Make sure we have some zeroes
      values.push_back(MzIntPair(mz - 0.3, 0));
      values.push_back(MzIntPair(mz - 0.25, 0));
      values.push_back(MzIntPair(mz + 0.25, 0));
      values.push_back(MzIntPair(mz + 0.3, 0));

      //make sure the data points are sorted by mz values
      std::sort( values.begin(), values.end() );

      //TODO: heavy code-duplication! Compare PeakPickerHiRes
      //convert into GeomTools Matrix structure
      std::vector<double> X;
      std::vector<double> Y;
      std::vector< std::pair<double, double> >::const_iterator it;
      for (it = values.begin(); it != values.end(); ++it)
      {
        X.push_back( it->first );//mz
        Y.push_back( it->second );//intensity
      }
      Wm5::IntpAkimaNonuniform1<double> peak_spline (values.size(), X.data(), Y.data());

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
        DoubleReal mid_int = peak_spline( mid );

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
        DoubleReal mid_int = peak_spline( mid );

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
    std::vector<double> mzsVec, fwhmsVec, intensitiesVec;
    {
      Size count = fwhms.size() / 2;
      std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> >::reverse_iterator it = fwhms.rbegin();
      for (; count && it != fwhms.rend(); --count, ++it)
      {
//        std::cout << it->get<1>() << ',' << it->get<2>() << ',' << it->get<0>() << std::endl;  // generates nice plots
        mzsVec.push_back(std::log(it->get<1>()));
        fwhmsVec.push_back(std::log(it->get<2>()));
        intensitiesVec.push_back(it->get<0>());
      }
    }

    Math::LinearRegression linreg;
    linreg.computeRegressionWeighted(0.95, mzsVec.begin(), mzsVec.end(), fwhmsVec.begin(), intensitiesVec.begin());

    return Result(linreg.getIntercept(), linreg.getSlope());
  }

}
