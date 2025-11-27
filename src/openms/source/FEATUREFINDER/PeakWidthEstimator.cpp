// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/FEATUREFINDER/PeakWidthEstimator.h>

namespace OpenMS
{

  PeakWidthEstimator::PeakWidthEstimator(const PeakMap & exp_picked, const std::vector<std::vector<PeakPickerHiRes::PeakBoundary> > & boundaries)
  {
    std::vector<double> peaks_mz;
    std::vector<double> peaks_width;
    PeakMap::ConstIterator it_rt;
    std::vector<std::vector<PeakPickerHiRes::PeakBoundary> >::const_iterator it_rt_boundaries;
    for (it_rt = exp_picked.begin(), it_rt_boundaries = boundaries.begin();
         it_rt < exp_picked.end() && it_rt_boundaries < boundaries.end();
         ++it_rt, ++it_rt_boundaries)
    {
      MSSpectrum::ConstIterator it_mz;
      std::vector<PeakPickerHiRes::PeakBoundary>::const_iterator it_mz_boundary;
      for (it_mz = it_rt->begin(), it_mz_boundary = it_rt_boundaries->begin();
           it_mz < it_rt->end() && it_mz_boundary < it_rt_boundaries->end();
           ++it_mz, ++it_mz_boundary)
      {
          peaks_mz.push_back(it_mz->getMZ());
          peaks_width.push_back((*it_mz_boundary).mz_max - (*it_mz_boundary).mz_min);
      }
    }

    mz_min_ = peaks_mz.front();
    mz_max_ = peaks_mz.back();
    bspline_ = new BSpline2d(peaks_mz, peaks_width, std::min(500.0, (mz_max_ - mz_min_)/2), BSpline2d::BC_ZERO_SECOND, 1);
      
    if (!(*bspline_).ok())
    {
      throw Exception::UnableToFit(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unable to fit B-spline to data.", "");
    }
  }
  
  PeakWidthEstimator::~PeakWidthEstimator()
  {
    delete bspline_;
  }
  
  double PeakWidthEstimator::getPeakWidth(double mz)
  {
    double width;

    if (mz < mz_min_)
    {
      width = (*bspline_).eval(mz_min_);
    }
    else if (mz > mz_max_)
    {
      width = (*bspline_).eval(mz_max_);
    }
    else
    {
      width = (*bspline_).eval(mz);
    }

    if (width < 0)
    {
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Estimated peak width is negative.", "");
    }

    return width;
  }

}
