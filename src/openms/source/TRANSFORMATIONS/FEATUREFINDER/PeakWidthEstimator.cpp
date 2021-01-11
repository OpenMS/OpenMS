// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>

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
