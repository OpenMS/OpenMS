// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKWIDTHESTIMATOR_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_PEAKWIDTHESTIMATOR_H


#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

#include <boost/tuple/tuple.hpp>

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
    class Result
    {
public:
      DoubleReal c0, c1;

      Result() :
        c0(0), c1(0)
      {}

      Result(const DoubleReal c0, const DoubleReal c1) :
        c0(c0), c1(c1)
      {}

      DoubleReal operator()(const DoubleReal mz) const
      {
        return std::exp(c0 + c1 * std::log(mz));
      }

    };

    static void estimateSpectrumFWHM(const MSSpectrum<Peak1D> & input, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> > & fwhms);
    static Result estimateFWHM(const MSExperiment<Peak1D> & input);
  };
}

#endif // OPENMS_TRANSFORMATIONS_PEAKWIDTHESTIMATOR_H
