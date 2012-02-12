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

          Result()
            : c0(0), c1(0)
          { }

          Result(const DoubleReal c0, const DoubleReal c1)
            : c0(c0), c1(c1)
          { }

          DoubleReal operator()(const DoubleReal mz) const
          {
            return std::exp(c0 + c1 * std::log(mz));
          }
      };

      static void estimateSpectrumFWHM(const MSSpectrum<Peak1D>& input, std::set<boost::tuple<DoubleReal, DoubleReal, DoubleReal> >& fwhms);
      static Result estimateFWHM(const MSExperiment<Peak1D>& input);
  };
}

#endif // OPENMS_TRANSFORMATIONS_PEAKWIDTHESTIMATOR_H
