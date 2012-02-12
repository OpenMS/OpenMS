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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------


#ifndef OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H
#define OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
  /**
    @brief LOWESS (locally weighted scatterplot smoothing).

    A smoothing technique that fits simple models (linear, quadratic) to localized subsets of the data, point by point.
    This is particularly useful for smoothing intensities in spectra or chromatograms. In this case, the window size for the smoothing
    should be setted proportional to the peak width (see LowessSmoothing parameters).

    @htmlinclude OpenMS_LowessSmoothing.parameters

    @ingroup SignalProcessing
  */
  class OPENMS_DLLAPI LowessSmoothing
    : public DefaultParamHandler
  {
  public:
    /// Default constructor
    LowessSmoothing();

    /// Destructor
    virtual ~LowessSmoothing();

    typedef std::vector<DoubleReal> DoubleVector;

    /// Smoothing method that receives x and y coordinates (e.g., RT and intensities) and computes smoothed intensities.
    void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&);

  protected:
    virtual void updateMembers_();

  private:
    DoubleReal window_size_;

    DoubleReal tricube_(DoubleReal, DoubleReal);
  };


} // namespace OpenMS
#endif // OPENMS_FILTERING_SMOOTHING_LOWESSSMOOTHING_H
