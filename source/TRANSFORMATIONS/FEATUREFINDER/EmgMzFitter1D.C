// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Rene Hussong $
// --------------------------------------------------------------------------

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/EmgMzFitter1D.h>
#include <OpenMS/MATH/STATISTICS/StatisticFunctions.h>

namespace OpenMS
{
    EmgMzFitter1D::EmgMzFitter1D()
    : EmgFitter1D()
    {
    }

    EmgMzFitter1D::EmgMzFitter1D(const EmgMzFitter1D& source)
    : EmgFitter1D (source)
    {
    }

    EmgMzFitter1D::~EmgMzFitter1D()
    {
    }

    EmgMzFitter1D& EmgMzFitter1D::operator = (const EmgMzFitter1D& source)
    {
        if (&source == this) return *this;

				EmgFitter1D::operator = (source);

        return *this;
    }
		
              
    void EmgMzFitter1D::setInitialParameters_(const RawDataArrayType& set)
    {
      // sum over all intensities
      CoordinateType sum = 0.0;
      for (UInt i=0; i<set.size(); ++i) sum += set[i].getIntensity();

      // calculate the median
      Int median = 0;
      Real count = 0.0;
      for ( UInt i = 0; i < set.size(); ++i )
      {
        count += set[i].getIntensity();
        if ( count <= sum / 2 ) median = i;
      }
 
      // calculate the height of the peak
      height_ = set[median].getIntensity();

      // calculate retention time
      retention_ = set[median].getPos();
      
      // default is an asymmetric peak
      symmetric_ = false;

			//The average MZ spacing is a good starting point
			DoubleReal mean=0;
			for (UInt i=0; i<set.size()-1; ++i)
			{
				mean += set[i+1].getMZ()-set[i].getMZ();
			};

			mean /= (DoubleReal) (set.size()-1.);

      symmetry_ = 3*mean;
      width_ = 3*mean;
    }
}
