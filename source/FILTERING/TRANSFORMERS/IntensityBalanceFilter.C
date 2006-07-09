// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: IntensityBalanceFilter.C,v 1.6 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityBalanceFilter.h>

using namespace std;

namespace OpenMS
{
  IntensityBalanceFilter::IntensityBalanceFilter()
    :FilterFunctor()
  { 
		name_ = IntensityBalanceFilter::getName();
  }

  IntensityBalanceFilter::IntensityBalanceFilter(const IntensityBalanceFilter& source )
    : FilterFunctor(source)
  {
  }
  
  IntensityBalanceFilter& IntensityBalanceFilter::operator=(const IntensityBalanceFilter& source)
  {
    FilterFunctor::operator=(source);
    return *this;
  }
  
  IntensityBalanceFilter::~IntensityBalanceFilter()
  {
  }

	/*
  vector<double> IntensityBalanceFilter::operator() (const ClusterSpectrum& cspec)
  {
    double bands = 10;
    multimap<double,uint> bandIntensity;
    double size = cspec.getSpec().getPrecursorPeak().getPosition()[0];
    uint j = 0;
    for (uint i = 0; i < bands; ++i)
    {
      double intensity = 0;
      //bern 2004 says to only check between 300 and size
      //but that seems inappropriate for small peptides (smallest is ca 450)
      while ( j < cspec.getSpec().size() && cspec.getSpec().getContainer()[j].getPosition()[0] < (size-300)/bands*(i+1) +300)
      {
        intensity += cspec.getSpec().getContainer()[j++].getIntensity();
      }
      bandIntensity.insert(make_pair(intensity,i));
    }
    j = 0;
    double totalIntensity = 0;
    double twobiggest = 0;
    double sevensmallest = 0;
    for (multimap<double,uint>::reverse_iterator mmrit = bandIntensity.rbegin(); mmrit != bandIntensity.rend(); ++mmrit,++j)
    {
      totalIntensity += mmrit->first;
      //take the two biggest
      if (j < 2)
      {
        twobiggest+=mmrit->first;
      }
      //take the seven smallest
      if ( j > 2 )
      {
        sevensmallest += mmrit->first;
      }
    }
    vector<double> result;
    result.push_back((twobiggest - sevensmallest)/totalIntensity);
    return result;
  }
	*/
}
