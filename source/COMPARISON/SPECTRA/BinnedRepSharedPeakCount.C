// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSharedPeakCount.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <cmath>

using namespace std;

namespace OpenMS
{

  BinnedRepSharedPeakCount::BinnedRepSharedPeakCount()
		: BinnedRepCompareFunctor()
  {
		setName(BinnedRepSharedPeakCount::getProductName());
		defaultsToParam_();
  }

  BinnedRepSharedPeakCount::BinnedRepSharedPeakCount(const BinnedRepSharedPeakCount& source)
    : BinnedRepCompareFunctor(source)
  {
  }

  BinnedRepSharedPeakCount& BinnedRepSharedPeakCount::operator= (const BinnedRepSharedPeakCount& source)
  {
		if (&source != this)
		{
	    BinnedRepCompareFunctor::operator=(source);
		}
    return *this;
  }
  
  BinnedRepSharedPeakCount::~BinnedRepSharedPeakCount()
  {
  }

	double BinnedRepSharedPeakCount::operator () (const BinnedRep& a) const
	{
		return operator () (a, a);
	}
	
  double BinnedRepSharedPeakCount::operator()(const BinnedRep& a, const BinnedRep& b)const
  {
    double cutoff = 0;
    //const BinnedRep& a = csa.getBinrep();
    //const BinnedRep& b = csb.getBinrep();
    //double filterfactor = filter(csa,csb);
    //if ( filterfactor < 1e-12 ) return 0;
    double similarity = 0;
    BinnedRep::const_iterator ait = a.begin();
    BinnedRep::const_iterator bit = b.begin();
    while (ait != a.end() && bit != b.end())
    {
      //we are at the same position (+- precision)
      if ( fabs( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  < 1e-8)
      {
        if ( *ait > cutoff && *bit > cutoff ) 
        {
          similarity++;
        }
        ait.hop();
        bit.hop();
      }
      //ait lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  <= 1e-8)
      {
        ait.hop();
      }
      //bit lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  >= 1e-8)
      {
        bit.hop();
      } 
    }
    return similarity;
  }

}
