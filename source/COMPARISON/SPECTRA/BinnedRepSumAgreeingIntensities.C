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
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  BinnedRepSumAgreeingIntensities::BinnedRepSumAgreeingIntensities()
		: BinnedRepCompareFunctor()
  {
		setName(BinnedRepSumAgreeingIntensities::getProductName());
		defaultsToParam_();
  }

  BinnedRepSumAgreeingIntensities::BinnedRepSumAgreeingIntensities(const BinnedRepSumAgreeingIntensities& source)
    : BinnedRepCompareFunctor(source)
  {
  }

  BinnedRepSumAgreeingIntensities& BinnedRepSumAgreeingIntensities::operator=(const BinnedRepSumAgreeingIntensities& source)
  {
		if (this != &source)
		{
    	BinnedRepCompareFunctor::operator=(source);
		}
    return *this;
  }
  
  BinnedRepSumAgreeingIntensities::~BinnedRepSumAgreeingIntensities()
  {
  }
	
	double BinnedRepSumAgreeingIntensities::operator () (const BinnedRep& a) const
	{
		return operator () (a, a);
	}

	
  double BinnedRepSumAgreeingIntensities::operator () (const BinnedRep& a, const BinnedRep& b) const
  {
    //const BinnedRep& a = csa.getBinrep();
    //const BinnedRep& b = csb.getBinrep();
    //double filterfactor = filter(csa,csb);
    //if ( filterfactor < 1e-12 ) return 0;
    double similarity = 0;
    double agreedIntensity = 0;
    double suma = 0;
    double sumb = 0;
    //go to the regions where both binreps overlap
    BinnedRep::const_iterator ait = a.begin();
    BinnedRep::const_iterator bit = b.begin();
    while (ait != a.end() && bit != b.end())
    {
      //we are at the same position (+- precision)
      if ( fabs( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  < 1e-8)
      {
        suma+= *ait;
        sumb+= *bit;
        agreedIntensity = (*ait+*bit)/2 - fabs(*ait-*bit);
        if ( agreedIntensity > 0 ) 
        {
          similarity += agreedIntensity;
        }
        ait.hop();
        bit.hop();
      }
      //ait lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  <= 1e-8)
      {
        suma+=*ait;
        ait.hop();
      }
      //bit lags
      else if ( ( ( ait.position()*a.getBinSize() + a.min() ) - ( bit.position()*b.getBinSize() + b.min() ) )  >= 1e-8)
      {
        sumb+=*bit;
        bit.hop();
      } 
    }                       
    //needed for normalization
    while ( bit != b.end() )
    {
      sumb+= *bit;
      bit.hop();
    }
    while (ait != a.end() )
    {
      suma+= *ait;
      ait.hop();
    }
    return 2/(suma+sumb)*similarity;
  }

}

