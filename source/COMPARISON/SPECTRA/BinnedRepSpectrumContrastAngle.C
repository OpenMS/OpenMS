//   -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/COMPARISON/SPECTRA/BinnedRepSpectrumContrastAngle.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <cmath>

using namespace std;

namespace OpenMS
{
  BinnedRepSpectrumContrastAngle::BinnedRepSpectrumContrastAngle()
		: BinnedRepCompareFunctor()
  {
		setName(BinnedRepSpectrumContrastAngle::getProductName());
  }

  BinnedRepSpectrumContrastAngle::BinnedRepSpectrumContrastAngle(const BinnedRepSpectrumContrastAngle& source)
    : BinnedRepCompareFunctor(source)
  {
  }

  BinnedRepSpectrumContrastAngle& BinnedRepSpectrumContrastAngle::operator= (const BinnedRepSpectrumContrastAngle& source)
  {
		if (this != &source)
		{
			BinnedRepCompareFunctor::operator=(source);
		}
    return *this;
  }
  
  BinnedRepSpectrumContrastAngle::~BinnedRepSpectrumContrastAngle()
  {
  }

	double BinnedRepSpectrumContrastAngle::operator () (const BinnedRep& a) const
	{
		return operator () (a, a);
	}
	
  double BinnedRepSpectrumContrastAngle::operator()(const BinnedRep& csa, const BinnedRep& csb) const
  {
    //all binreps do this, check if sizes are similar, if bindimensions are compatible...
    double similarity(0), sum_a(0), sum_b(0);
    
    BinnedRep a = csa;
		a.normalize();
    BinnedRep b = csb;
		b.normalize();
    
		BinnedRep::const_iterator ait = a.begin();
    BinnedRep::const_iterator bit = b.begin();
    
    //go to overlapping region of both binreps
    while ( (a.min() - (b.min() + bit.position() * b.getBinSize())) > 1e-12 && bit != b.end() ) // == a.min() > b.min() + bpos*b.getBinSize()
    {
      bit.hop();
    }
    while ( (b.min() -  (a.min() + ait.position() * a.getBinSize())) > 1e-12 && ait != a.end() )
    {
      ait.hop();
    }

    //bins at same positions with high heights => score
    while ( bit != b.end() && ait != a.end())
    {
      if ( fabs ( (a.min() + ait.position() * a.getBinSize()) - ( b.min() + bit.position() * b.getBinSize())) < 1e-12 ) 
      {
        similarity += *bit * *ait;
				sum_a += *ait;
				sum_b += *bit;
        ait.hop();
        bit.hop();
      }
      else if ( (a.min() + ait.position() * a.getBinSize()) - ( b.min() + bit.position() * b.getBinSize() ) > 0 ) 
      {
        bit.hop();
      }
      else if ( (a.min() + ait.position() * a.getBinSize()) - ( b.min() + bit.position() * b.getBinSize() ) < 0 ) 
      {
        ait.hop();
      }
      else 
      {
				/// @todo another Exeption?! (Andreas)
        //throw CanNotHappen(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
    }
		//similarity /= sqrt(sum_a * sum_b);
    return similarity;
  }

}
