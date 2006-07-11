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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/NLargest.h>

using namespace std;
namespace OpenMS
{

  NLargest::NLargest()
    : PreprocessingFunctor()
  {
		name_ = NLargest::getName();
    defaults_.setValue("n", 200);
		param_ = defaults_;
  }

	NLargest::NLargest(Size n)
		: PreprocessingFunctor()
	{
		name_ = NLargest::getName();
		defaults_.setValue("n", 200);
		param_ = defaults_;
		param_.setValue("n", (SignedInt)n);
	}

  NLargest::NLargest(const NLargest& source)
    : PreprocessingFunctor(source)
  {
		defaults_ = source.defaults_;
		param_ = source.param_;
		name_ = source.getName();
  }

  NLargest::~NLargest()
  {
  }

  NLargest& NLargest::operator=(const NLargest& source)
  {
    PreprocessingFunctor::operator=(source);
    return *this;
  }
/*
  void NLargest::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    // sort
    map<double,int> peakssorted;
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end();++it )
    {
      peakssorted[it->getIntensity()] = 0;
    } 
    
    // rank
    int count = (int)param_.getValue("n");
    for(map<double,int>::reverse_iterator rmit = peakssorted.rbegin(); rmit != peakssorted.rend(); ++rmit)
    {
      if ( --count > 0 ) peakssorted[rmit->first] = count;
      else break;
    }
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); )
    {
      double rank = peakssorted[it->getIntensity()];
      if (rank > 0 )
      {
        ++it;
      }
      else 
      {
        it = spec.getContainer().erase(it);
      }
    }                                                
  }
*/
}
