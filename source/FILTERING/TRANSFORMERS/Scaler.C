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
// $Id: Scaler.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/Scaler.h>

using namespace std;
namespace OpenMS
{

  //const String Scaler::info_ = "Scales the peaks";
  
  Scaler::Scaler()
    : PreprocessingFunctor()
  {
		name_ = PreprocessingFunctor::getName();
  }

  Scaler::Scaler(const Scaler& source)
    : PreprocessingFunctor(source)
  {
  }

  Scaler::~Scaler()
  {
  }

  Scaler& Scaler::operator=(const Scaler& source)
  {
    PreprocessingFunctor::operator=(source);
    return *this;
  }

/*
  void Scaler::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    map<double,int> peakssorted;
    int count = 0;
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end();++it )
    {
      peakssorted[it->getIntensity()] = 0;
      ++count;
    } 
    for(map<double,int>::reverse_iterator rmit = peakssorted.rbegin(); rmit != peakssorted.rend(); ++rmit)
    {
      if ( --count > 0 ) peakssorted[rmit->first] = count;
    }
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); )
    {
      double newIntensity = peakssorted[it->getIntensity()];
      if (newIntensity > 0 )
      {
        it->getIntensity() = peakssorted[it->getIntensity()];
        ++it;
      }
      else 
      {
        it = spec.getContainer().erase(it);
      }
    }                                                
  }

  String Scaler::info() const
  {
    return info_;
  }
	*/

}
