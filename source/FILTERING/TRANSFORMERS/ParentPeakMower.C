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
// $Id: ParentPeakMower.C,v 1.7 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>

using namespace std;

namespace OpenMS
{

  //const String ParentPeakMower::info_ = "scales down large ( greater [x]*mean(Peaks) ) Peaks which have the same (+-[windowsize] m/z as the precursorpeak";
  
  ParentPeakMower::ParentPeakMower()
    : PreprocessingFunctor()
  {
		name_ = ParentPeakMower::getName();
    defaults_.setValue("windowsize", 35);
		param_ = defaults_;
  }

  ParentPeakMower::ParentPeakMower(const ParentPeakMower& source)
    : PreprocessingFunctor(source)
  {
		name_ = source.getName();
  }

  ParentPeakMower::~ParentPeakMower()
  {
  }

  ParentPeakMower& ParentPeakMower::operator = (const ParentPeakMower& source)
  {
    PreprocessingFunctor::operator = (source);
    return *this;
  }

/*
  String ParentPeakMower::info() const
  {
    return info_;
  }

  void ParentPeakMower::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    double window = (double)param_.getValue("windowsize");
    double mean = 0;
    
    spec.getContainer().sortByPosition();
    
    //calculate mean
    for (MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); ++it)
    {
      mean += it->getIntensity();
    }
    mean /= spec.size();
    
    // assumed position of precursorpeak
    double pppos = spec.getPrecursorPeak().getPosition()[0] / spec.getPrecursorPeak().getCharge();
    MSSpectrum< DPeak<1> >::iterator it = spec.end();
    if ( it == spec.begin() ) return;
    do
    {
      --it;
      if (it->getPosition()[0] <= pppos + window && it->getPosition()[0] >= pppos - window)
      {
        if( it->getIntensity() > mean ) it->setIntensity(mean);
      }
      else if (it->getPosition()[0] < pppos - window)
      {
        break;
      }
    }while (it != spec.begin());
  }
*/
}
