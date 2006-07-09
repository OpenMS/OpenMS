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
// $Id: BernNorm.C,v 1.5 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/BernNorm.h>

using namespace std;
namespace OpenMS
{
  BernNorm::BernNorm()
    : PreprocessingFunctor()
  {
		name_ = BernNorm::getName();
    // values from the paper
    // they should be good for GoodDiff and Complements
    // IsotopeDiffs needs lower peaks
    defaults_.setValue("C1", 28.0f);
    defaults_.setValue("C2", 400.0f);
    defaults_.setValue("threshold", 0.1f); // i.e. what is a significant peak
		param_ = defaults_;
  }

  BernNorm::BernNorm(const BernNorm& source)
    : PreprocessingFunctor(source)
  {
  }

  BernNorm::~BernNorm()
  {
  }

  BernNorm& BernNorm::operator = (const BernNorm& source)
  {
    PreprocessingFunctor::operator = (source);
    return *this;
  }

/*
  void BernNorm::operator()(MSSpectrum< DPeak<1> >& spec) const
  {
    double c1 = (double)param_.getValue("C1");
    double c2 = (double)param_.getValue("C2");
    double treshold = (double)param_.getValue("threshold");

    spec.getContainer().sortByPosition();

    // find highest peak and ranking
    double maxint = 0;
    map<double,uint> peakranks;
    for ( MSSpectrum< DPeak<1> >::iterator it = spec.begin(); it != spec.end(); ++it )
    {
      peakranks[it->getIntensity()] = 0;
      if ( it->getIntensity() > maxint )
      {
        maxint = it->getIntensity();
      }
    }
    uint rank = 0;
    for ( map<double,uint>::reverse_iterator mit = peakranks.rbegin(); mit != peakranks.rend(); ++mit )
    {
      mit->second = ++rank;
    }

    // find maxmz i.e. significant ( > threshold * maxpeak ) peak with highest m/z
    double maxmz = 0;
    for ( int i = spec.size() -1 ; i >= 0 ; --i )
    {
      if ( spec.getContainer()[i].getIntensity() > maxint * treshold )
      {
        maxmz = spec.getContainer()[i].getPosition()[0];
        break;
      }
    }

    // rank
    for ( MSSpectrum< DPeak<1> >::iterator it = spec.begin() ; it != spec.end(); )
    {
      double newint = c1-(c2/maxmz)*peakranks[it->getIntensity()];
      if ( newint < 0 )
      {
        it = spec.getContainer().erase(it);
      }
      else
      {
        it->getIntensity() = newint;
        ++it;
      }
    }
  }
*/
}
