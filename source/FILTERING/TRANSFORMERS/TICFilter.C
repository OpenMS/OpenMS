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
// $Id: TICFilter.C,v 1.7 2006/06/09 13:52:51 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#include <OpenMS/FILTERING/TRANSFORMERS/TICFilter.h>
#include <stdexcept>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cmath>
#include <ctime>

using namespace std;

namespace OpenMS
{
  TICFilter::TICFilter()
    :FilterFunctor()
  {
		name_ = TICFilter::getName();
    defaults_.setValue("window", 5);
		param_ = defaults_;
  }

  TICFilter::TICFilter(const TICFilter& source)
    :FilterFunctor(source)
  {
  }

  TICFilter& TICFilter::operator=(const TICFilter& source)
  {
    FilterFunctor::operator=(source);
    return *this;
  }
  
  TICFilter::~TICFilter()
  {
  }

	/*
  vector<double> TICFilter::operator()(const ClusterSpectrum& cspec)
  {
    vector<double> result;
    double TIC = 0;
    double PTIC = 0;
    double window = (double)param_.getValue("window");
    double parentpeak = cspec.getParentMass()/ cspec.getParentionCharge();
    for (MSSpectrum< DPeak<1> >::const_iterator it = cspec.getSpec().begin(); it != cspec.getSpec().end();++it )
    {
      TIC += it->getIntensity();
      if ( fabs( it->getPosition()[0] - parentpeak )  <  window )
      {
        PTIC += it->getIntensity();
      }
    }
    result.push_back(TIC);
    result.push_back(PTIC);
    return result;
  }
	*/
}
