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
#include <OpenMS/FILTERING/TRANSFORMERS/PeakPosBins.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <algorithm>

using namespace std;

namespace OpenMS
{
  PeakPosBins::PeakPosBins()
    :FilterFunctor() 
  {
		setName(PeakPosBins::getProductName());
		defaults_.setValue("bins", 10);
		defaultsToParam_();
  }

  PeakPosBins::PeakPosBins(const PeakPosBins& source)
    : FilterFunctor(source)
  {
  }
    
  PeakPosBins& PeakPosBins::operator = (const PeakPosBins& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator=(source);
		}
    return *this;
  }

  PeakPosBins::~PeakPosBins()
  {
  }

  vector<double> PeakPosBins::operator() ( const ClusterSpectrum& cspec)
  {
    vector<double> result = vector<double>((unsigned int)param_.getValue("bins"));
		double min = std::numeric_limits<double>::max();
    double max = 0.0;
    for ( MSSpectrum< DPeak<1> >::const_iterator dit = cspec.getSpec().begin(); dit != cspec.getSpec().end(); ++dit )
    {
      min = std::min(min, dit->getPosition()[0]);
      max = std::max(max, dit->getPosition()[0]);
    }
    for (MSSpectrum< DPeak<1> >::const_iterator dit = cspec.getSpec().begin(); dit != cspec.getSpec().end(); ++dit)
    {
      uint pos = (uint)(((dit->getPosition()[0]-min)/(max-min)) * result.size());
      // todo: is this necessary?
      if (pos == result.size()) --pos;
      result[pos] += dit->getIntensity();
    } 
    return result;
  }
}  
