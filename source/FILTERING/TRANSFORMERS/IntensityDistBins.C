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
#include <OpenMS/FILTERING/TRANSFORMERS/IntensityDistBins.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <cstdlib>

using namespace std;

namespace OpenMS
{
  IntensityDistBins::IntensityDistBins()
    : FilterFunctor()
  {
		setName(IntensityDistBins::getProductName());
    defaults_.setValue("bins", 10, "The number of bins.");
		defaultsToParam_();
  }
  
  IntensityDistBins::IntensityDistBins(const IntensityDistBins& source)
    : FilterFunctor(source)
  {
  }

  IntensityDistBins& IntensityDistBins::operator = (const IntensityDistBins& source)
  {
		if (this != &source)
		{
    	FilterFunctor::operator = (source);
		}
    return *this;
  }
  
  IntensityDistBins::~IntensityDistBins()
  {
  }

  vector<double> IntensityDistBins::operator() ( const ClusterSpectrum& cspec)
  {
    uint bins = (unsigned int)param_.getValue("bins");
    Real min = 1e258;
    Real max = 0;
    uint count = 0;
    for (MSSpectrum<Peak1D >::const_iterator dit = cspec.getSpec().begin(); dit != cspec.getSpec().end(); ++dit)
    {
      min = std::min(min,dit->getIntensity());
      max = std::max(max,dit->getIntensity());
      ++count;
    }

    vector<double> result = vector<double>(bins);
    
    for (MSSpectrum<Peak1D >::const_iterator dit = cspec.getSpec().begin(); dit != cspec.getSpec().end(); ++dit)
    {
      //uint bin = ( ( max - log(dit->getIntensity()) )/ (max-min) )  * bins ;
      uint bin = (uint)((max - dit->getIntensity()) / (max-min)) * bins ;
      if (bin == result.size())
			{
				--bin;
			}
      result.at(bin)++;
    }

    for (vector<double>::iterator vit = result.begin(); vit != result.end(); ++vit )
    {
      *vit /= count;
    }

    return result;
  }
}
