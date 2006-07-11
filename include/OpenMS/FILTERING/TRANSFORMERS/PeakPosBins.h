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
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKPOSBINS_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKPOSBINS_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <map>
#include <vector>
#include <string>

namespace OpenMS
{
  /**
  	@brief PeakPosBins sums the intensity in <i>bins</i> regions<br>
  
  	\param bins number of regions
  */
  class PeakPosBins : public FilterFunctor
  {
  public:
    /// standard constructor
    PeakPosBins();

    /// copy constructor
    PeakPosBins(const PeakPosBins& source);

    /// assignment operator
    PeakPosBins& operator=(const PeakPosBins& source);

    /// destructor
    ~PeakPosBins();

    static FactoryProduct* create() { return new PeakPosBins();}

    std::vector<double> operator()(const ClusterSpectrum& spec);

		static const String getName()
		{
			return "PeakPosBins";
		}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_PEAKPOSBINS_H
