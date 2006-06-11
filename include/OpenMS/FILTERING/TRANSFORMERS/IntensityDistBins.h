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
// $Id: IntensityDistBins.h,v 1.3 2006/03/28 12:53:13 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_INTENSITYDISTBINS_H
#define OPENMS_FILTERING_TRANSFORMERS_INTENSITYDISTBINS_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>
#include <OpenMS/COMPARISON/CLUSTERING/ClusterSpectrum.h>

#include <map>
#include <vector>
#include <string>

namespace OpenMS
{
  /**
	  @brief IntensityDistBins divides the intensity range into <i>bins</i> regions and counts
	         the peaks that fall into each bin
	  
	  \param bins number of regions
  */
  class IntensityDistBins : public FilterFunctor
  {
  public:
    /// standard constructor
    IntensityDistBins();

    /// copy constructor
    IntensityDistBins(const IntensityDistBins& source);

    /// assignment operator
    IntensityDistBins& operator=(const IntensityDistBins& source);

    /// destructor
    ~IntensityDistBins();

    static FactoryProduct* create() { return new IntensityDistBins();}

    std::vector<double> operator()(const ClusterSpectrum& spec);

    String info() const;

		static const String getName()
		{
			return "IntensityDistBins";
		}

  private:
    static const String info_;
    //all (unique) aminoacid masses
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_INTENSITYDISTBINS_H
