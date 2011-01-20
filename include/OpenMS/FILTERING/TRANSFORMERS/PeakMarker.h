// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H

#include <map>

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

namespace OpenMS
{
	/**
  	@brief PeakMarker marks peaks that seem to fulfill some criterion

  */
  class OPENMS_DLLAPI PeakMarker
    : public DefaultParamHandler
  {
  public:

    /// default constructor
    PeakMarker() ;

    /// copy constructor
    PeakMarker(const PeakMarker& source);

    /// destructor
    virtual ~PeakMarker();

    /// assignment operator
    PeakMarker& operator = (const PeakMarker& source);

		/// method to mark peaks
		template <typename SpectrumType> void apply(std::map<double, bool>& /* marked */, SpectrumType& /* spectrum */) {}

		///
		static const String getProductName()
		{
			return "PeakMarker";
		}
  };

}
#endif // OPENMS_FILTERING_TRANSFORMERS_PEAKMARKER_H
