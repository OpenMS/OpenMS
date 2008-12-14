// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKDENSITYFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKDENSITYFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

namespace OpenMS
{
  /**
  	@brief PeakDensityfilter calculates peak density

		@ingroup SpectraFilter
  */
  class OPENMS_DLLAPI PeakDensityFilter : public FilterFunctor
  {
  public:

		// @name Constructors and Destructors
		// @{
    /// default constructor
    PeakDensityFilter();

    /// copy constructor
    PeakDensityFilter(const PeakDensityFilter& source);

		/// destructor
		virtual ~PeakDensityFilter();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    PeakDensityFilter& operator = (const PeakDensityFilter& source);
		// @}

		// @name Accessors
		// @{
		///
    static FilterFunctor* create() { return new PeakDensityFilter();}

		template <typename SpectrumType> double apply(SpectrumType& spectrum)
		{
	    int nrpeaks = spectrum.size();
			double size = spectrum.getPrecursorPeak().getPosition()[0];
	    double density = nrpeaks/size;
	    return density;
		}

		///
		static const String getProductName()
		{
			return "PeakDensityFilter";
		}
		// @}

  };
}
#endif //OPENMS_FILTERING_TRANSFORMERS_PEAKDENSITYFILTER_H

