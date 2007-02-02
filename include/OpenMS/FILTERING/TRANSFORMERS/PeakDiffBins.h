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
#ifndef OPENMS_FILTERING_TRANSFORMERS_PEAKDIFFBINS_H
#define OPENMS_FILTERING_TRANSFORMERS_PEAKDIFFBINS_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <vector>

namespace OpenMS
{
	class ClusterSpectrum;
	
  /**
  	@brief PeakDiffBins calculates all differences between pairs and returns the fraction of the total intensity in the specified regions

		@ingroup SpectraFilter
  */
  class PeakDiffBins : public FilterFunctor
  {
  public:

		// @name
		// @{
    /// default constructor
    PeakDiffBins();

    /// copy constructor
    PeakDiffBins(const PeakDiffBins& source);

		/// destructor
		virtual ~PeakDiffBins();
		// @}

		// @name Operators
		// @{
    /// assignment operator
    PeakDiffBins& operator = (const PeakDiffBins& source);
		// @}

		// @name Accessors
		// @{
		///
		static FilterFunctor* create() { return new PeakDiffBins(); }

		///
    std::vector<double> operator () (const ClusterSpectrum& spec);

		///
		static const String getProductName()
		{
			return "PeakDiffBins";
		}

    /// change layout of bins/regions
    void setmask(std::vector<double>& newmask);
		// @}

  private:
	
    /**
    current layout of bins/regions <br>
    standard is 1-187, size 1 <br>
    */
    std::map<double, int> mask_;
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_PEAKDIFFBINS_H
