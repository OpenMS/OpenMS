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

#ifndef OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H
#define OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H

#include <OpenMS/COMPARISON/SPECTRA/BinnedRepCompareFunctor.h>

namespace OpenMS
{

	/** 
		@brief calculates the Shared Peak Count for two binned spectra

		the shared peak count was defined in ????

		@ingroup SpectraComparison
	*/

  class BinnedRepSharedPeakCount
    : public BinnedRepCompareFunctor
  {
  	public:

			// @name Constructors and Destructors
			// @{
			/// default constructor
    	BinnedRepSharedPeakCount();
			
			/// copy constructor
    	BinnedRepSharedPeakCount(const BinnedRepSharedPeakCount& source);

			/// destructor
    	virtual ~BinnedRepSharedPeakCount();
			// @}

			// @name Operators
			// @{
			/// assignment operator
    	BinnedRepSharedPeakCount& operator=(const BinnedRepSharedPeakCount& source);
	
			/// 
			double operator () (const BinnedRep& csa, const BinnedRep& csb) const;

			double operator () (const BinnedRep& a) const;
			// @}
			
			// @name Accessors
			// @{
			///
    	static BinnedRepCompareFunctor* create() {return new BinnedRepSharedPeakCount();}

			///
			static const String getProductName()
			{
				return "BinnedRepSharedPeakCount";
			}
			// @}
  };

}

#endif //OPENMS_COMPARISON_SPECTRA_BINNEDREPSHAREDPEAKCOUNT_H

