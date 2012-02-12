// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Sandro Andreotti $
// $Authors: Sandro Andreotti, Andreas Bertsch$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_DENOVOPOSTSCORING_H
#define OPENMS_ANALYSIS_DENOVO_DENOVOPOSTSCORING_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for ion scoring implementation for de novo algorithms


    
		@ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI DeNovoPostScoring
  	: public DefaultParamHandler
  {
  	public:
		
			/** @name Constructors and destructors
			*/
			//@{
	  	/// default constructor
	  	DeNovoPostScoring();
  	
			/// destructor
			virtual ~DeNovoPostScoring();

  		/// copy constructor
  		DeNovoPostScoring(const DeNovoPostScoring& rhs);
			//@}
  		
			/// assignment operator
			DeNovoPostScoring& operator = (const DeNovoPostScoring& rhs);

			virtual void apply(std::vector<PeptideIdentification>& identifications, const RichPeakMap& exp) = 0;

			virtual void apply(PeptideIdentification& identification, const RichPeakSpectrum& spec) = 0;

  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DENOVO_DENOVOPOSTSCORING_H
