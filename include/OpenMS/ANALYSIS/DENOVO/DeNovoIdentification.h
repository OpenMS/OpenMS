// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Sandro Andreotti, Andreas Bertsch$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H
#define OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <vector>

namespace OpenMS
{
  /**
    @brief Base class for de novo identification


    
		@ingroup Analysis_DeNovo
  */
  class OPENMS_DLLAPI DeNovoIdentification
  	: public DefaultParamHandler
  {
  	public:
		
			/** @name Constructors and destructors
			*/
			//@{
	  	/// default constructor
	  	DeNovoIdentification();
  	
			/// destructor
			virtual ~DeNovoIdentification();

  		/// copy constructor
  		DeNovoIdentification(const DeNovoIdentification& rhs);
			//@}
  		
			/// assignment operator
			DeNovoIdentification& operator = (const DeNovoIdentification& rhs);
		
      /// performs an ProteinIdentification run on a RichPeakMap
			virtual void getIdentifications(std::vector<PeptideIdentification>& ids, const RichPeakMap& exp) = 0;

			/// performs an ProteinIdentification run on a PeakSpectrum
			virtual void getIdentification(PeptideIdentification& id, const RichPeakSpectrum& spectrum) = 0;

  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_DENOVO_DENOVOIDENTIFICATION_H
