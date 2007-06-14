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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_ID_CONSENSUSID_H
#define OPENMS_ANALYSIS_ID_CONSENSUSID_H

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates a consensus ID from several ID runs
    
    This class combines several ID runs using one of the following algorithms:
    - Marge
    - Ranked
    - Average

		@ingroup Analysis_ID

  	@ref ConsensusID_Parameters are explained on a separate page.
  */
  class ConsensusID
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	ConsensusID();
  		
  		/**
  			@brief Calculates the consensus ID for a set of PeptideIdentification instances of the same spectrum
  			
  			@note Make sure that the score orientation (PeptideIdentification::isHigherScoreBetter())is set properly!
  		*/
  		void apply(std::vector<PeptideIdentification>& ids) throw (Exception::InvalidValue);
  		
  	private:
  		///Hidden and not implemented copy constructor
  		ConsensusID(const ConsensusID&);
  		
			///Hidden and not implemented assignment operator
			ConsensusID& operator = (const ConsensusID&);
			
			/// Merge algorithm
			void merge_(std::vector<PeptideIdentification>& ids);
			/// Ranked algorithm
			void ranked_(std::vector<PeptideIdentification>& ids);
			/// Average score algorithm
			void average_(std::vector<PeptideIdentification>& ids);
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSID_H
