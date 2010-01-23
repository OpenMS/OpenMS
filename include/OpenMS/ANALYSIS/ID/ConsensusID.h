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
// $Authors: Andreas Bertsch, Marc Sturm, Sven Nahnsen $
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
    
    This class combines several ID runs using one of several, available algorithms.
		
		@htmlinclude OpenMS_ConsensusID.parameters

		@ingroup Analysis_ID
  */
  class OPENMS_DLLAPI ConsensusID
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	ConsensusID();
  		
  		/**
  			@brief Calculates the consensus ID for a set of PeptideIdentification instances of the same spectrum
  			
  			@note Make sure that the score orientation (PeptideIdentification::isHigherScoreBetter())is set properly!
  		*/
  		void apply(std::vector<PeptideIdentification>& ids);
  		
  	private:
  		///Not implemented
  		ConsensusID(const ConsensusID&);
  		
			///Not implemented
			ConsensusID& operator = (const ConsensusID&);
			
			/// Merge algorithm
			void merge_(std::vector<PeptideIdentification>& ids);
			
			/// Ranked algorithm
			void ranked_(std::vector<PeptideIdentification>& ids);
			
			/// Average score algorithm
			void average_(std::vector<PeptideIdentification>& ids);
      
			/// Probability based algorithm
			void probability_(std::vector<PeptideIdentification>& ids);

			/// Majority vote
			void majority_(std::vector<PeptideIdentification>& ids);
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSID_H
