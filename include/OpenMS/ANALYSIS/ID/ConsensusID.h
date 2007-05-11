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
    
    Available algorithms are:
    - Merge -- merges the runs with respect to their score. The score is not modified.
    					 Make sure to use PeptideIdentifications with the same score type only!
    - Ranked -- reorders the hits according to a consensus score computed from the ranks in the input runs.
                The score is normalized to the interval (0,100).
                The PeptideIdentifications do not need to have the same score type.
    - Average -- reorders the hits according to the average score of the input runs.
    						 Make sure to use PeptideIdentifications with the same score type only!
		
		The following parameters can be given:
		- Algorithm -- see above.
		- ConsideredHits -- the number of top hits to use as input.
		- NumberOfRuns -- the number of runs used as input. This information is used in 'Ranked' and 'Average' to
		                  compute the new scores. If not given, the number of input identifications is taken. 
  */
  class ConsensusID
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	ConsensusID();
  		
  		/**
  			@brief Calculates the consensus ID for a set of PeptideIdentification
  			
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
