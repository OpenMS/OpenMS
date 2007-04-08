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
#include <OpenMS/METADATA/Identification.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief Calculates a consensus ID from several ID runs
    
    Available algorithms are:
    - Merge -- merges the runs with respect to their score. The score is not modified.
    - Ranked -- reorders the hits according to a consensus score computed from the ranks in the input runs.
                The score is normalized to the interval (0,100).
    - Average -- reorders the hits according to the average score of the input runs.
		
		The following parameters can be given:
		- Algorithm -- see above.
		- ConsideredHits -- the number of top hits to use as input.
		- NumberOfRuns -- the number of runs used as input. This information is used in 'Ranked' and 'Average' to
		                  compute the new scores.
		- InverseOrder -- indicates that the lower scores are better scores.
		- MinOutputScore -- the minimum score a hit has to have to be reported as output.
  */
  class ConsensusID
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	ConsensusID();
  		
  		///Calculates the consensus ID for a Feature
  		void apply(std::vector<Identification>& ids) throw (Exception::InvalidValue);
  		
  	private:
  		///Hidden and not implemented copy constructor
  		ConsensusID(const ConsensusID&);
  		
			///Hidden and not implemented assignment operator
			ConsensusID& operator = (const ConsensusID&);
			
			/// Merge algorithm
			void merge_(std::vector<Identification>& ids);
			/// Ranked algorithm
			void ranked_(std::vector<Identification>& ids);
			/// Average score algorithm
			void average_(std::vector<Identification>& ids);
			
			/// This flag indicates that a bigger score is better
			bool inverse_order_;
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSID_H
