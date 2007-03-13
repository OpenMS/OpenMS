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

namespace OpenMS
{
	class Feature;
	
  /**
    @brief Calculates a consensus ID from several ID runs
    
		
    
  */
  class ConsensusID
  	: public DefaultParamHandler
  {
  	public: 
	  	///Default constructor
	  	ConsensusID();
  		
  		///Calculates the consensus ID for a Feature
  		void apply(Feature& feature) throw (Exception::InvalidValue);
  		
  	private:
  		///Hidden and not implemented copy constructor
  		ConsensusID(const ConsensusID&);
  		
			///Hidden and not implemented assignment operator
			ConsensusID& operator = (const ConsensusID&);
			
			/// Merge algorithm
			void merge_(Feature& feature);
			/// Ranked algorithm
			void ranked_(Feature& feature);
			/// Average score algorithm
			void average_(Feature& feature);
  };
 
} // namespace OpenMS

#endif // OPENMS_ANALYSIS_ID_CONSENSUSID_H
