// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
 
  /**
		@brief The base class of all element group finding algorithms.

		This class defines the basic interface for all element group finding algorithms.
				
		All derived algorithms take one or several consensus maps and find corresponding features
		across the maps (or within one map). They return one consensus map containing the found consensus features.
		
		The element indices of the result consensus features are the container access indices of the
		input maps. The map indices of the result consensus features are are the indices in the input map vector.
  */
  class OPENMS_DLLAPI BaseGroupFinder
  	: public DefaultParamHandler,
			public ProgressLogger
  {
	  public:
	    /// Default constructor
	    BaseGroupFinder();
	
			/// Destructor
	    virtual ~BaseGroupFinder();
	
			/**
				@brief Run the algorithm
			
				@exception Exception::IllegalArgument is thrown if the input data is not valid.
			*/
			virtual void run(const std::vector<ConsensusMap>& input, ConsensusMap& result) = 0;
	
	    /// Register all derived classes here
	    static void registerChildren();
		
		 protected:
		 	
		 	/**
		 		@brief Checks if all file descriptions have disjoint map identifiers
		 		
		 		@exception Exception::IllegalArgument Is thrown if a file id is found twice
		  */	
		 	void checkIds_(const std::vector<ConsensusMap>& maps) const;
			
		 private:
	
	    /// Copy constructor intentionally not implemented
	    BaseGroupFinder(const BaseGroupFinder&);
			
	    /// Assignment operator intentionally not implemented
	    BaseGroupFinder & operator=(const BaseGroupFinder&);
		
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEGROUPFINDER_H
