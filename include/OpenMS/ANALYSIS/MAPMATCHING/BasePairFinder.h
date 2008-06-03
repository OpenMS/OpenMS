// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/CONCEPT/FactoryProduct.h>

#include <utility>
#include <fstream>

namespace OpenMS
{
 
  /**
		@brief The base class of all element pair finding algorithms.
			
		This class defines the basic interface for all element pair finding
		algorithms. It works on two consensus maps.
  */
  class BasePairFinder
  	: public FactoryProduct
  {
	  public:
	    /// Default constructor
	    BasePairFinder();
	
			/// Destructor
	    virtual ~BasePairFinder();
	
			/// Run the algorithm
			virtual void run(const std::vector<ConsensusMap>&, ConsensusMap&)
			{
			}
	
	    /// Register all derived classes here
	    static void registerChildren();
	
			
		 private:
	
	    /// Copy constructor intentionally not implemented
	    BasePairFinder(const BasePairFinder&);
			
	    /// Assignment operator intentionally not implemented
	    BasePairFinder & operator=(const BasePairFinder&);
		
	};

} // namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHING_BASEPAIRFINDER_H
