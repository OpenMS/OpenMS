// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H
#define OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H


#include <OpenMS/DATASTRUCTURES/IsotopeCluster.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{

	// forward declaration
	class FeatureFinder;
	
	/**@brief The purpose of this struct is to provide definitions of classes and typedefs which are used throughout all FeatureFinder classes.  */
	struct OPENMS_DLLAPI FeatureFinderDefs 
	{	
		/// Index to peak consisting of two UInts (scan index / peak index)	
		typedef IsotopeCluster::IndexPair IndexPair;
		
		/// Index to peak consisting of two UInts (scan index / peak index) with charge information
		typedef IsotopeCluster::ChargedIndexSet ChargedIndexSet;
		
		/// A set of peak indices
		typedef IsotopeCluster::IndexSet IndexSet;
		
		/// Flags that indicate if a peak is already used in a feature
		enum Flag { UNUSED, USED };

		/// Exception that is thrown if a method an invalid IndexPair is given
		class OPENMS_DLLAPI NoSuccessor :
			public Exception::BaseException
		{
			public:
				NoSuccessor(const char* file, int line, const char* function, const IndexPair& index) 
					:	BaseException(file, line, function, "NoSuccessor", "no successor/predecessor"), 
					index_(index)
				{
					what_ = String("there is no successor/predecessor for the given Index: ") + index_.first + "/" + index_.second;
					Exception::globalHandler.setMessage(what_);
				}
				virtual ~NoSuccessor() throw()
				{
				}
			protected:
				IndexPair index_;  // index without successor/predecessor
		};
	};
}
#endif // OPENMS_TRANSFORMATIONS_FEATUREFINDER_FEATUREFINDERDEFS_H
