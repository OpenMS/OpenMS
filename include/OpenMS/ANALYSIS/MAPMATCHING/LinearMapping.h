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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_LINEARMAPPING_H
#define OPENMS_ANALYSIS_MAPMATCHING_LINEARMAPPING_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>

namespace OpenMS
{
	// forward declaration to eliminate need for including the class header file
	class ConsensusMap;

	/**
		 @brief This class represents a linear coordinate transformation.
		 
		 @todo Replace by generic mapping that stores name and parameters. GenericMapping is used to create the transformation function (Clemens) 
	*/
	class LinearMapping
	{
	 public:
		/// Constructor
		LinearMapping();
		/// Destructor
	  ~LinearMapping();
      
    /// Copy constructor 
		LinearMapping(const LinearMapping& source);
		/// assignment operator
    LinearMapping& operator = (const LinearMapping& source);
	  /// Equality operator
    bool operator==(const LinearMapping& rhs) const;
	  /// Equality operator
    bool operator!=(const LinearMapping& rhs) const;

		/// Apply transformation to coordinate
    void apply(DPosition<1>& pos) const;

		/// Apply transformation to coordinate
    void apply( DoubleReal & pos) const;

		/// Apply transformation to retention time of each consensus feature
		template< typename ContainerT >
    void applyRT( ContainerT & container ) const
		{
			for ( typename ContainerT::iterator iter = container.begin();
						iter != container.end(); ++iter )
			{ 
				DoubleReal pos = iter->getRT();
				apply(pos);
				iter->setRT(pos);
			}
			return;
		}
			
		/// Apply transformation to retention time of each consensus feature
		template< typename ContainerT >
    void applyMZ( ContainerT & container ) const
		{
			for ( typename ContainerT::iterator iter = container.begin();
						iter != container.end(); ++iter )
			{ 
				DoubleReal pos = iter->getMZ();
				apply(pos);
				iter->setMZ(pos);
			}
			return;
		}
			
		/// Non-mutable access to slope
		inline DoubleReal getSlope() const
		{ 
			return slope_;
		}
		/// Set slope
		inline void setSlope(DoubleReal sl) 
		{ 
		  slope_ = sl; 
		}
		
		/// Non-mutable access to intercept
		inline DoubleReal getIntercept() const 
		{ 
			return intercept_;
		}
		/// Set intercept
		inline void setIntercept(DoubleReal in) 
		{ 
		  intercept_ = in; 
		}

	 protected:
		/// slope of the transform
		DoubleReal slope_;
		/// intercept 
		DoubleReal intercept_;	
	};
	
	/// output to a stream, for debugging
	std::ostream& operator<< (std::ostream& os, LinearMapping const & rhs);

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_LINEARMAPPING_H

