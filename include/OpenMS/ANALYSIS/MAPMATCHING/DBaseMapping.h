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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DBASEMAPPING_H
#define OPENMS_ANALYSIS_MAPMATCHING_DBASEMAPPING_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/Param.h>

namespace OpenMS
{
	/**
		 @brief This base class represents an coordinate transformation
		 that can be applied to a DPosition. 
	  
	*/
	template <UnsignedInt D, typename Traits = KernelTraits>		
	class DBaseMapping
	{
	 public:
		
		/// Constructor
		DBaseMapping() {}
		
		/// Destructor
		virtual ~DBaseMapping() {}
			
		/// Copy constructor 
		DBaseMapping(const DBaseMapping& source)
		{
			setParam(source.param_);
		}
			
		/// Assignment operator
		DBaseMapping& operator = (const DBaseMapping& rhs)
		{
			if (this==&rhs) return *this;
			
			param_ = rhs.param_;
			return *this;
		}		
		
		/// Equality operator
		bool operator == (const DBaseMapping& rhs)
		{
			return (param_ == rhs.param_);
		}	
		
		/// Inequality operator
		bool operator != (const DBaseMapping& rhs)
		{
			return !(param_ == rhs.param_);
		}	
		
		/// Set parameters
		virtual void setParam(const Param& p) { param_ = p; } 
		/// Get parameters
		virtual const Param& getParam() const { return param_; }
						
		/// Apply the transform to a feature
		virtual void apply(DPosition<D>& ) const = 0;

		/// Apply the transformation
		virtual void apply( typename Traits::RealType & pos) const = 0;
	
		/// Return the name of this transformation
		virtual const String getName() = 0;
					
	 protected:		
		/// Parameters defining the transformation
		Param param_;				
	}; // end of DBaseMapping
	
	///Print the contents to a stream.
	template <Size D, typename Traits>
	std::ostream& operator << (std::ostream& os, const DBaseMapping<D, Traits>& mapping)
	{
		os << mapping.getParam();
		return os;
	}
					
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DBASEMAPPING_H

