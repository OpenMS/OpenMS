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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_BASEMAPPING_H
#define OPENMS_ANALYSIS_MAPMATCHING_BASEMAPPING_H

#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/FORMAT/Param.h>

namespace OpenMS
{
	/**
		 @brief This base class represents an coordinate transformation
		 that can be applied to a 2D-Position or a DoubleReal. 
	  
	*/
	class BaseMapping
	{
	 public:
		
		/// Constructor
		BaseMapping() {}
		
		/// Destructor
		virtual ~BaseMapping() {}
			
		/// Copy constructor 
		BaseMapping(const BaseMapping& source)
		{
			setParam(source.param_);
		}
			
		/// Assignment operator
		BaseMapping& operator = (const BaseMapping& rhs)
		{
			if (this==&rhs) return *this;
			
			param_ = rhs.param_;
			return *this;
		}		
		
		/// Set parameters
    virtual void setParam(const Param& p);
		/// Get parameters
    virtual const Param& getParam() const; 
    						
		/// Apply the transform to a feature
    virtual void apply(DPosition<1>& ) const = 0;

		/// Apply the transformation
    virtual void apply( DoubleReal& pos) const = 0;
	
		/// Return the name of this transformation
		virtual const String getName() = 0;
					
	 protected:		
		/// Parameters defining the transformation
		Param param_;				
	}; // end of BaseMapping
					
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DBASEMAPPING_H

