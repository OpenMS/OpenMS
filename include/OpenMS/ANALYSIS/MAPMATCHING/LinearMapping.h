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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_LINEARMAPPING_H
#define OPENMS_ANALYSIS_MAPMATCHING_LINEARMAPPING_H

#include <OpenMS/ANALYSIS/MAPMATCHING/BaseMapping.h>

namespace OpenMS
{
	/**
		 @brief This class represents a linear coordinate transformation that can
		 be applied to a pair of features.
	  
		 We assume that the coordinates of features in different maps are related
		 by an (affine) transformation. The parameters of the transformation are
		 estimated by an instance of the base class MapMatcher and it is applied
		 in class MapDewarper.
	  
	*/
	class LinearMapping
		: public BaseMapping
	{
	 public:
		/// Constructor
		LinearMapping() 
			: BaseMapping(),
				slope_(1), intercept_(0) {}
			
		/// Easy Constructor
		LinearMapping(DoubleReal slope, DoubleReal intercept) 
			: BaseMapping(),
				slope_(slope), intercept_(intercept)
		{
			this->param_.setValue("LinearMapping:slope", (double) slope);
			this->param_.setValue("LinearMapping:intercept", (double) intercept);
		}
		
		/// Destructor
	  ~LinearMapping() {}
      
    /// Copy constructor 
		LinearMapping(const LinearMapping& source) 
			: BaseMapping(source),
				slope_(source.slope_), 
				intercept_(source.intercept_)
		{
		}
			
		/// assignment operator
    LinearMapping& operator = (const LinearMapping& source);
		
    void setParam(DoubleReal sl, DoubleReal in);
			
    void setParam(const Param& pa);
		
    void apply(DPosition<1>& pos) const;
			
    void apply( DoubleReal & pos) const;
			
		inline const String getName()
		{
			return "LinearMapping";
		}
		
		/// Non-mutable access to slope
		inline const DoubleReal& getSlope() const { return slope_; }
		/// Set slope
		inline void setSlope(const DoubleReal& sl) { slope_ = sl; }
		
		/// Non-mutable access to slope
		inline const DoubleReal& getIntercept() const { return intercept_; }
		/// Set slope
		inline void setIntercept(const DoubleReal& in) { intercept_ = in; }
			
	 protected:		
		/// slope of the transform
		DoubleReal slope_;
		/// intercept 
		DoubleReal intercept_;	
	};
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_LINEARMAPPING_H

