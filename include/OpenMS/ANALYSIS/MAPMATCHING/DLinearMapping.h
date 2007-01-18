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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DLINEARMAPPING_H
#define OPENMS_ANALYSIS_MAPMATCHING_DLINEARMAPPING_H

#include <OpenMS/KERNEL/KernelTraits.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

namespace OpenMS
{
	class Param;
	
	/**
		 @brief This class represents a linear coordinate transformation that can
		 be applied to a pair of features.
	  
		 We assume that the coordinates of features in different maps are related
		 by an (affine) transformation. The parameters of the transformation are
		 estimated by an instance of the base class MapMatcher and it is applied
		 in class MapDewarper.
	  
	*/
	template <UnsignedInt D, typename Traits = KernelTraits>
	class DLinearMapping
		: public DBaseMapping<D, Traits>
	{
	 public:
		typedef typename Traits::RealType SlopeType;
		typedef typename Traits::RealType InterceptType;
				
		/// Constructor
		DLinearMapping() 
			: DBaseMapping<D, Traits>(),
				slope_(0), intercept_(0) {}
			
		/// Easy Constructor
		DLinearMapping(SlopeType sl, InterceptType in) 
			: DBaseMapping<D, Traits>(),
				slope_(sl), intercept_(in)
		{
			this->param_.setValue("DLinearMapping:slope", (float) sl);
			this->param_.setValue("DLinearMapping:intercept", (float) in);
		}
		
		/// Destructor
		virtual ~DLinearMapping() {}
			
		/// Copy constructor 
		DLinearMapping(const DLinearMapping& source) 
			: DBaseMapping<D, Traits>(source),
				slope_(source.slope_), 
				intercept_(source.intercept_)
		{
		}
			
		/// assignment operator
		DLinearMapping& operator = (const DLinearMapping& source)
		{
			if (this==&source) return *this;
			
			DBaseMapping<D, Traits>::operator = (source);
			slope_     = source.slope_;
			intercept_ = source.intercept_;
			return *this;
		}		
		
		/// equality operator
		bool operator == (const DLinearMapping& rhs)
		{
			return (DBaseMapping<D, Traits>::operator == (rhs) 
							&& slope_     == rhs.slope_ 
							&& intercept_ == rhs.intercept_);
			
		}	
		
		/// equality operator
		bool operator != (const DLinearMapping& rhs)
		{
			return !(*this == rhs);
		}	
			
		void setParam(SlopeType sl, InterceptType in)
		{
			slope_ = sl;
			this->param_.setValue("DLinearMapping:slope", (double) sl);
			intercept_ = in;
			this->param_.setValue("DLinearMapping:intercept", (double) in);
		}
			
		void setParam(const Param& pa)
		{
			slope_       = pa.getValue("DLinearMapping:slope");
			intercept_   = pa.getValue("DLinearMapping:intercept");
			
			this->param_ = pa;
		}
		
		// ok for D==1 but otherwise ????
		// Who cares about D>1 anyway ? :-)
		void apply(DPosition<D>& pos) const
		{
			for (UnsignedInt i=0; i < D; ++i)
			{
				pos[i] = intercept_ + slope_ * pos[i];
			}
		}
			
		void apply( typename Traits::RealType & pos) const
		{
			pos *= slope_;
			pos += intercept_;
		}
			
		const String getName()
		{
			return "DLinearMapping";
		}
		
		/// Non-mutable access to slope
		const SlopeType& getSlope() const { return slope_; }
		/// Mutable access to slope
		SlopeType& getSlope() { return slope_; }
		/// Set slope
		void setSlope(const SlopeType sl) { slope_ = sl; }
		
		/// Non-mutable access to slope
		const InterceptType& getIntercept() const { return intercept_; }
		/// Mutable access to slope
		InterceptType& getIntercept() { return intercept_; }
		/// Set slope
		void setIntercept(const InterceptType in) { intercept_ = in; }
			
	 protected:		
		/// slope of the transform
		SlopeType slope_;
		/// intercept 
		InterceptType intercept_;
			
	};
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DLINEARMAPPING_H

