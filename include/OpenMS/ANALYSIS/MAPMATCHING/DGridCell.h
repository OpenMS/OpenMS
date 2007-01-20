// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2: 
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

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DGRIDCELL_H
#define OPENMS_ANALYSIS_MAPMATCHING_DGRIDCELL_H

#include <OpenMS/KERNEL/KernelTraits.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DBaseMapping.h>

#include <vector>

namespace OpenMS
{
	
	/** @brief D-dimensional grid cell over a map.
	 */
	template <Size D, typename Traits = KernelTraits>
	class DGridCell : public DRange<D, Traits>
	{
		public:
		enum { DIMENSION = D };
		typedef Traits TraitsType;
		typedef typename Traits::CoordinateType  CoordinateType;
		typedef DPosition<D,Traits> PositionType;
		typedef DBaseMapping<1,Traits> MappingType;
		typedef std::vector<MappingType*> MappingVector;
		
		/** @name Constructors and Destructor
		 */
		//@{
		/// Default constructor		
		DGridCell() : 
			DRange<D,TraitsType>(),
			mappings_()			
			{	}

		/// "Convenience contructor"
		DGridCell(CoordinateType i, CoordinateType j, 
							CoordinateType k, CoordinateType l) : 
			DRange<D,TraitsType>(i, j, k, l),	
			mappings_()
			{	}  

		/// Copy constructor
		DGridCell(const DGridCell& gc) 
			: DRange<D,TraitsType>(gc)
		{
			setMappings(gc.mappings_);
		}		

		/// Destructor
		virtual ~DGridCell() 
		{			
		}
		//@}
		
		/// assignment operator
    DGridCell& operator = (const DGridCell& rhs)
		{
			if (&rhs==this) return *this;			
			DRange<D,TraitsType>::operator = (rhs);
			mappings_ = rhs.mappings_;
			return *this;
		}
		
		/// Test for equality
   	bool operator == (const DGridCell& rhs) const
		{
			return (DRange<D,TraitsType>::operator == (rhs) 
							&& mappings_ == rhs.mappings_);			
		}
		
		/// Test for inequality
   	bool operator != (const DGridCell& rhs) const
		{
			return ! (*this == rhs);
		}
		
		/** @name Accesssor methods
		*/
		//@{	
		/// Set transform 
    	void setMappings(MappingVector const & m) { mappings_ = m; }
		/// Get transform 
   		MappingVector& getMappings() { return mappings_; }
    	/// Get transform (non-mutable)
    	const MappingVector& getMappings() const { return mappings_; }
		 //@} 
				 	
		protected:
		/// We estimate a different mapping for each coordinate and
		/// therefore store a vector of transformations
		MappingVector mappings_;
		
	}; // end of class DGridCell

} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DGRIDCELL_H
