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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_DGRID_H
#define OPENMS_ANALYSIS_MAPMATCHING_DGRID_H

#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class represents a D-dimensional grid over a LC/MS map.	
	
	*/
	template<Size D> 
	class DGrid : public std::vector<DGridCell<D> >
	{
	 public:
			
		/**	
			 @name Type definitions
		*/
		//@{
		typedef DGridCell<D> CellType;
		typedef std::vector<CellType> Base;
		typedef std::vector<CellType> ContainerType;
		typedef typename ContainerType::iterator Iterator;
		typedef typename ContainerType::const_iterator ConstIterator;
		typedef typename ContainerType::reverse_iterator ReverseIterator;
		typedef typename ContainerType::const_reverse_iterator ConstReverseIterator;
		typedef CellType& Reference;
		typedef const CellType& ConstReference;
			
		// STL compatibility
		typedef CellType value_type;
		typedef CellType* pointer;
		typedef const CellType* const_pointer;
		typedef Reference reference;
		typedef ConstReference const_reference;
		typedef typename ContainerType::size_type size_type;
		typedef typename ContainerType::difference_type difference_type;
		typedef Iterator iterator;
		typedef ConstIterator const_iterator;
		typedef ReverseIterator reverse_iterator;
		typedef ConstReverseIterator const_reverse_iterator;
		//@}
			
		/**	
			 @name Constructors and Destructor
		*/
		//@{
		/// Default constructor
		DGrid() {}
		/// Copy constructor
		DGrid(const DGrid& grid)
			: Base(grid)
		{}
		/// Destructor
		virtual ~DGrid() {}
		//@}
			
		/// Assignment operator
		DGrid& operator = (const DGrid& rhs)
		{
			if (&rhs==this) return *this;
				
			Base::operator=(rhs);
								
			return *this;
		}

		/// Equality operator
		bool operator == (const DGrid& rhs) const
		{
			return	std::operator==(*this, rhs)	;				
		}
			
		/// Equality operator
		bool operator != (const DGrid& rhs) const
		{
			return !(operator==(rhs));
		}
						
	}; // end of class DGrid
  
  ///Print the contents to a stream.
  template < Size D >
  std::ostream& operator << (std::ostream& os, const DGrid<D>& grid)
  {
    os << "---------- Grid BEGIN -----------------\n";
    for (typename DGrid<D>::const_iterator it = grid.begin(); it != grid.end(); ++it)
    {
      os  << "GridCell: " <<  *it << std::endl;
    }
    os << "---------- Grid END -----------------\n";
    return os;
  }
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DGRIDCELL_H
