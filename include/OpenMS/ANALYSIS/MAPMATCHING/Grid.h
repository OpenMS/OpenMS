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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#ifndef OPENMS_ANALYSIS_MAPMATCHING_GRID_H
#define OPENMS_ANALYSIS_MAPMATCHING_GRID_H

#include <OpenMS/ANALYSIS/MAPMATCHING/GridCell.h>

#include <vector>

namespace OpenMS
{
	/**
		@brief This class represents a D-dimensional grid over a LC/MS map.	
	
	*/
	class Grid : public std::vector<GridCell>
	{
	 public:
			
		/**	
			 @name Type definitions
		*/
		//@{
		typedef std::vector<GridCell> Base;
		typedef Base::iterator Iterator;
		typedef Base::const_iterator ConstIterator;
		typedef Base::reverse_iterator ReverseIterator;
		typedef Base::const_reverse_iterator ConstReverseIterator;
		typedef GridCell& Reference;
		typedef const GridCell& ConstReference;
			
		// STL compatibility
		typedef GridCell value_type;
		typedef GridCell* pointer;
		typedef const GridCell* const_pointer;
		typedef Reference reference;
		typedef ConstReference const_reference;
		typedef Base::size_type size_type;
		typedef Base::difference_type difference_type;
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
		Grid() {}
		/// Copy constructor
		Grid(const Grid& grid)
			: Base(grid)
		{}
		/// Destructor
		virtual ~Grid() {}
		//@}
			
		/// Assignment operator
   Grid& operator = (const Grid& rhs);
	}; // end of class Grid
	
} // end of namespace OpenMS

#endif  // OPENMS_ANALYSIS_MAPMATCHER_DGRIDCELL_H
