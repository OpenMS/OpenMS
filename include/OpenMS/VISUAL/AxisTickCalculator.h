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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_AXISTICKCALCULATOR_H
#define OPENMS_VISUAL_AXISTICKCALCULATOR_H

// STL
#include <vector>
#include <OpenMS/CONCEPT/Types.h>

namespace OpenMS
{
	/**
		@brief Calculates ticks for a given value range.
		
		It has only static methods, that's by the constructor is private.
	
		@ingroup Visual
	*/
  class AxisTickCalculator 
	{
		public:
	 
	  /// Typedef for the grid vector
		typedef std::vector<std::vector<double> > GridVector;

		/**
			 @brief Returns a GridVector with ticks for linear scales.
			 
			 @param x1 minimum value
			 @param x2 maximum value
			 @param levels numbers of different tick levels (maximum is 3)
			 @param grid the grid_vector to fill
			 @param max_num_big
			 @param max_num_small 
			 @param grid_line_dist the distance of the gridlines
		*/
		static void calcGridLines(double x1, double x2, int levels, GridVector& grid, UnsignedInt max_num_big, UnsignedInt max_num_small, double& grid_line_dist);
		
		/**
			 @brief Returns a GridVector with ticks for logarithmic scales.
			 
			 @param x1 minimum value
			 @param x2 maximum value
			 @param grid the grid_vector to fill
		*/
		static void calcLogGridLines(double x1, double x2, GridVector& grid);
			
		private: 
		
		///Constructor: only static methods
		AxisTickCalculator();
	};
}
#endif
