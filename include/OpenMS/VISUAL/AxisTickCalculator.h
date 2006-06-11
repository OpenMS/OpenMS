// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: AxisTickCalculator.h,v 1.4 2006/03/03 18:40:07 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_AXISTICKCALCULATOR_H
#define OPENMS_VISUAL_AXISTICKCALCULATOR_H

// STL
#include <vector>
#include <math.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>
namespace OpenMS
{
	/**
		@brief Calculates ticks for a given value range.
		
		It has only static methods, that's by the constructor is private.
	
		@ingroup Visual
		
		@todo write tests (Cornelia)
	*/
  class AxisTickCalculator 
	{
		
		public:
	    ///
	    enum {TOP, BOTTOM, LEFT, RIGHT} ALIGNMENT_ENUM;
	    ///
			typedef std::vector<std::vector<double> > GridVector;
			
			/**
				@brief Returns a GridVector with ticks for linear scales.
				
				@param x1 minimum value
				@param x2 maximum value
				@param level numbers of different tick levels (maximum is 3)
			*/
			static GridVector calcGridLines_(double x1, double x2, int level);
			/**
				@brief Returns a GridVector with ticks for logarithmic scales.
				
				@param x1 minimum value
				@param x2 maximum value
			*/
			static GridVector calcLogGridLines_(double x1, double x2);
			
		private: 
			///Constructor: only static methods
			AxisTickCalculator();

			///Destructor
			~AxisTickCalculator();
	};
}
#endif
