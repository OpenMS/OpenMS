// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


// STL
#include <iostream>
#include <cmath>

// OpenMS
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

#include <boost/math/special_functions/fpclassify.hpp>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	void AxisTickCalculator::calcGridLines(DoubleReal x1, DoubleReal x2, Int levels, GridVector& grid, UInt max_num_big, UInt max_num_small)
	{		
		grid.clear();
		
		if (boost::math::isnan(x1) || boost::math::isnan(x2)) return;
		
		if (x1 > -0.0001 && x1 < 0.0001) 
		{
			x1 = 0.0001; 
		}

		DoubleReal dx = x2 - x1;

		if (dx < 0.0000001)
		{
			return ;
		}
		DoubleReal epsilon = dx/200;

		DoubleReal sDecPow = floor(log10(dx));
		DoubleReal sDec = pow(10.0,sDecPow);
		std::vector<DoubleReal> big;
		DoubleReal currGL = ceilDecimal(x1, (UInt)sDecPow);
		while (currGL < (x2+epsilon) )
		{			
			big.push_back(currGL);
			currGL += sDec;
		}
		grid.push_back(big);
		if (big.size() < max_num_big && levels>=2) 
		{
			std::vector<DoubleReal> small;
			currGL = grid[0][0]-sDec/2;
			while(currGL<(x2+epsilon))
			{
				if(currGL>x1)
				{
					small.push_back(currGL);
				}
				currGL +=sDec;
			}
			
			grid.push_back(small);

			if(big.size() < max_num_small && levels==3)
			{	
				std::vector<DoubleReal> smaller;
				currGL=grid[0][0]-0.75*sDec;
				while(currGL<(x2+epsilon))
				{
					if(currGL>x1)
					{
						smaller.push_back(currGL);
					}
					currGL +=sDec/2;
				}
				grid.push_back(smaller);
			}
		}
	}
	
	void AxisTickCalculator::calcLogGridLines(DoubleReal x1, DoubleReal x2, GridVector& grid)
	{
		grid.clear();		
		DoubleReal scalValues[8];
		scalValues[0]=log10(2.0);
		scalValues[1]=log10(3.0);
		scalValues[2]=log10(4.0);
		scalValues[3]=log10(5.0);
		scalValues[4]=log10(6.0);
		scalValues[5]=log10(7.0);
		scalValues[6]=log10(8.0);
		scalValues[7]=log10(9.0);
		DoubleReal dx = x2-x1;

		if(dx<0.00000001)
		{
			return;
		}
		
		Int x1ceil = (Int)floor(x1);
		Int x2floor = (Int)ceil(x2);
		std::vector<DoubleReal> big;
		for(Int i = x1ceil;i!=x2floor;++i)
		{
			big.push_back(i);
		}
		grid.push_back(big);
		std::vector<DoubleReal> small;
		for (Size i = 0;i!=grid[0].size();++i)
		{
			DoubleReal currGL =grid[0][i];
			for(Int j = 0;j!=8;++j)
			{
				if(currGL + scalValues[j]>x2)
				{
					break;
				}
				small.push_back(currGL + scalValues[j]);
			}
		}
		grid.push_back(small);
	}

} //namespace


