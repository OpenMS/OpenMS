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
// $Id: AxisTickCalculator.C,v 1.6 2006/05/31 14:56:42 cfriedle Exp $
// $Author: cfriedle $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------


// STL
#include <iostream>
#include "math.h"
#include <OpenMS/VISUAL/AxisTickCalculator.h>
#include <OpenMS/KERNEL/DSpectrum.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	AxisTickCalculator::GridVector AxisTickCalculator::calcGridLines_(double x1, double x2, int level)
	{
		GridVector gridlines;		
		if (x1 > -0.0001 && x1 < 0.0001) 
			{
				x1 = 0.0001; 
			}
		double dx = x2 - x1;
		if (dx < 0.0000001)
			{
				std::cerr << "Error: grid line intervall too small! Line: " << __LINE__ << " in File " << __FILE__ << std::endl;
				return gridlines;
			}
		double epsilon = dx/200;
		double sDecPow = (int)floor(log10(dx));
		if (sDecPow<0) sDecPow = 0;
		double sDec = pow(10.0,sDecPow);
		std::vector<double> big;
		double currGL = ceil_decimal(x1, (UnsignedInt)sDecPow);
		while (currGL < (x2+epsilon) )
			{
				
				big.push_back(currGL);
				//cout<<"big"<<currGL<<endl;
				currGL += sDec;
			}
		gridlines.push_back(big);
		std::vector<double> small;
		currGL = gridlines[0][0]-sDec/2;
		while(currGL<(x2+epsilon))
			{
				if(currGL>x1)
					small.push_back(currGL);
				//	cout<<"small"<<currGL<<endl;
				currGL +=sDec;
			}
		gridlines.push_back(small);
		if(big.size() <6 && level==3)
			{
					std::vector<double> smaller;
					currGL=gridlines[0][0]-0.75*sDec;
					while(currGL<(x2+epsilon))
						{
							if(currGL>x1)
								smaller.push_back(currGL);
							
							//	cout<<"smaller"<<currGL<<endl;
							currGL +=sDec/2;
						}
					gridlines.push_back(smaller);
			}
		return gridlines;
	}
	
	AxisTickCalculator::GridVector AxisTickCalculator::calcLogGridLines_(double x1, double x2)
	{
		GridVector grid_lines;
		double scalValues[8];
		scalValues[0]=log10(2.0);
		scalValues[1]=log10(3.0);
		scalValues[2]=log10(4.0);
		scalValues[3]=log10(5.0);
		scalValues[4]=log10(6.0);
		scalValues[5]=log10(7.0);
		scalValues[6]=log10(8.0);
		scalValues[7]=log10(9.0);
		double dx = x2-x1;
		if(dx<0.00000001)
			{
				std::cerr<<"Error: grid line intevall too small:"<<__LINE__<<"in File:"<<__FILE__<<std::endl;
				return grid_lines;
			}
		
		int x1ceil = (int)floor(x1);
		int x2floor = (int)ceil(x2);
		std::vector<double> big;
		for(int i = x1ceil;i!=x2floor;++i)
			{
				big.push_back(i);
			}
		grid_lines.push_back(big);
		
		std::vector<double> small;
		for(UnsignedInt i = 0;i!=grid_lines[0].size();++i)
			{
				double currGL =grid_lines[0][i];
				for(int j = 0;j!=8;++j)
					{
						if(currGL + scalValues[j]>x2)break;
						small.push_back(currGL + scalValues[j]);
					}
			}
		grid_lines.push_back(small);
		return grid_lines;
	}

}
