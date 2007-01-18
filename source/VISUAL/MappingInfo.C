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
// $Maintainer: Marc Sturm $
// -------------------------------------------------------------------------- 

#include <OpenMS/VISUAL/MappingInfo.h>
#include <OpenMS/FORMAT/Param.h>

namespace OpenMS
{
	using namespace std;
	
	Param MappingInfo::getParam()
	{
		Param	param;
		std::string value;
	
		isMzToXAxis() ? value = "X-Axis" : value = "Y-Axis";
		param.setValue("MappingOfMzTo", value);
	
		isXAxisAsc() ? value = "Ascending" : value = "Descending";
		param.setValue("X-Axis-Orientation", value);
	
		isYAxisAsc() ? value = "Ascending" : value = "Descending";
		param.setValue("Y-Axis-Orientation", value);
	
		return param;
	}
	
	void MappingInfo::setParam(const Param& p) 
	{
		std::string value;
	
		value = p.getValue("MappingOfMzTo").toString();
		if (value == "Y-Axis") setMzToYAxis();
		else if (value == "X-Axis") setMzToXAxis();
		else 
		{
			cerr << "Corrected wrong Axis mapping"<<endl;
			setMzToXAxis();
		}
		
		value = p.getValue("X-Axis-Orientation").toString();
		if (value == "Ascending") setYAxisDesc();
		else if (value == "Descending") setXAxisDesc();
		else 
		{
			setYAxisDesc();
			cerr << "Corrected wrong x axis orientation"<<endl;
		}
		
		value = p.getValue("Y-Axis-Orientation").toString();
		if (value == "Ascending") setYAxisAsc();
		else if (value == "Descending") setYAxisDesc();
		else
		{
			setYAxisAsc();
			cerr << "Corrected wrong y axis orientation"<<endl;
		}
	}

}//Namespace

