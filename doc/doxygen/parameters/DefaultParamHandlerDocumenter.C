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

#include <fstream>
#include <OpenMS/ANALYSIS/ID/ConsensusID.h>

using namespace std;
using namespace OpenMS;

//**********************************************************************************
//Helper method - do not change anything here
//**********************************************************************************
void writeParameters(std::ofstream& f, const String& class_name, const Param& param)
{
	f << "/**" << endl << " @page " << class_name << "_Parameters " << class_name << " algorithm parameters" << endl;
	
	String type, description;
	for(map<String,DataValue>::const_iterator it = param.begin(); it != param.end();++it)
	{
		if (it->second.valueType()==DataValue::INTVALUE || it->second.valueType()==DataValue::LONVALUE || it->second.valueType()==DataValue::SHOVALUE  )
		{
			type = "int";
		}
		if (it->second.valueType()==DataValue::FLOVALUE || it->second.valueType()==DataValue::DOUVALUE )
		{
			type = "float";
		}
		if (it->second.valueType()==DataValue::STRVALUE )
		{
			type = "string";
		}
		description = param.getDescription(it->first);
		description.substitute("\n","@n ");
		f <<" - @b "<< it->first << " (" << type << "): " << description << endl;
	}
	f << "*/" << endl;
	f << endl;
}

//**********************************************************************************
//Main method - add your class here
//**********************************************************************************
int main (int, char**)
{
	ofstream f;
	f.open("DefaultParameters.doxygen");
	
	//ConsensusID
	writeParameters(f,"ConsensusID",ConsensusID().getParameters());

  return 0;
}


