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

#include <OpenMS/FORMAT/FileHandler.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
	const std::string FileHandler::NamesOfTypes[] = {"Unknown", "DTA", "DTA2D", "mzData", "mzXML", "FeatureFile", "FeaturePairs", "ANDIMS" };


	FileHandler::Type FileHandler::getTypeByFileName(const String& filename)
	{
		String tmp;
		try
		{
			tmp = filename.suffix('.');
		}
		// no '.' => unknown type
		catch (Exception::ElementNotFound<char>)
		{
			return UNKNOWN;
		}
		tmp.toUpper();
		if (tmp == "MZDATA")
		{
			return MZDATA;
		}
		else if (tmp == "DTA")
		{
			return DTA;
		}
		else if (tmp == "DTA2D")
		{
			return DTA2D;
		}
		else if (tmp == "MZXML")
		{
			return MZXML;
		}
		else if (tmp == "CDF")
		{
			return ANDIMS;
		}
		else if (tmp == "FEAT")
		{
			return FEATURE;
		}
		else if (tmp == "PAIRS")
		{
			return FEATURE_PAIRS;
		}

		return UNKNOWN;

	}

	FileHandler::Type FileHandler::nameToType(const String& name)
	{
		String tmp = name;
		tmp.toUpper();
		String tmp2;

		for (int i=0; i < SIZE_OF_TYPE; ++i)
		{
			tmp2 = NamesOfTypes[i];
			tmp2.toUpper();
			if (tmp == tmp2)
			{
				return (Type)i;
			}
		}

		return UNKNOWN;
	}

	String FileHandler::typeToName(Type type)
	{
		return NamesOfTypes[type];
	}

	bool FileHandler::isSupported(Type type)
	{
		switch (type)
		{
		case DTA:
			return true;
		case DTA2D:
			return true;
		case MZXML:
			return true;
		case MZDATA:
			return true;
		case FEATURE:
			return true;
		case FEATURE_PAIRS:
			return true;
#ifdef ANDIMS_DEF
		case ANDIMS:
			return true;
#endif
		default:
			return false;
		}
	}

	FileHandler::Type FileHandler::getTypeByContent(const String& filename) throw (Exception::FileNotFound)
	{
		ifstream is(filename.c_str());
    if (!is)
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    //load first 5 lines
    String one, two, three, four, five;
    getline(is,one,'\n');
    one.trim();
    getline(is,two,'\n');
    two.trim();
    getline(is,three,'\n');
    three.trim();
    getline(is,four,'\n');
   	four.trim();
    getline(is,five,'\n');
    five.trim();
    // concatenate trimmed lines
    String two_five = two + ' ' + three + ' ' + four + ' ' + five;
    // replace tabs by spaces
    two_five.substitute('\t',' ');
    
		//mzXML (all lines)
    if ((one + ' ' + two_five).find("mzXML")!=string::npos) return MZXML;
    
    //mzData (all lines)
    if ((one + ' ' + two_five).find("mzData")!=string::npos) return MZDATA;
    
    //feature map (all lines)
    if ((one + ' ' + two_five).find("featureMap")!=string::npos) return FEATURE;

    //feature pairs (all lines)
    if ((one + ' ' + two_five).find("featurePairs")!=string::npos) return FEATURE_PAIRS;
        
    //ANDIMS (first line)
    if (one.find("CDF")!=string::npos) return ANDIMS;

		//tokenize lines two to five
		vector<String> parts;
		two_five.split(' ',parts);
		
		//DTA
		if (parts.size()==8)
		{
			bool conversion_error = false;
			try
			{
				for (UnsignedInt i=0; i<8; ++i)
				{
					parts[i].toFloat();
				}
			}
			catch (Exception::ConversionError)
			{
				conversion_error = true;
			}
			if (!conversion_error) return DTA;
		}
		
		//DTA2D
		if (parts.size()==12)
		{
			bool conversion_error = false;
			try
			{
				for (UnsignedInt i=0; i<12; ++i)
				{
					parts[i].toFloat();
				}
			}
			catch (Exception::ConversionError)
			{
				conversion_error = true;
			}
			if (!conversion_error) return DTA2D;
		}
		
		return UNKNOWN;
	}

} // namespace OpenMS
