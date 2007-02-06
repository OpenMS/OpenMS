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
	const std::string FileHandler::NamesOfTypes[] = {"Unknown", "DTA", "DTA2D", "mzData", "mzXML", "FeatureFile", "ANDIMS" };


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
    String input, tmp;
    for (UnsignedInt i=0; i<5; ++i)
    {
			getline(is,tmp,'\n');
			input += tmp;
    }

		//Search for strings
    if (input.find("mzXML")!=string::npos) return MZXML;
    if (input.find("mzData")!=string::npos) return MZDATA;
    if (input.find("featureMap")!=string::npos) return FEATURE;

		return UNKNOWN;
	}

} // namespace OpenMS
