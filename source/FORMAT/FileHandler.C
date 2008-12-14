// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
	const std::string FileHandler::NamesOfTypes[] = {"Unknown", "DTA", "DTA2D", "mzData", "mzXML", "FeatureXML", "cdf", "IdXML", "ConsensusXML", "mgf", "Param", "TrafoXML", "mzML", "ms2"};


	FileHandler::Type FileHandler::getType(const String& filename)
	{
		Type type = getTypeByFileName(filename);
		if (type==UNKNOWN)
		{
			type = getTypeByContent(filename);
		}
		return type;
	}

	FileHandler::Type FileHandler::getTypeByFileName(const String& filename)
	{
		String tmp;
		try
		{
			tmp = filename.suffix('.');
		}
		// no '.' => unknown type
		catch (Exception::ElementNotFound&)
		{
			return UNKNOWN;
		}
		tmp.toUpper();
		if (tmp == "MZDATA")
		{
			return MZDATA;
		}
		else if (tmp == "MZML")
		{
			return MZML;
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
		else if (tmp == "NETCDF")
		{
			return ANDIMS;
		}
		else if (tmp == "FEATUREXML")
		{
			return FEATUREXML;
		}
		else if (tmp == "IDXML")
		{
			return IDXML;
		}
		else if (tmp == "CONSENSUSXML")
		{
			return CONSENSUSXML;
		}
		else if (tmp == "MGF")
		{
			return MGF;
		}
		else if (tmp == "INI")
		{
			return PARAM;
		}
		else if (tmp == "TRAFOXML")
		{
			return TRANSFORMATIONXML;
		}
		else if (tmp == "MS2")
		{
			return MS2;
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
		case MZML:
			return true;
		case MZDATA:
			return true;
		case FEATUREXML:
			return true;
#ifdef USE_ANDIMS
		case ANDIMS:
			return true;
#endif
		case IDXML:
			return true;
		case CONSENSUSXML:
			return true;
		case MGF:
			return true;
		case PARAM:
			return true;
		case TRANSFORMATIONXML:
			return true;
		case MS2:
			return true;
		default:
			return false;
		}
	}

	FileHandler::Type FileHandler::getTypeByContent(const String& filename)
	{
    //load first 5 lines
    TextFile file(filename,true,5);
    file.resize(5); // in case not enough lines are in the file
    String two_five = file[1] + ' ' + file[2] + ' ' + file[3] + ' ' + file[4];
    two_five.substitute('\t',' ');
		String all_simple = file[0] + ' ' + two_five;

		//mzXML (all lines)
    if (all_simple.hasSubstring("<mzXML")) return MZXML;
    
    //mzData (all lines)
    if (all_simple.hasSubstring("<mzData")) return MZDATA;

    //mzML (all lines)
    if (all_simple.hasSubstring("<mzML")) return MZML;
    
    //feature map (all lines)
    if (all_simple.hasSubstring("<featureMap")) return FEATUREXML;

    //ANDIMS (first line)
    if (file[0].hasSubstring("CDF")) return ANDIMS;

    //IdXML (all lines)
    if (all_simple.hasSubstring("<IdXML")) return IDXML;

    //ConsensusXML (all lines)
    if (all_simple.hasSubstring("<consensusXML")) return CONSENSUSXML;

    //mzData (all lines)
    if (all_simple.hasSubstring("<PARAMETERS")) return PARAM;


    //mzData (all lines)
    if (all_simple.hasSubstring("<TrafoXML")) return TRANSFORMATIONXML;

		//tokenize lines two to five
		vector<String> parts;
		two_five.split(' ',parts);
		
		//DTA
		if (parts.size()==8)
		{
			bool conversion_error = false;
			try
			{
				for (UInt i=0; i<8; ++i)
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
				for (UInt i=0; i<12; ++i)
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
		
		// MGF (Mascot Generic Format)
		if (two_five.hasSubstring("BEGIN IONS"))
		{
			return MGF;
		}
		else
		{
			ifstream is(filename.c_str());
			String line;
			while (getline(is, line))
			{
				if (line.hasSubstring("BEGIN IONS"))
				{
					return MGF;
				}
			}
		}

		// MS2 file format
		if (all_simple.hasSubstring("CreationDate"))
		{
			if (all_simple.size() > 0 && all_simple[0] == 'H')
			{
				return MS2;
			}
		}
		return UNKNOWN;
	}

	PeakFileOptions& FileHandler::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& FileHandler::getOptions() const
  {
  	return options_;
  }
} // namespace OpenMS

