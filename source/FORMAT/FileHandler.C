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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/FORMAT/GzipIfstream.h>
#include <OpenMS/FORMAT/Bzip2Ifstream.h>

#include <fstream>

using namespace std;

namespace OpenMS
{

	const std::string FileHandler::NamesOfTypes[] = {"Unknown", "DTA", "DTA2D", "mzData", "mzXML", "FeatureXML", "cdf", "IdXML", "ConsensusXML", "mgf", "ini", "TrafoXML", "mzML", "ms2", "pepXML", "mzIdentML", "GelML", "TraML", "MSP", "OMSSAXML", "PNG"};

	FileTypes::Type FileHandler::getType(const String& filename)
	{
		FileTypes::Type type = getTypeByFileName(filename);
		if (type==FileTypes::UNKNOWN)
		{
			type = getTypeByContent(filename);
		}
		return type;
	}

	FileTypes::Type FileHandler::getTypeByFileName(const String& filename)
	{
		String tmp;
		try
		{
			tmp = filename.suffix('.');
		}
		// no '.' => unknown type
		catch (Exception::ElementNotFound&)
		{
			return FileTypes::UNKNOWN;
		}
		tmp.toUpper();
		if (tmp == "MZDATA")
		{
			return FileTypes::MZDATA;
		}
		else if (tmp == "MZML")
		{
			return FileTypes::MZML;
		}
		else if (tmp == "DTA")
		{
			return FileTypes::DTA;
		}
		else if (tmp == "DTA2D")
		{
			return FileTypes::DTA2D;
		}
		else if (tmp == "MZXML")
		{
			return FileTypes::MZXML;
		}
		else if (tmp == "CDF")
		{
			return FileTypes::ANDIMS;
		}
		else if (tmp == "NETCDF")
		{
			return FileTypes::ANDIMS;
		}
		else if (tmp == "FEATUREXML")
		{
			return FileTypes::FEATUREXML;
		}
		else if (tmp == "IDXML")
		{
			return FileTypes::IDXML;
		}
		else if (tmp == "CONSENSUSXML")
		{
			return FileTypes::CONSENSUSXML;
		}
		else if (tmp == "MGF")
		{
			return FileTypes::MGF;
		}
		else if (tmp == "INI")
		{
			return FileTypes::INI;
		}
		else if (tmp == "TRAFOXML")
		{
			return FileTypes::TRANSFORMATIONXML;
		}
		else if (tmp == "MS2")
		{
			return FileTypes::MS2;
		}
		else if (tmp == "PEPXML") 
		{
			return FileTypes::PEPXML;
		}
		else if (tmp == "MZIDENTML")
		{
			return FileTypes::MZIDENTML;
		}
		else if (tmp == "GELML")
		{
			return FileTypes::GELML;
		}
		else if (tmp == "TRAML")
		{
			return FileTypes::TRAML;
		}
		else if (tmp == "MSP")
		{
			return FileTypes::MSP;
		}
		else if (tmp == "PNG")
		{
			return FileTypes::PNG;
		}
		else if (tmp == "BZ2" || tmp == "ZIP" || tmp == "GZ")
		{
			return getTypeByContent(filename);
		}

		return FileTypes::UNKNOWN;

	}

	FileTypes::Type FileHandler::nameToType(const String& name)
	{
		String tmp = name;
		tmp.toUpper();
		String tmp2;

		for (int i=0; i < FileTypes::SIZE_OF_TYPE; ++i)
		{
			tmp2 = NamesOfTypes[i];
			tmp2.toUpper();
			if (tmp == tmp2)
			{
				return (FileTypes::Type)i;
			}
		}

		return FileTypes::UNKNOWN;
	}

	String FileHandler::typeToName(FileTypes::Type type)
	{
		return NamesOfTypes[type];
	}

	bool FileHandler::isSupported(FileTypes::Type type)
	{
		switch (type)
		{
		case FileTypes::DTA:
			return true;
		case FileTypes::DTA2D:
			return true;
		case FileTypes::MZXML:
			return true;
		case FileTypes::MZML:
			return true;
		case FileTypes::MZDATA:
			return true;
		case FileTypes::FEATUREXML:
			return true;
#ifdef USE_ANDIMS
		case FileTypes::ANDIMS:
			return true;
#endif
		case FileTypes::IDXML:
			return true;
		case FileTypes::CONSENSUSXML:
			return true;
		case FileTypes::MGF:
			return true;
		case FileTypes::INI:
			return true;
		case FileTypes::TRANSFORMATIONXML:
			return true;
		case FileTypes::MS2:
			return true;
		case FileTypes::MZIDENTML:
			return true;
		case FileTypes::PEPXML:
			return true;
		case FileTypes::GELML:
			return true;
		case FileTypes::OMSSAXML:
			return true;
		case FileTypes::PNG:
			return true;
		case FileTypes::TRAML:
		default:
			return false;
		}
	}

	FileTypes::Type FileHandler::getTypeByContent(const String& filename)
	{
		String first_line;
		String two_five;
		String all_simple;

		// only the first five lines will be set for compressed files
		// so far, compression is only supported for XML files
		vector<String> complete_file;


		// test whether the file is compressed (bzip2 or gzip)
    ifstream compressed_file(filename.c_str());
    char bz[2];
    compressed_file.read(bz,2);
    char g1 = 0x1f;
    char g2 = 0;
    g2 |= 1 << 7;
    g2 |= 1 <<3;
    g2  |=1  <<1;
    g2 |=1 <<0;
		compressed_file.close();
    if(bz[0] == 'B' && bz[1] =='Z' ) // bzip2
		{
			Bzip2Ifstream bzip2_file(filename.c_str());
			char buffer[1024];
			bzip2_file.read(buffer, 1024);
			String buffer_str(buffer);
      vector<String> split;
      buffer_str.split('\n', split);
      split.resize(5);
      first_line = split[0];
      two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
			all_simple = first_line + ' ' + two_five;
      complete_file = split;
		}
		else if (bz[0] == g1 && bz[1] == g2) // gzip
    {
			GzipIfstream gzip_file(filename.c_str());
			char buffer[1024];
			gzip_file.read(buffer, 1024);
			String buffer_str(buffer);
			vector<String> split;
			buffer_str.split('\n', split);
			split.resize(5);
			first_line = split[0];
			two_five = split[1] + ' ' + split[2] + ' ' + split[3] + ' ' + split[4];
			all_simple = first_line + ' ' + two_five;
			complete_file = split;
    }		
		else // uncompressed
		{
    	//load first 5 lines
    	TextFile file(filename, true, 5);
    	file.resize(5); // in case not enough lines are in the file
    	two_five = file[1] + ' ' + file[2] + ' ' + file[3] + ' ' + file[4];
    	two_five.substitute('\t',' ');
			all_simple = file[0] + ' ' + two_five;
			first_line = file[0];
			complete_file = file;
		}

		//mzXML (all lines)
    if (all_simple.hasSubstring("<mzXML")) return FileTypes::MZXML;

    //mzData (all lines)
    if (all_simple.hasSubstring("<mzData")) return FileTypes::MZDATA;

    //mzML (all lines)
    if (all_simple.hasSubstring("<mzML")) return FileTypes::MZML;

		//analysisXML (all lines)
		if (all_simple.hasSubstring("<mzIdentML")) return FileTypes::MZIDENTML;

		//pepXML (all lines)
		if (all_simple.hasSubstring("xmlns=\"http://regis-web.systemsbiology.net/pepXML\"")) return FileTypes::PEPXML;

    //feature map (all lines)
    if (all_simple.hasSubstring("<featureMap")) return FileTypes::FEATUREXML;

    //ANDIMS (first line)
    if (first_line.hasSubstring("CDF")) return FileTypes::ANDIMS;

    //IdXML (all lines)
    if (all_simple.hasSubstring("<IdXML")) return FileTypes::IDXML;

    //ConsensusXML (all lines)
    if (all_simple.hasSubstring("<consensusXML")) return FileTypes::CONSENSUSXML;

    //mzData (all lines)
    if (all_simple.hasSubstring("<PARAMETERS")) return FileTypes::INI;

    //mzData (all lines)
    if (all_simple.hasSubstring("<TrafoXML")) return FileTypes::TRANSFORMATIONXML;

		//GelML (all lines)
		if (all_simple.hasSubstring("<GelML")) return FileTypes::GELML;

		//traML (all lines)
		if (all_simple.hasSubstring("<TraML")) return FileTypes::TRAML;
	
		//OMSSAXML file
		if (all_simple.hasSubstring("<MSResponse")) return FileTypes::OMSSAXML;

		// PNG file (to be really correct, the first eight bytes of the file would
		// have to be checked; see e.g. the wikipedia article)
		if (first_line.substr(1, 3) == "PNG") return FileTypes::PNG;

		//MSP (all lines)
		for (Size i = 0; i != complete_file.size(); ++i)
		{
			if (complete_file[i].hasPrefix("Name: ") && complete_file[i].hasSubstring("/"))
			{
				return FileTypes::MSP;
			}
			if (complete_file[i].hasPrefix("Num peaks: "))
			{
				return FileTypes::MSP;
			}
		}

		//tokenize lines two to five
		vector<String> parts;
		two_five.split(' ',parts);

		//DTA
		if (parts.size()==8)
		{
			bool conversion_error = false;
			try
			{
				for (Size i=0; i<8; ++i)
				{
					parts[i].toFloat();
				}
			}
			catch (Exception::ConversionError)
			{
				conversion_error = true;
			}
			if (!conversion_error) return FileTypes::DTA;
		}

		//DTA2D
		if (parts.size()==12)
		{
			bool conversion_error = false;
			try
			{
				for (Size i=0; i<12; ++i)
				{
					parts[i].toFloat();
				}
			}
			catch (Exception::ConversionError)
			{
				conversion_error = true;
			}
			if (!conversion_error) return FileTypes::DTA2D;
		}

		// MGF (Mascot Generic Format)
		if (two_five.hasSubstring("BEGIN IONS"))
		{
			return FileTypes::MGF;
		}
		else
		{
			for (Size i = 0; i != complete_file.size(); ++i)
			{
				if (complete_file[i].hasSubstring("BEGIN IONS"))
				{
					return FileTypes::MGF;
				}
			}
		}

		// MS2 file format
		if (all_simple.hasSubstring("CreationDate"))
		{
			if (all_simple.size() > 0 && all_simple[0] == 'H')
			{
				return FileTypes::MS2;
			}
		}

		return FileTypes::UNKNOWN;
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

