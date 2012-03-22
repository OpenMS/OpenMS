// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/EDTAFile.h>
#include <cmath>

using namespace std;

namespace OpenMS
{

	EDTAFile::EDTAFile()
	{
	}

	EDTAFile::~EDTAFile()
	{
	}

  DoubleReal EDTAFile::checkedToDouble_(const std::vector<String> &parts, Size index, DoubleReal def)
  {
    if (index < parts.size())
    {
      return parts[index].toDouble();
    }
    return def;
  }

  /**
    * Check if column exists and convert String into Int.
    */
  Int EDTAFile::checkedToInt_(const std::vector<String> &parts, Size index, Int def)
  {
    if (index < parts.size())
    {
      return parts[index].toInt();
    }
    return def;
  }

  /**
    @brief Loads a EDTA file into a consensusXML.
 				
 		The content of the file is stored in @p features.

		@exception Exception::FileNotFound is thrown if the file could not be opened
		@exception Exception::ParseError is thrown if an error occurs during parsing
  */
  void EDTAFile::load(const String& filename, ConsensusMap& consensus_map)
  {
    // load input
		TextFile input(filename);
		
    // reset map
    consensus_map = ConsensusMap();
    consensus_map.setUniqueId();

    char separator = ' ';
    if (input[0].hasSubstring("\t")) separator = '\t';
    else if (input[0].hasSubstring(" ")) separator = ' ';
    else if (input[0].hasSubstring(",")) separator = ',';

		// parsing header line
		std::vector<String> headers;
		input[0].split(separator, headers);
		int offset = 0;
		for (Size i=0; i<headers.size(); ++i)
		{
			headers[i].trim();
		}
		String header_trimmed = input[0];
		header_trimmed.trim();

    enum
    {
      TYPE_UNDEFINED,
      TYPE_OLD_NOCHARGE,
      TYPE_OLD_CHARGE,
      TYPE_CONSENSUS
    }
    input_type = TYPE_UNDEFINED;
    Size input_features = 1;

    DoubleReal rt = 0.0;
    DoubleReal mz = 0.0;
    DoubleReal it = 0.0;
    Int ch = 0;

    if (headers.size() <= 2)
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Failed parsing in line 1: not enough columns! Expected at least 3 columns!\nOffending line: '") + header_trimmed + "'  (line 1)\n");
    }
    else if (headers.size() == 3) input_type = TYPE_OLD_NOCHARGE;
    else if (headers.size() == 4) input_type = TYPE_OLD_CHARGE;

		// see if we have a header
		try
		{
      // try to convert... if not: thats a header
      rt = headers[0].toDouble();
			mz = headers[1].toDouble();
			it = headers[2].toDouble();
		}
		catch (Exception::BaseException&)
		{
			offset=1;
			LOG_INFO << "Detected a header line.\n";
		}
        
    if (headers.size() >= 5)
    {
      if (String(headers[4].trim()).toUpper() == "RT1") input_type = TYPE_CONSENSUS;
      else input_type = TYPE_OLD_CHARGE;
    }
    if (input_type == TYPE_CONSENSUS)
    {
      // Every consensus style line includes features with four columns.
      // The remainder is meta data
      input_features = headers.size() / 4;
    }

    if (offset==0 && (input_type==TYPE_OLD_CHARGE || input_type==TYPE_CONSENSUS))
    {
      throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Failed parsing in line 1: No HEADER provided. This is only allowed for three columns. You have more!\nOffending line: '") + header_trimmed + "'  (line 1)\n");
    }

    ConsensusMap::FileDescription desc;
    desc.filename = filename;
    desc.size = input.size() - offset;
    consensus_map.getFileDescriptions()[0] = desc;

    // parsing features
    consensus_map.reserve(input.size());

		for (Size i=offset; i<input.size(); ++i)
		{
			//do nothing for empty lines
			String line_trimmed = input[i];
			line_trimmed.trim();
			if (line_trimmed=="")
			{
				if (i<input.size()-1) LOG_WARN << "Notice: Empty line ignored (line " << (i+1) << ").";
				continue;
			}
					
			//split line to tokens
			std::vector<String> parts;
			input[i].split(separator, parts);

			//abort if line does not contain enough fields
			if (parts.size()<3)
			{
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Failed parsing in line ") + String(i+1) + ": At least three columns are needed! (got  " + String(parts.size()) + ")\nOffending line: '" + line_trimmed + "'  (line " + (i+1) + ")\n");
			}
					
      ConsensusFeature cf;
      cf.setUniqueId();

      try
      {
        // Convert values. Will return -1 if not available.
        rt = checkedToDouble_(parts, 0);
        mz = checkedToDouble_(parts, 1);
        it = checkedToDouble_(parts, 2);
        ch = checkedToInt_(parts, 3);

        cf.setRT(rt);
        cf.setMZ(mz);
        cf.setIntensity(it);
        if (input_type != TYPE_OLD_NOCHARGE) cf.setCharge(ch);

        // Check all features in one line
        for (Size j = 1; j < input_features; ++j)
        {
          Feature f;
          f.setUniqueId();

          // Convert values. Will return -1 if not available.
          rt = checkedToDouble_(parts, j * 4 + 0);
          mz = checkedToDouble_(parts, j * 4 + 1);
          it = checkedToDouble_(parts, j * 4 + 2);
          ch = checkedToInt_(parts, j * 4 + 3);

          // Only accept features with at least RT and MZ set
          if (rt != -1 && mz != -1)
          {
            f.setRT(rt);
            f.setMZ(mz);
            f.setIntensity(it);
            f.setCharge(ch);

            cf.insert(j-1, f);
          }
        }
      }
      catch (Exception::BaseException&)
      {
        throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Failed parsing in line") + String(i + 1) + ": Could not convert the first three columns to float! Is the correct separator specified?\nOffending line: '" + line_trimmed + "'  (line " + (i + 1) + ")\n");
      }

 			//parse meta data
      for (Size j = input_features * 4; j < parts.size(); ++j)
			{
				String part_trimmed = parts[j];
				part_trimmed.trim();
				if (part_trimmed!="")
				{
					//check if column name is ok
					if (headers.size()<=j || headers[j]=="")
					{
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "",
							String("Error: Missing meta data header for column ") + (j+1) + "!"
              + String("Offending header line: '") + header_trimmed + "'  (line 1)");
					}
					//add meta value
					cf.setMetaValue(headers[j],part_trimmed);
				}
      }

      //insert feature to map
      consensus_map.push_back(cf);
		}
    
    // register FileDescriptions
    ConsensusMap::FileDescription fd;
    fd.filename = filename;
    fd.size = consensus_map.size();
    Size maps = std::max(input_features-1, Size(1)); // its either a simple feature or a consensus map
                                               // (in this case the 'input_features' includes the centroid, which we do not count)
    for (Size i = 0; i < maps; ++i)
    {
      fd.label = String("EDTA_Map ") + String(i);
      consensus_map.getFileDescriptions()[i] = fd;
    }

  }

  /**
    @brief Stores a ConsensusMap as an enhanced DTA file.
      	
    NOT IMPLEMENTED

		@exception Exception::UnableToCreateFile is thrown if the file could not be created
  */
  void EDTAFile::store(const String& filename, const ConsensusMap& map) const
  {
    std::cerr << "Store() for EDTAFile not implemented. Filename was: " << filename << ", CM of size " << map.size() << "\n";
    throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
  }
} // namespace OpenMS

