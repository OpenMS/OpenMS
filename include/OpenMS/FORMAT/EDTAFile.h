// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_FORMAT_EDTAFILE_H
#define OPENMS_FORMAT_EDTAFILE_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for Enhanced DTA files.
 		
    Input text file containing tab, space or comma separated columns.
    The separator between columns is checked in the first line in this order.

    It supports three variants of this format.

    - Columns are: RT, MZ, Intensity. Header is optional.

    - Columns are: RT, MZ, Intensity, Charge, Meta. Header is mandatory.

      @code
      RT m/z Intensity charge mymeta
      321 405.233 24543534 2 lala
      321 406.207 4343344  2 blubb
      @endcode

    - Columns are: (RT, MZ, Intensity, Charge){1,}, Meta. Header is mandatory.
      First quadruplet is the consensus. All following quadruplets describes the features.

      @code
      RT MZ INT CHARGE RT1 MZ1 INT1 CHARGE1 RT2 MZ2 INT2 CHARGE2
      321 405 100 2 321 405 100 2 321 406 50 2
      323 406 200 2 323 406 200 2 323 407 100 2 323 407 50 2
      @endcode

  	@ingroup FileIO
  */
  class OPENMS_DLLAPI EDTAFile
  {
    public:
      /// Default constructor
      EDTAFile();
			/// Destructor
      virtual ~EDTAFile();
      
    private:
      /**
       * Check if column exists and convert String into DoubleReal.
       */
      DoubleReal checkedToDouble_(const std::vector<String> &parts, Size index, DoubleReal def = -1)
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
      Int checkedToInt_(const std::vector<String> &parts, Size index, Int def = -1)
      {
        if (index < parts.size())
        {
          return parts[index].toInt();
        }
        return def;
      }

    public:
      /**
        @brief Loads a EDTA file into a consensusXML.
 				
 				The content of the file is stored in @p features.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      void load(const String& filename, ConsensusMap& consensus_map)
      {
        // load input
				TextFile input(filename);
		
        // reset map
        ConsensusMap cmap;
        consensus_map = cmap;
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
          if (headers[4].trim() == "RT1") input_type = TYPE_CONSENSUS;
          if (headers[4].trim() != "RT1") input_type = TYPE_OLD_CHARGE;
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
      }

      /**
      	@brief Stores a featureXML as an enhanced DTA file.
      	
        NOT IMPLEMENTED

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum) const
      {
        std::cerr << "Store() for EDTAFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
        throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_EDTAFILE_H

