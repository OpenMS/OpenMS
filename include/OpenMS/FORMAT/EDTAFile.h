// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <fstream>
#include <vector>

namespace OpenMS
{
 	/**
 		@brief File adapter for Enhanced DTA files.
 		
  	Input text file containing the following columns: RT m/z intensity.");
		Additionally meta data columns may follow.
		If meta data is used, meta data column names have to be specified in a header line, e.g.
@code
    RT m/z Intensity charge mymeta
    321 405.233 24543534 2 lala
    321 406.207 4343344  2 blubb
@endcode

    The separator between columns is checked in the first line in this order:
    Tab, Space, Comma

		If a meta column named 'charge' with numeric data exists, the charge of the features will be set accordingly.
    Every subsequent line is a feature.
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI EDTAFile
  {
    public:
      /// Default constructor
      EDTAFile();
			/// Destructor
      virtual ~EDTAFile();
      
      /**
 				@brief Loads a EDTA file into a featureXML.
 				
 				The content of the file is stored in @p features.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      template <typename FeatureMapType>
      void load(const String& filename, FeatureMapType& feature_map)
      {
        // load input
				TextFile input(filename);
		
				// reset map
        FeatureMapType fmap;
				feature_map = fmap;
			
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
				DoubleReal rt = 0.0;
				DoubleReal mz = 0.0;
				DoubleReal it = 0.0;
				// see if we have a header
				try
				{
          if (headers.size()>3) throw Exception::BaseException(); // there is meta-data, so these must be their names
          else if (headers.size()<3) throw Exception::BaseException(); // not enough data columns in first line...
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

				// parsing features
				feature_map.reserve(input.size());
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
					
					//convert coordinate columns to doubles
					try
					{
						rt = parts[0].toDouble();
						mz = parts[1].toDouble();
						it = parts[2].toDouble();
					}
					catch (Exception::BaseException&)
					{
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "",
              String("Failed parsing in line") + String(i+1) + ": Could not convert the first three columns to float! Is the correct separator specified?\nOffending line: '" + line_trimmed + "'  (line " + (i+1) + ")\n");
					}
					Feature f;
					f.setMZ(mz);
					f.setRT(rt);
					f.setIntensity(it);

					//parse meta data
					for (Size j=3; j<parts.size(); ++j)
					{
						String part_trimmed = parts[j];
						part_trimmed.trim();
						if (part_trimmed!="")
						{
							//check if column name is ok
							if (headers.size()<=j || headers[j]=="")
							{
                throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "",
								  String("Error: Missing meta data header for column ") + (j+i) + "!"
                  + String("Offending header line: '") + header_trimmed + "'  (line 1)");
							}
							//add meta value
							f.setMetaValue(headers[j],part_trimmed);
							if (headers[j] == "charge")
							{
								try
								{
									f.setCharge(part_trimmed.toInt());
								}
								catch (...)
								{
									LOG_WARN << "Failed to convert metavalue 'charge' into integer (line '" << (i+1) << ")";
								}
							}
						}

					}
					
					//insert feature to map
					feature_map.push_back(f);
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
        throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_EDTAFILE_H

