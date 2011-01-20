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

#ifndef OPENMS_FORMAT_MSINSPECTFILE_H
#define OPENMS_FORMAT_MSINSPECTFILE_H

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
 		@brief File adapter for MsInspect files.
 		
    Lines with "#" are comments and are ignored.

  	The first non-comment line is the header and contains the column names:<br>
    scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description

    Every subsequent line is a feature.
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI MsInspectFile
  {
    public:
      /// Default constructor
      MsInspectFile();
			/// Destructor
      virtual ~MsInspectFile();
      
      /**
 				@brief Loads a MsInspect file into a featureXML.
 				
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

				bool first_line = true;
				for (Size i=1; i<input.size(); ++i)
				{
					String line = input[i];
					
					//ignore comment lines
					if (line.empty() || line[0]=='#') continue;
			
					//skip leader line
					if (first_line)
					{
						first_line = false;
						continue;
					}
					
					//split lines: scan	time	mz	accurateMZ	mass	intensity	charge	chargeStates	kl	background	median	peaks	scanFirst	scanLast	scanCount	totalIntensity	sumSquaresDist	description
					std::vector< String > parts;
					line.split('\t', parts);
					
          if (parts.size()<18) 
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"", String("Failed to convert line ")  + String(i+1) + "not enough columns (expected 18 or more, got " + String(parts.size()) + ")");
          }

					//create feature
					Feature f;
					Size column_to_convert=0;
					try
					{
						column_to_convert = 1;
						f.setRT(parts[1].toDouble());
						column_to_convert = 2;
						f.setMZ(parts[2].toDouble());
						column_to_convert = 5;
						f.setIntensity(parts[5].toDouble());
						column_to_convert = 6;
						f.setCharge(parts[6].toInt());
						column_to_convert = 8;
						f.setOverallQuality(parts[8].toDouble());

						column_to_convert = 3;
						f.setMetaValue("accurateMZ",parts[3]);
						column_to_convert = 4;
						f.setMetaValue("mass",parts[4].toDouble());
						column_to_convert = 7;
						f.setMetaValue("chargeStates",parts[7].toInt());
						column_to_convert = 9;
						f.setMetaValue("background",parts[9].toDouble());
						column_to_convert = 10;
						f.setMetaValue("median",parts[10].toDouble());
						column_to_convert = 11;
						f.setMetaValue("peaks",parts[11].toInt());
						column_to_convert = 12;
						f.setMetaValue("scanFirst",parts[12].toInt());
						column_to_convert = 13;
						f.setMetaValue("scanLast",parts[13].toInt());
						column_to_convert = 14;
						f.setMetaValue("scanCount",parts[14].toInt());
						column_to_convert = 15;
						f.setMetaValue("totalIntensity",parts[15].toDouble());
						column_to_convert = 16;
						f.setMetaValue("sumSquaresDist",parts[16].toDouble());
					}
					catch (Exception::BaseException /*&e*/)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"", String("Failed to convert value in column ") + String(column_to_convert+1) + "into a number (line '" + (i+1) + ")");
					}
					f.setMetaValue("description",parts[17]);
					feature_map.push_back(f);
				}
			
      }

      /**
      	@brief Stores a featureXML as a MsInspect file.
      	
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

#endif // OPENMS_FORMAT_MSINSPECTFILE_H

