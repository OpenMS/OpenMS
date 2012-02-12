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

#ifndef OPENMS_FORMAT_SPECARRAYFILE_H
#define OPENMS_FORMAT_SPECARRAYFILE_H

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
 		@brief File adapter for SpecArray (.pepList) files.
 		
  	The first line is the header and contains the column names:<br>
           m/z	     rt(min)	       snr	      charge	   intensity

    Every subsequent line is a feature.
    Entries are separated by Tab (\\t).
    
  	
  	@ingroup FileIO
  */
  class OPENMS_DLLAPI SpecArrayFile
  {
    public:
      /// Default constructor
      SpecArrayFile();
			/// Destructor
      virtual ~SpecArrayFile();
      
      /**
 				@brief Loads a SpecArray file into a featureXML.
 				
 				The content of the file is stored in @p features.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
      */
      template <typename FeatureMapType>
      void load(const String& filename, FeatureMapType& feature_map)
      {
        // load input
				TextFile input(filename,false);
		
				// reset map
        FeatureMapType fmap;
				feature_map = fmap;
			
				for (Size i=1; i<input.size(); ++i)
				{
					String line = input[i];

					std::vector< String > parts;
					line.split('\t', parts);
					
          if (parts.size()<5) 
          {
            throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__,"", String("Failed to convert line")  + String(i+1) + "not enough columns (expected 5 or more, got " + String(parts.size()) + ")");
          }

					Feature f;
					try
					{
						f.setMZ(parts[0].toDouble());
						f.setRT(parts[1].toDouble() *60.0);
						f.setMetaValue("s/n",parts[2].toDouble());
						f.setCharge(parts[3].toInt());
						f.setIntensity(parts[4].toDouble());
					}
					catch (Exception::BaseException /*&e*/)
					{
						throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, "", String("Failed to convert value into a number (line '") + (i+1) + ")");
					}
					feature_map.push_back(f);
				}
      }

      /**
      	@brief Stores a featureXML as a SpecArray file.
      	
        NOT IMPLEMENTED

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
      */
      template <typename SpectrumType>
      void store(const String& filename, const SpectrumType& spectrum) const
      {
        std::cerr << "Store() for SpecArrayFile not implemented. Filename was: " << filename << ", spec of size " << spectrum.size() << "\n";
        throw Exception::NotImplemented (__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
  };
} // namespace OpenMS

#endif // OPENMS_FORMAT_SPECARRAYFILE_H

