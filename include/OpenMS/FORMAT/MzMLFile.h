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

#ifndef OPENMS_FORMAT_MZMLFILE_H
#define OPENMS_FORMAT_MZMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>

namespace OpenMS
{
	/**
		@brief File adapter for MzML files
		
		This is a beta implementation that does not support all functionality of MzML:
		- Currently only reading is supported
  	- Missing features are:
  	  - chromatograms
  	  - zlib compression of base64 data
  	  - base64 integer data
  	  - base64 16 bit data
		- Meta information that does not fit into the %OpenMS object model is ignored.
		
		If a critical error occurs, Exception::NotImplemented is thrown.
		
		@todo Implement mzML semantic validator (Hiwi)
		@todo Add tests for PeakFileOptions (Hiwi)
		@todo Implement and use Base64 integers and 16 bit (Hiwi)
		@todo Implement and use zlib support (Hiwi)
		@todo Implement all extensions to our data model that are necessary for mzML (Hiwi, Marc)
		@todo Implement 'id' attribute - stored in ExperimentalSettings (Clemens, Chris, Marc)
		
		@ingroup FileIO
	*/
	class MzMLFile 
		:	public Internal::XMLFile,
			public ProgressLogger
	{
		public:
			///Default constructor
			MzMLFile();
			///Destructor
			~MzMLFile();
			
      /// Mutable access to the options for loading/storing 
      PeakFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing 
      const PeakFileOptions& getOptions() const;

			/**
				@brief Loads a map from a MzML file.

				@p map has to be a MSExperiment or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			template <typename MapType>
			void load(const String& filename, MapType& map)
			{
				map.reset();
				
				Internal::MzMLHandler<MapType> handler(map,filename,schema_version_,*this);
				handler.setOptions(options_);
				parse_(filename, &handler);
			}

			/**
				@brief Stores a map in a MzML file.

				@p map has to be a MSExperiment or have the same interface.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			template <typename MapType>
			void store(const String& filename, const MapType& map) const
			{
				Internal::MzMLHandler<MapType> handler(map,filename,schema_version_,*this);
				handler.setOptions(options_);
				save_(filename, &handler);
			}
			
		private:
			
			/// Options for loading / storing
			PeakFileOptions options_;
	};

} // namespace OpenMS

#endif // OPENMS_FOMAT_MZMLFILE_H
