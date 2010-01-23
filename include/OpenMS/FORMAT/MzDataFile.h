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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZDATAFILE_H
#define OPENMS_FORMAT_MZDATAFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/FORMAT/HANDLERS/MzDataHandler.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

namespace OpenMS
{
	/**
		@brief File adapter for MzData files
		
		@ingroup FileIO
	*/
	class OPENMS_DLLAPI MzDataFile
		:	public Internal::XMLFile,
			public ProgressLogger
	{

		public:

			///Default constructor
			MzDataFile();
			///Destructor
			~MzDataFile();

      /// Mutable access to the options for loading/storing
      PeakFileOptions& getOptions();

      /// Non-mutable access to the options for loading/storing
      const PeakFileOptions& getOptions() const;

			/**
				@brief Loads a map from a MzData file.

				@p map has to be a MSExperiment or have the same interface.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
			template <typename MapType>
			void load(const String& filename, MapType& map)
			{
				map.reset();

				//set DocumentIdentifier
				map.setLoadedFileType(filename);
				map.setLoadedFilePath(filename);

				Internal::MzDataHandler<MapType> handler(map,filename,schema_version_,*this);
				handler.setOptions(options_);
				parse_(filename, &handler);
			}

			/**
				@brief Stores a map in a MzData file.

				@p map has to be a MSExperiment or have the same interface.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			template <typename MapType>
			void store(const String& filename, const MapType& map) const
			{
				Internal::MzDataHandler<MapType> handler(map,filename,schema_version_,*this);
				handler.setOptions(options_);
				save_(filename, &handler);
			}

			/**
				@brief Checks if a file is valid with respect to the mapping file and the controlled vocabulary.

				@param filename File name of the file to be checked.
				@param errors Errors during the validation are returned in this output parameter.
				@param warnings Warnings during the validation are returned in this output parameter.

				@exception Exception::FileNotFound is thrown if the file could not be opened
			*/
			bool isSemanticallyValid(const String& filename, StringList& errors, StringList& warnings);

		private:

			/// Options for loading / storing
			PeakFileOptions options_;
	};

} // namespace OpenMS

#endif // OPENMS_FOMAT_MZXMLFILE_H
