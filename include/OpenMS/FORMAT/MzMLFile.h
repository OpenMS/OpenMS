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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_MZMLFILE_H
#define OPENMS_FORMAT_MZMLFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/MzMLHandler.h>
#include <OpenMS/FORMAT/PeakFileOptions.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/METADATA/DocumentIdentifier.h>

namespace OpenMS
{
	/**
		@brief File adapter for MzML files

		This implementation does currently not support the whole functionality of MzML.
		Some minor features are still missing:
			- chromatograms
			
		@todo Implement chromatograms (Andreas)
		
		@ingroup FileIO
	*/
	class OPENMS_DLLAPI MzMLFile
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

				//set DocumentIdentifier
				map.setLoadedFileType(filename);
				map.setLoadedFilePath(filename);

				Internal::MzMLHandler<MapType> handler(map,filename,schema_version_,*this);
				handler.setOptions(options_);
				//handler can throw parse error and other errors - catch those here - they are the cause for a parse error - report accordingly
				try
				{
				parse_(filename, &handler);
			}
				catch (Exception::BaseException& e)
				{
					std::string expr;
					expr.append(e.getFile());
					expr.append("@");
					std::stringstream ss;
					ss << e.getLine(); // we need c++11!! maybe in 2012?
					expr.append(ss.str());
					expr.append("-");
					expr.append(e.getFunction());
					std::string mess = "- due to that error of type ";
					mess.append(e.getName());
					throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,expr, mess);
				}

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

			/**
				@brief Checks if a file validates against the XML schema.

		  	@exception Exception::FileNotFound is thrown if the file cannot be found.
			*/
			bool isValid(const String& filename, std::ostream& os = std::cerr);

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

			/// Location of indexed mzML schema
			String indexed_schema_location_;
	};

} // namespace OpenMS

#endif // OPENMS_FOMAT_MZMLFILE_H
