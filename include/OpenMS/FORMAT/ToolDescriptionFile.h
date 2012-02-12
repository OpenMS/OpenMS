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
// $Maintainer:  Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H
#define OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H

#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/DATASTRUCTURES/ToolDescription.h>

namespace OpenMS
{
	/**
		@brief File adapter for ToolDescriptor files

		If a critical error occurs due to the missing functionality, Exception::NotImplemented is thrown.

		@ingroup FileIO
	*/
	class OPENMS_DLLAPI ToolDescriptionFile
		:	public Internal::XMLFile,
			public ProgressLogger
	{
		public:
			///Default constructor
			ToolDescriptionFile();
			///Destructor
			virtual ~ToolDescriptionFile();

			/**
				@brief Loads a map from a ToolDescriptor file.

				@exception Exception::FileNotFound is thrown if the file could not be opened
				@exception Exception::ParseError is thrown if an error occurs during parsing
			*/
      void load(const String& filename, std::vector <Internal::ToolDescription>& tds);

			/**
				@brief Stores a map in a ToolDescriptor file.

				@exception Exception::UnableToCreateFile is thrown if the file could not be created
			*/
			void store(const String& filename, const std::vector <Internal::ToolDescription>& tds) const;

		private:

	};

} // namespace OpenMS

#endif // OPENMS_FORMAT_TOOLDESCRIPTIONFILE_H
