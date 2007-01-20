// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Thomas Kadauke $
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_SCHEMAFILE_H
#define OPENMS_FORMAT_SCHEMAFILE_H

// OpenMS includes
#include <OpenMS/CONCEPT/Exception.h>

namespace OpenMS
{
	// forward declarations
	class String;
	
	namespace Internal
	{
		class SchemaHandler;
	
		///Base class for loading/storing XML files that have a handler derived from SchemaHandler.
		class SchemaFile
		{
			public:
				///Default constructor
				SchemaFile();
				///Destructor
				~SchemaFile();
	
			protected:
				/**
					Parses the XML file given by @p filename using the handler given by @p handler.
				*/
				void parse_(const String& filename, SchemaHandler* handler) throw (Exception::FileNotFound, Exception::ParseError);
	
				/**
					Stores the contents of the XML handler given by @p handler in the file given by @p filename.
				*/
				void save_(const String& filename, SchemaHandler* handler) const throw (Exception::UnableToCreateFile);
		};
	}
} // namespace OpenMS

#endif // OPENMS_FOMAT_SchemaFile_H
