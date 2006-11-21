// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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

#include <string>

#ifndef OPENMS_SYSTEM_FILE_H
#define OPENMS_SYSTEM_FILE_H

namespace OpenMS
{
	/**
		@brief Basic file handling operations.
		
		@ingroup System
	*/
	class File
	{
		public:
			
			/// Method used to test if a @p file exists.
			static bool exists(const std::string& file);
		
			/// Return true if the file does not exist of the file is empty
			static bool empty(const std::string& file);
		
			/// Removes a file (if it exists). Returns true if the file was successfully deleted (or if it did not exist).
			static bool remove(const std::string& file);

			/// Replaces the relative path in the argument with the absolute path.
			static void absolutePath(std::string& file);
			
			/// Return true if the file exists and is readable
			static bool readable(const std::string& file);

			/// Return true if the file is writable
			static bool writable(const std::string& file);
	};

}

#endif // OPENMS_SYSTEM_FILE_H
