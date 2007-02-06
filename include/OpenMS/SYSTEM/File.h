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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <string>
#include <vector>

#ifndef OPENMS_SYSTEM_FILE_H
#define OPENMS_SYSTEM_FILE_H

namespace OpenMS
{
	class String;
	
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
		
			/**
				@brief Removes a file (if it exists). 
			
				@return Returns true if the file was successfully deleted (or if it did not exist).
			*/
			static bool remove(const std::string& file);

			/// Replaces the relative path in the argument with the absolute path.
			static void absolutePath(std::string& file);

			/// Returns the basename of the file (without the path).
			static String basename(const std::string& file);

			/// Returns the path of the file (without the file name).
			static String path(const std::string& file);
			
			/// Return true if the file exists and is readable
			static bool readable(const std::string& file);

			/// Return true if the file is writable
			static bool writable(const std::string& file);
			
			/**
				@brief Looks up the location of @p filename
				
				First the directories in @p directories are cheched, 
				then the 'data' directory of the environment variable $OPENMS_PATH is checked,
				at last the 'data' directory of the OpenMS built directory is checked.
				
				If the file is not found there, an empty string is returned.
			*/
			static String find(const String& filename, std::vector<String> directories = std::vector<String>());
			
			/**
				@brief Method for getting list of files matching @p file_pattern in directory @p dir
				
				@return If there are matching files
			*/
			static bool fileList(const std::string& dir, const std::string& file_pattern, std::vector<String>& output);
	};

}

#endif // OPENMS_SYSTEM_FILE_H
