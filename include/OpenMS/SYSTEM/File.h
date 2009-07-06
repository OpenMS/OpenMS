// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
// --------------------------------------------------------------------------


#ifndef OPENMS_SYSTEM_FILE_H
#define OPENMS_SYSTEM_FILE_H

#include <OpenMS/DATASTRUCTURES/StringList.h>
#include <OpenMS/config.h>

#ifdef OPENMS_WINDOWSPLATFORM  
	#undef   _WIN32_WINNT        // avoid warning
	#define  _WIN32_WINNT 0x0500
	#include <Windows.h>
#endif

namespace OpenMS
{
	class String;
	
	/**
		@brief Basic file handling operations.
		
		@ingroup System
	*/
	class OPENMS_DLLAPI File
	{
		public:
			
			/// Method used to test if a @p file exists.
			static bool exists(const String& file);
		
			/// Return true if the file does not exist or the file is empty
			static bool empty(const String& file);
		
			/**
				@brief Removes a file (if it exists). 
			
				@return Returns true if the file was successfully deleted (or if it did not exist).
			*/
			static bool remove(const String& file);

			/// Replaces the relative path in the argument with the absolute path.
			static String absolutePath(const String& file);

			/// Returns the basename of the file (without the path).
			static String basename(const String& file);

			/// Returns the path of the file (without the file name).
			static String path(const String& file);

			/**
				Returns the file name without the extension
				
				The extension is the suffix of the string upto and including the last dot.
				
				If no extension is found, the whole file name is returned
			*/
			static String removeExtension(const String& file);
			
			/// Return true if the file exists and is readable
			static bool readable(const String& file);

			/// Return true if the file is writable
			static bool writable(const String& file);

			/// Return true if the given path specifies a directory
			static bool isDirectory(const String& path);

			/**
				@brief Looks up the location of the file @p filename
				
				The following locations are checked in this order:
				- the directories in @p directories
				- the directory contained in the environment variable $OPENMS_DATA_PATH
				- the 'share/OpenMS/' directory of the OpenMS install directory
				
				@exception FileNotFound is thrown, if the file is not found
			*/
			static String find(const String& filename, StringList directories = StringList());
			
			/**
				@brief Retrieves a list of files matching @p file_pattern in directory @p dir (returns filenames without paths unless @p full_path is true)
				
				@return true => there are matching files
			*/
			static bool fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path = false);

			/// Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
			static String getUniqueName();
			
			/// Returns the OpenMS data path (environment variable overwrites the default installation path)
			static String getOpenMSDataPath();
	};

}

#endif // OPENMS_SYSTEM_FILE_H
