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


#ifndef OPENMS_SYSTEM_FILE_H
#define OPENMS_SYSTEM_FILE_H

#include <vector>
#include <OpenMS/CONCEPT/Types.h>
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
			
			/// Return true if the file exists and is readable
			static bool readable(const String& file);

			/// Return true if the file is writable
			static bool writable(const String& file);

			/**
				@brief Looks up the location of the file @p filename
				
				The following locations are checked in this order:
				- the directories in @p directories
				- the directory contained in the environment variable $OPENMS_DATA_PATH
				- the 'share/OpenMS/' directory of the OpenMS install directory
				
				@exception FileNotFound is thrown, if the file is not found
			*/
			static String find(const String& filename, std::vector<String> directories = std::vector<String>());
			
			/**
				@brief Retrieves a list of files matching @p file_pattern in directory @p dir
				
				@return true => there are matching files
			*/
			static bool fileList(const String& dir, const String& file_pattern, std::vector<String>& output);

			/// Returns a string, consisting of date, time, hostname, process id, and a incrementing number.  This can be used for temporary files.
			static String getUniqueName();

			/**
				@brief Creates a sparse file @p filename of size @p filesize bytes
				
      	Creates a sparse* file @p filename (*requires Filesystem support!) of size @p filesize bytes using platform specific fileIO
      	The function is using 64-bit fileoffsets automatically (and is therefore independent of compiler flags)
			*/
      static bool createSparseFile(const String& filename, const Int64& filesize);

      /**
				@brief Extends a sparse file with handle @p hFile to size @p filesize bytes
			
				Extends a sparse* file with handle @p hFile (*requires Filesystem support!) to size @p filesize bytes using platform specific fileIO
      	The function is using 64-bit fileoffsets automatically (and is therefore independent of compiler flags)
			*/
			#ifdef OPENMS_WINDOWSPLATFORM
			static bool extendSparseFile(const HANDLE& hFile, const Int64& filesize);
			#else
			static bool extendSparseFile(const int& hFile, const Int64& filesize);
			#endif
						    
			/**
				@brief get a handle to a sparse file @p filename with size @p filesize to used for swap
				
				The file can be created (@p create) if necessary
				  
      	@return handle to a file (which is created if necessary)
      	
				@throws Exception::UnableToCreateFile or Exception::FileNotFound on failure to acquire the handle (to make cross platform error handling easy)
      	
				@note implementation is platform dependent, as handles in Windows are void* vs. int in Unix
			*/
      #ifdef OPENMS_WINDOWSPLATFORM
      static HANDLE getSwapFileHandle(const String& filename, const Int64& filesize, const bool& create);
      #else
      static    int getSwapFileHandle(const String& filename, const Int64& filesize, const bool& create);
      #endif
  	
			/**
				@brief close handle to a swap file
			*/
      #ifdef OPENMS_WINDOWSPLATFORM
      static void closeSwapFileHandle(const HANDLE & f_handle);
      #else
      static void closeSwapFileHandle(const int & f_handle);
      #endif
    
	};

}

#endif // OPENMS_SYSTEM_FILE_H
