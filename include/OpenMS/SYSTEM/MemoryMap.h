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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_SYSTEM_MEMORYMAP_H
#define OPENMS_SYSTEM_MEMORYMAP_H

#include <OpenMS/CONCEPT/Types.h>

#define OPENMS_MUNMAP_FAILURE (-1)

#ifdef OPENMS_WINDOWSPLATFORM
  #include <Windows.h>
  #define MAP_FAILED ((void *) -1) /* from mman.h (linux)      */
//	#define off64_t __int64  ##del this
#else
  #include <unistd.h>       //TODO: remove?
  #include <sys/mman.h>            /* WARNING: use #undef MAP_TYPE when done!! see bottom of file! */
#endif  




namespace OpenMS
{  


	/**
		@brief Cross platform memory mapping.
		
		@ingroup System
	*/
	class MemoryMap
	{
    public:
    // recreate the functions already offered by POSIX  
  
    static std::size_t OpenMS_getpagesize (void)
    {
      static std::size_t g_pagesize = 0;
      if (! g_pagesize) 
      {
        #ifdef OPENMS_WINDOWSPLATFORM
          SYSTEM_INFO system_info;
          GetSystemInfo (&system_info);
          g_pagesize = system_info.dwPageSize;
        #else
          g_pagesize = getpagesize();
        #endif  
      }
      return g_pagesize;    
    }
    
/*    long getregionsize (void) 
    {
        static long g_regionsize = 0;
        if (! g_regionsize) {
            SYSTEM_INFO system_info;
            GetSystemInfo (&system_info);
            g_regionsize = system_info.dwAllocationGranularity;
        }
        return g_regionsize;
    }*/
    
  
    #ifdef OPENMS_WINDOWSPLATFORM
    /* mmap for windows */
      static void* OpenMS_mmap (std::size_t size, HANDLE handle, Offset64Int file_offset)
      {
        LARGE_INTEGER iTmp;
        iTmp.QuadPart = file_offset;
        DWORD hi = iTmp.HighPart;
        DWORD lo = iTmp.LowPart;
    
        //std::cout << "in mmap: " << size << " " << file_offset << "\n";
    
        LPVOID map = MapViewOfFile(
                                    handle,									// A file handle
                                    FILE_MAP_ALL_ACCESS,		// A read/write view of the file is mapped
                                    hi,
                                    lo, 
                                    size
        );
      
        if (map == NULL) //FAILED
        { // convert to UNIX return value
          map = MAP_FAILED;	// ((void *) -1)
        }
        return map;
    
      }
    #else
      static void* OpenMS_mmap (std::size_t size, long fileHandle, Offset64Int file_offset)
      {
        
        void* map =  mmap64(0,
                            size,
                            PROT_READ | PROT_WRITE,
                            MAP_SHARED,
                            fileHandle,
                            file_offset
                          );
        return map;
      }
    #endif
  
  
    /// undo memory mapping at position @p p and size @p bytes 
    /// returns OPENMS_MUNMAP_FAILURE on failure
    static int OpenMS_unmap (void* p, std::size_t bytes)
    {
      #ifdef OPENMS_WINDOWSPLATFORM
        int result = UnmapViewOfFile(p);  
        if (result == 0) result = OPENMS_MUNMAP_FAILURE; // 0 indicates failure in Windows
      #else
        int result = munmap(p, bytes);
      #endif
      return result;
    }
  
  
  
  }; // end MemoryMap class


} //namespace OpenMS

#undef MAP_TYPE //this is a MACRO defined in mmap.h and conflicts with CGAL, which uses "MAP_TYPE" in an internal class as TEMPLATE parameter!  


#endif // OPENMS_SYSTEM_MEMORYMAP_H
