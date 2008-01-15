// -*- Mode: C++; tab-width: 2; -*-
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SYSTEM/MemoryMap.h>
///////////////////////////

#include <OpenMS/SYSTEM/File.h>


using namespace OpenMS;
using namespace std;

START_TEST(MemoryMap, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MemoryMap* ptr = 0;
CHECK(MemoryMap())
{
	ptr = new MemoryMap();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~MemoryMap())
{
	delete ptr;
}
RESULT

CHECK((static std::size_t OpenMS_getpagesize(void)))
{
  std::size_t page;
	page = MemoryMap::OpenMS_getpagesize();
	// architecture dependent, usually its 4KB
}
RESULT

void* mapping = 0;

#ifdef OPENMS_WINDOWSPLATFORM
HANDLE h;
CHECK((static void* OpenMS_mmap (std::size_t size, HANDLE handle, Offset64Int file_offset)))
#else
long h;
CHECK((static void* OpenMS_mmap(std::size_t size, long fileHandle, Offset64Int file_offset)))
#endif
{
	String filename;
	NEW_TMP_FILE(filename);
  h = File::getSwapFileHandle(filename, 1000LL, true);
  
	mapping = MemoryMap::OpenMS_mmap(1000, h, 0);
  
	
}
RESULT

CHECK((static int OpenMS_unmap(void *p, std::size_t bytes)))
{
  int r = MemoryMap::OpenMS_unmap(mapping, 1000);
	TEST_NOT_EQUAL(r, OPENMS_MUNMAP_FAILURE);
	
	File::closeSwapFileHandle(h);
}
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



