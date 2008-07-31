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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/SYSTEM/MemoryMap.h>
///////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>


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

CHECK((static std::size_t OpenMS_getFileBlocksize(void)))
{
  std::size_t page;
	page = MemoryMap::OpenMS_getFileBlocksize();
	// architecture dependent, usually its 4KB on unix and 64KB on windows
	NOT_TESTABLE
}
RESULT

void* mapping = 0;

#ifdef OPENMS_WINDOWSPLATFORM
HANDLE h;
CHECK([EXTRA](static void* OpenMS_mmap (const std::size_t& size, const HANDLE& handle, const Offset64Int& file_offset)))
#else
long h;
CHECK([EXTRA](static void* OpenMS_mmap (const std::size_t& size, const int& fileHandle, const Offset64Int& file_offset)))
#endif
{
	String filename;
	NEW_TMP_FILE(filename);
  h = File::getSwapFileHandle(filename, 1000LL, true);
  
	mapping = MemoryMap::OpenMS_mmap(1000, h, 0);
  
	TEST_NOT_EQUAL(mapping, 0);
}
RESULT

CHECK((static int OpenMS_unmap (void* p, const std::size_t& bytes)))
{
  int r = MemoryMap::OpenMS_unmap(mapping, 1000);
	TEST_NOT_EQUAL(r, OPENMS_MUNMAP_FAILURE);
	
	File::closeSwapFileHandle(h);
}
RESULT



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



