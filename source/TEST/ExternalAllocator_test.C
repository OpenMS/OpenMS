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
#include <OpenMS/SYSTEM/ExternalAllocator.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExternalAllocator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExternalAllocator<int>* ptr = 0;
CHECK(ExternalAllocator())
{
	ptr = new ExternalAllocator<int>();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~ExternalAllocator())
{
	delete ptr;
}
RESULT

CHECK((pointer address(reference value) const ))
{
	int i = 123;
	ExternalAllocator<int> extalloc;
  TEST_EQUAL(extalloc.address(i), &i);
}
RESULT

CHECK((const_pointer address(const_reference value) const))
{
	const int i = 123;
	ExternalAllocator<int> extalloc;
  TEST_EQUAL(extalloc.address(i), &i);
}
RESULT

CHECK((ExternalAllocator(const String &filename=File::getUniqueName(), const Offset64Int &filesize=1)))
{
	// this should work
  ExternalAllocator<int> extalloc;
	
	// this should NOT work
	TEST_EXCEPTION(Exception::UnableToCreateFile, ExternalAllocator<int> extalloc("this/file/does/not/exist", 10000) )
}
RESULT

CHECK((ExternalAllocator(const ExternalAllocator &rhs)))
{
  ExternalAllocator<int> extalloc;
	ExternalAllocator<int> extalloc2(extalloc);
	NOT_TESTABLE
}
RESULT

CHECK((template <class U> ExternalAllocator(const ExternalAllocator< U > &rhs)))
{
  ExternalAllocator<double> extalloc;
	ExternalAllocator<int> extalloc2(extalloc);
	NOT_TESTABLE
}
RESULT

CHECK((size_type max_size() const))
{
	String filename;
	NEW_TMP_FILE(filename);
  ExternalAllocator<int> extalloc(filename, 10000);
	TEST_EQUAL(extalloc.max_size(), 10000/sizeof(int))
}
RESULT


String filename;
NEW_TMP_FILE(filename);
ExternalAllocator<int> extalloc(filename, 10000);
ExternalAllocator<int>::pointer p;

CHECK((pointer allocate(size_type num, const void *=0)))
{
  p = extalloc.allocate(4);
	TEST_EQUAL(extalloc.getMappingSize(), (Offset64Int) MemoryMap::OpenMS_getFileBlocksize())
}
RESULT

CHECK((void construct(pointer p, const T &value)))
{
	extalloc.construct(p  , 123456);
	extalloc.construct(p+1, 23456);
	extalloc.construct(p+2, 3456);
	extalloc.construct(p+3, 456);
	// now check if it worked
	TEST_EQUAL(*p    , 123456)
	TEST_EQUAL(*(p+1), 23456)
	TEST_EQUAL(*(p+2), 3456)
	TEST_EQUAL(*(p+3), 456)
  
}
RESULT

CHECK((void destroy(pointer)))
{
	extalloc.destroy(p);
	extalloc.destroy(p+1);
	extalloc.destroy(p+2);
	extalloc.destroy(p+3);
	NOT_TESTABLE	
}
RESULT

CHECK(Offset64Int getMappingSize())
{
	TEST_EQUAL(extalloc.getMappingSize(), (Offset64Int) MemoryMap::OpenMS_getFileBlocksize())
}
RESULT

CHECK((void deallocate(pointer p, size_type num)))
{
  extalloc.deallocate(p, 4);
	TEST_EQUAL(extalloc.getMappingSize(), 0);
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



