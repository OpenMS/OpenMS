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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/METADATA/SourceFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SourceFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SourceFile* ptr = 0;
CHECK(SourceFile())
	ptr = new SourceFile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~SourceFile())
	delete ptr;
RESULT

CHECK(const String& getFileType() const)
  SourceFile tmp;
  TEST_EQUAL(tmp.getFileType(),"");
RESULT

CHECK(void setFileType(const String& file_type))
  SourceFile tmp;
	tmp.setFileType("mp3");
  TEST_EQUAL(tmp.getFileType(),"mp3");
RESULT

CHECK(const String& getNameOfFile() const)
  SourceFile tmp;
  TEST_EQUAL(tmp.getNameOfFile(),"");
RESULT

CHECK(void setNameOfFile(const String& name_of_file))
  SourceFile tmp;
  tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
  TEST_EQUAL(tmp.getNameOfFile(),"The White Stripes - Ball and Biscuit");
RESULT

CHECK(const String& getPathToFile() const)
  SourceFile tmp;
  TEST_EQUAL(tmp.getPathToFile(),"");
RESULT

CHECK(void setPathToFile(const String& path_path_to_file))
  SourceFile tmp;
  tmp.setPathToFile("/misc/sturm/mp3/");
  TEST_EQUAL(tmp.getPathToFile(),"/misc/sturm/mp3/");
RESULT

CHECK(SourceFile(const SourceFile& source))
	SourceFile tmp;
	tmp.setFileType("mp3");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	
	SourceFile tmp2(tmp);
	TEST_EQUAL(tmp2.getFileType(),"mp3");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
RESULT

CHECK(SourceFile& operator= (const SourceFile& source))
	SourceFile tmp;
	tmp.setFileType("mp3");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	
	//normal assignment
	SourceFile tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getFileType(),"mp3");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
	
	//assignment of empty object
	tmp2 = SourceFile();
	TEST_EQUAL(tmp2.getFileType(),"");
	TEST_EQUAL(tmp2.getNameOfFile(),"");
	TEST_EQUAL(tmp2.getPathToFile(),"");
RESULT

CHECK(bool operator!= (const SourceFile& rhs) const)
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setFileType("mp3");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp==tmp2, false);
RESULT

CHECK(bool operator== (const SourceFile& rhs) const)
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setFileType("mp3");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp!=tmp2, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



