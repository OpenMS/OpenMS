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
// $Maintainer: $
// $Authors: Marc Sturm $
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
START_SECTION((SourceFile()))
	ptr = new SourceFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~SourceFile()))
	delete ptr;
END_SECTION

START_SECTION((Real getFileSize() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getFileSize(),0);
END_SECTION

START_SECTION((void setFileSize(Real file_size)))
  SourceFile tmp;
	tmp.setFileSize(1.667f);
  TEST_REAL_SIMILAR(tmp.getFileSize(),1.667f);
END_SECTION

START_SECTION((const String& getFileType() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getFileType(), "");
END_SECTION

START_SECTION((void setFileType(const String& file_type)))
  SourceFile tmp;
	tmp.setFileType("PEAKDATA");
  TEST_EQUAL(tmp.getFileType(), "PEAKDATA");
END_SECTION

START_SECTION((const String& getNameOfFile() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getNameOfFile(),"");
END_SECTION

START_SECTION((void setNameOfFile(const String& name_of_file)))
  SourceFile tmp;
  tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
  TEST_EQUAL(tmp.getNameOfFile(),"The White Stripes - Ball and Biscuit");
END_SECTION

START_SECTION((const String& getPathToFile() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getPathToFile(),"");
END_SECTION

START_SECTION((void setPathToFile(const String& path_path_to_file)))
  SourceFile tmp;
  tmp.setPathToFile("/misc/sturm/mp3/");
  TEST_EQUAL(tmp.getPathToFile(),"/misc/sturm/mp3/");
END_SECTION

START_SECTION((const String& getChecksum() const))
  SourceFile tmp;
  TEST_EQUAL(tmp.getChecksum(), "");
END_SECTION

START_SECTION(ChecksumType getChecksumType() const)
  SourceFile tmp;
  TEST_EQUAL(tmp.getChecksumType(), SourceFile::UNKNOWN_CHECKSUM);
END_SECTION

START_SECTION((void setChecksum(const String& checksum, ChecksumType type)))
  SourceFile tmp;
  tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12",SourceFile::SHA1);
  TEST_EQUAL(tmp.getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
  TEST_EQUAL(tmp.getChecksumType(), SourceFile::SHA1);
END_SECTION

START_SECTION((const String& getNativeIDType() const))
	SourceFile tmp;
	TEST_STRING_EQUAL(tmp.getNativeIDType(), "");
END_SECTION

START_SECTION((void setNativeIDType(const String& type)))
  SourceFile tmp;
  tmp.setNativeIDType("bla");
  TEST_STRING_EQUAL(tmp.getNativeIDType(), "bla");
END_SECTION

START_SECTION((SourceFile(const SourceFile& source)))
	SourceFile tmp;
	tmp.setFileType("CALIBRATIONINFO");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::MD5);
	tmp.setMetaValue("bla",4.0);
	
	SourceFile tmp2(tmp);
	TEST_EQUAL(tmp2.getFileType(), "CALIBRATIONINFO");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
	TEST_EQUAL(tmp2.getChecksum(), "2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::MD5);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("bla"), 4.0);
END_SECTION

START_SECTION((SourceFile& operator= (const SourceFile& source)))
	SourceFile tmp;
	tmp.setFileType("PUBLICATION");
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	tmp.setPathToFile("/misc/sturm/mp3/");
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::MD5);
	tmp.setMetaValue("bla",4.0);
	
	//normal assignment
	SourceFile tmp2;
	tmp2 = tmp;
	TEST_EQUAL(tmp2.getFileType(),"PUBLICATION");
	TEST_EQUAL(tmp2.getNameOfFile(),"The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp2.getPathToFile(),"/misc/sturm/mp3/");
	TEST_EQUAL(tmp2.getChecksum(),"2fd4e1c67a2d28fced849ee1bb76e7391b93eb12");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::MD5);
	TEST_REAL_SIMILAR(tmp2.getMetaValue("bla"), 4.0);
	
	//assignment of empty object
	tmp2 = SourceFile();
	TEST_EQUAL(tmp2.getFileType(), "");
	TEST_EQUAL(tmp2.getNameOfFile(),"");
	TEST_EQUAL(tmp2.getPathToFile(),"");
	TEST_EQUAL(tmp2.getChecksum(),"");
	TEST_EQUAL(tmp2.getChecksumType(), SourceFile::UNKNOWN_CHECKSUM);
	TEST_EQUAL(tmp2.metaValueExists("bla"), false);
END_SECTION

START_SECTION((bool operator== (const SourceFile& rhs) const))
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setFileType("PARAMETERSFILE");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setChecksum("", SourceFile::MD5);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("bla",4.0);
	TEST_EQUAL(tmp==tmp2, false);	
	
	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION((bool operator!= (const SourceFile& rhs) const))
	SourceFile tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setFileType("MISC");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setNameOfFile("The White Stripes - Ball and Biscuit");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setChecksum("2fd4e1c67a2d28fced849ee1bb76e7391b93eb12", SourceFile::UNKNOWN_CHECKSUM);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("bla",4.0);
	TEST_EQUAL(tmp!=tmp2, true);	

	tmp2 = tmp;	
	tmp.setPathToFile("/misc/sturm/mp3/");
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



