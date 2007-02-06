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

#include <OpenMS/CONCEPT/ClassTest.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

CHECK((bool exists(const std::string &file)))
	TEST_EQUAL(File::exists("does_not_exists.txt"), false)
	TEST_EQUAL(File::exists("data/File_test_text.txt"), true)
RESULT

CHECK((bool empty(const std::string &file)))
	TEST_EQUAL(File::empty("does_not_exists.txt"), true)
	TEST_EQUAL(File::empty("data/File_test_empty.txt"), true)	
	TEST_EQUAL(File::empty("data/File_test_text.txt"), false)
RESULT

CHECK((bool remove(const std::string &file)))
	//deleteing non-existing file
	TEST_EQUAL(File::remove("does_not_exists.txt"), true)
	
	//deleteing existing file
	String filename;
	NEW_TMP_FILE(filename);
	ofstream os;
	os.open (filename.c_str(), ofstream::out);
	os << "File_test dummy file to delete" << endl;
	os.close();
	TEST_EQUAL(File::remove(filename), true)	
RESULT

CHECK((bool readable(const std::string &file)))
	TEST_EQUAL(File::readable("does_not_exists.txt"), false)
	TEST_EQUAL(File::readable("data/File_test_empty.txt"), true)	
	TEST_EQUAL(File::readable("data/File_test_text.txt"), true)
RESULT

CHECK((bool writable(const std::string &file)))
	TEST_EQUAL(File::writable("/this/file/cannot/be/written/does_not_exists.txt"), false)
	String filename;
	NEW_TMP_FILE(filename);
	TEST_EQUAL(File::writable(filename), true)
RESULT

CHECK((String find(const String &filename, std::vector< String > directories=std::vector< String >())))
	TEST_EQUAL(File::find("File.h"),"");
	vector<String> vec;
	vec.push_back(OPENMS_PATH"/include/OpenMS/SYSTEM/");
	TEST_NOT_EQUAL(File::find("File.h",vec),"");
RESULT

CHECK((void absolutePath(std::string &file)))
	// not testable
RESULT

CHECK((String path(const std::string &file)))
	// not testable
RESULT

CHECK((String basename(const std::string &file)))
	TEST_EQUAL(File::basename("/souce/config/bla/bluff.h"),"bluff.h");
RESULT

CHECK((bool fileList(const std::string &dir, const std::string &file_pattern, std::vector< String > &output)))
	vector<String> vec;
	TEST_EQUAL(File::fileList("data/","*.bliblaluff",vec),false);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
