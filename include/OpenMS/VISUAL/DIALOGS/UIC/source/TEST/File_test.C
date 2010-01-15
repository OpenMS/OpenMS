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

/////////////////////////////////////////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id: File_test.C 6114 2009-10-13 21:43:07Z marc_sturm $")

/////////////////////////////////////////////////////////////

START_SECTION((static bool exists(const String &file)))
	TEST_EQUAL(File::exists("does_not_exists.txt"), false)
	TEST_EQUAL(File::exists(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), true)
END_SECTION

START_SECTION((static bool empty(const String &file)))
	TEST_EQUAL(File::empty("does_not_exists.txt"), true)
	TEST_EQUAL(File::empty(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)	
	TEST_EQUAL(File::empty(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), false)
END_SECTION

START_SECTION((static bool remove(const String &file)))
	//deleting non-existing file
	TEST_EQUAL(File::remove("does_not_exists.txt"), true)
	
	//deleting existing file
	String filename;
	NEW_TMP_FILE(filename);
	ofstream os;
	os.open (filename.c_str(), ofstream::out);
	os << "File_test dummy file to delete" << endl;
	os.close();
	TEST_EQUAL(File::remove(filename), true)	
END_SECTION

START_SECTION((static bool readable(const String &file)))
	TEST_EQUAL(File::readable("does_not_exists.txt"), false)
	TEST_EQUAL(File::readable(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)	
	TEST_EQUAL(File::readable(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), true)
END_SECTION

START_SECTION((static bool writable(const String &file)))
	TEST_EQUAL(File::writable("/this/file/cannot/be/written.txt"), false)
	TEST_EQUAL(File::writable(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)
	TEST_EQUAL(File::writable(OPENMS_GET_TEST_DATA_PATH("File_test_imaginary.txt")), true)
		
	String filename;
	NEW_TMP_FILE(filename);
	TEST_EQUAL(File::writable(filename), true)
END_SECTION

START_SECTION((static String find(const String &filename, StringList directories=StringList())))
	TEST_EXCEPTION(Exception::FileNotFound,File::find("File.h"))

	TEST_NOT_EQUAL(File::find("OpenMS_DB.sql"),"");
END_SECTION

START_SECTION((static String absolutePath(const String &file)))
	NOT_TESTABLE
END_SECTION

START_SECTION((static String path(const String &file)))
	NOT_TESTABLE
END_SECTION

START_SECTION((static String basename(const String &file)))
	TEST_EQUAL(File::basename("/souce/config/bla/bluff.h"),"bluff.h");
END_SECTION

START_SECTION((static bool fileList(const String &dir, const String &file_pattern, StringList &output, bool full_path=false)))
StringList vec;
TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "*.bliblaluff", vec), false);
TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "File_test_text.txt", vec), true);
TEST_EQUAL(vec[0], "File_test_text.txt");
TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "File_test_text.txt", vec, true), true);
// can't use "path + separator + filename", because sep. might be "/" or "\\"
TEST_EQUAL(vec[0].hasPrefix(OPENMS_GET_TEST_DATA_PATH("")), true);
TEST_EQUAL(vec[0].hasSuffix("File_test_text.txt"), true);
END_SECTION

START_SECTION((static String getUniqueName()))
	String unique_name = File::getUniqueName();
	
	// test if the string consists of three parts
	StringList split;
	unique_name.split('_', split);
	TEST_EQUAL(split.size() >= 4, true) // if name of machine also contains '_' ...
END_SECTION

START_SECTION((static String getOpenMSDataPath()))
	NOT_TESTABLE
END_SECTION

START_SECTION((static String removeExtension(const String& file)))
	TEST_STRING_EQUAL(File::removeExtension(""),"")
	TEST_STRING_EQUAL(File::removeExtension("/home/doe/file"),"/home/doe/file")
	TEST_STRING_EQUAL(File::removeExtension("/home/doe/file.txt"),"/home/doe/file")
	TEST_STRING_EQUAL(File::removeExtension("/home/doe/file.txt.tgz"),"/home/doe/file.txt")
END_SECTION

START_SECTION((static bool isDirectory(const String& path)))
	TEST_EQUAL(File::isDirectory(""),false)
	TEST_EQUAL(File::isDirectory("."),true)
	TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("")),true)
	TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("does_not_exists.txt")),false)
	TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")),false)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
