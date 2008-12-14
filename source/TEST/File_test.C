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

#include <OpenMS/CONCEPT/ClassTest.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static bool exists(const String &file)))
	TEST_EQUAL(File::exists("does_not_exists.txt"), false)
	TEST_EQUAL(File::exists("data/File_test_text.txt"), true)
END_SECTION

START_SECTION((static bool empty(const String &file)))
	TEST_EQUAL(File::empty("does_not_exists.txt"), true)
	TEST_EQUAL(File::empty("data/File_test_empty.txt"), true)	
	TEST_EQUAL(File::empty("data/File_test_text.txt"), false)
END_SECTION

START_SECTION((static bool remove(const String &file)))
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
END_SECTION

START_SECTION((static bool readable(const String &file)))
	TEST_EQUAL(File::readable("does_not_exists.txt"), false)
	TEST_EQUAL(File::readable("data/File_test_empty.txt"), true)	
	TEST_EQUAL(File::readable("data/File_test_text.txt"), true)
END_SECTION

START_SECTION((static bool writable(const String &file)))
	TEST_EQUAL(File::writable("/this/file/cannot/be/written/does_not_exists.txt"), false)
	String filename;
	NEW_TMP_FILE(filename);
	TEST_EQUAL(File::writable(filename), true)
END_SECTION

START_SECTION((static String find(const String &filename, std::vector< String > directories=std::vector< String >())))
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

START_SECTION((static bool fileList(const String &dir, const String &file_pattern, std::vector< String > &output)))
	vector<String> vec;
	TEST_EQUAL(File::fileList("data/","*.bliblaluff",vec),false);
END_SECTION

START_SECTION((static String getUniqueName()))
	String unique_name = File::getUniqueName();
	
	// test if the string consists of three parts
	vector<String> split;
	unique_name.split('_', split);
	TEST_EQUAL(split.size() >= 4, true) // if name of machine also contains '_' ...
END_SECTION

START_SECTION((static bool createSparseFile(const String &filename, const Int64 &filesize)))
	String filename;
	NEW_TMP_FILE(filename);
  
  //create sparse file
	TEST_EQUAL(File::createSparseFile(filename, 10000), true)	
  
  //delete file
  TEST_EQUAL(File::remove(filename), true)	
END_SECTION

	
#ifdef OPENMS_WINDOWSPLATFORM
START_SECTION((static int getSwapFileHandle(const String &filename, const Int64 &filesize, const bool &create)))
#else
START_SECTION((static int getSwapFileHandle(const String &filename, const Int64 &filesize, const bool &create)))
#endif
	String filename;
	NEW_TMP_FILE(filename);
  
  //create sparse file with 300GB
	File::closeSwapFileHandle(File::getSwapFileHandle(filename, 10000, true));
  STATUS("filename:" + filename);
  //delete file (this will fail if handle is not closed on Windows)
  //TEST_EQUAL(File::remove(filename), true)	
  
  // test failure if create-flag is false and file does not exist
  TEST_EXCEPTION(Exception::FileNotFound, File::getSwapFileHandle("this/file/does/not/exist", 1000, false) )
  
END_SECTION


#ifdef OPENMS_WINDOWSPLATFORM
START_SECTION((static bool extendSparseFile(const int &hFile, const Int64 &filesize)))
#else
START_SECTION((static bool extendSparseFile(const int &hFile, const Int64 &filesize)))
#endif
	String filename;
	NEW_TMP_FILE(filename);
  
  //create sparse file with 200GB
	#ifdef OPENMS_WINDOWSPLATFORM
	HANDLE hFile = File::getSwapFileHandle(filename, 0, true);
	#else
	int hFile = File::getSwapFileHandle(filename, 0, true);
	#endif
	
	TEST_EQUAL(File::extendSparseFile(hFile, 999), true);
	File::closeSwapFileHandle(hFile);
	
  //delete file
  //TEST_EQUAL(File::remove(filename), true)	

END_SECTION


#ifdef OPENMS_WINDOWSPLATFORM
START_SECTION((static void closeSwapFileHandle(const int &f_handle)))
#else
START_SECTION((static void closeSwapFileHandle(const int &f_handle)))
#endif 

  String filename;
	NEW_TMP_FILE(filename);
  
  File::closeSwapFileHandle(File::getSwapFileHandle(filename, 1000, true));

  //delete file
  //TEST_EQUAL(File::remove(filename), true)	
  
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
