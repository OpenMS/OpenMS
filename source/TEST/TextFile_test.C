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

///////////////////////////

#include <OpenMS/FORMAT/TextFile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

TextFile* ptr = 0;

START_SECTION((TextFile()))
	ptr = new TextFile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((~TextFile()))
	delete ptr;
END_SECTION

START_SECTION((void load(const String& filename, bool trim_lines=false, Int first_n=-1) ))
	TextFile file;
	
	TEST_EXCEPTION(Exception::FileNotFound, file.load("FileDoesNotExist.txt"))	
	
	file.load("data/TextFile_test_infile.txt");
	TEST_EQUAL(file.size(), 11)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	
	//trimmed
	file.load("data/TextFile_test_infile.txt",true);
	TEST_EQUAL(file.size(), 11)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	TEST_EQUAL(file[5].trim() == "space_line", true)
	TEST_EQUAL(file[6].trim() == "tab_line", true)
	TEST_EQUAL(file[7].trim() == "back_space_line", true)
	TEST_EQUAL(file[8].trim() == "back_tab_line", true)
	
	//only first few
	file.load("data/TextFile_test_infile.txt",true,1);
	TEST_EQUAL(file.size(), 1)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	
	file.load("data/TextFile_test_infile.txt",true,3);
	TEST_EQUAL(file.size(), 3)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[1].trim() == "", true)
	TEST_EQUAL(file[2].trim() == "", true)
	
	file.load("data/TextFile_test_infile.txt",true,4);
	TEST_EQUAL(file.size(), 4)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[1].trim() == "", true)
	TEST_EQUAL(file[2].trim() == "", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
END_SECTION

START_SECTION((void store(const String& filename) ))
	TextFile file;

	TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/does/not/exist/FileDoesNotExist.txt"))	

	file.push_back("line1");
	file.push_back("line2\n");
	file.push_back("line3\r\n");
	String filename;
	NEW_TMP_FILE(filename);
	file.store(filename);
	file.load(filename);
	TEST_EQUAL(file[0] == "line1",true);
	TEST_EQUAL(file[1] == "line2",true);
	TEST_EQUAL(file[2] == "line3",true);
END_SECTION

START_SECTION((TextFile(const String& filename, bool trim_lines=false, Int first_n=-1) ))
	TextFile file("data/TextFile_test_infile.txt");
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	TEST_EQUAL(file.size(), 11)
	
	TextFile file2("data/TextFile_test_empty_infile.txt");
	TEST_EQUAL(file2.size(), 0)
END_SECTION

START_SECTION((Iterator search(const Iterator& start, const String& text, bool trim=false)))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.search(file.begin(),"first_line") == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search(file.begin(),"space_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin(),"last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()," ") == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"\t") == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin()+9,"\t") == (file.end()), true)
	
	//trim
	TEST_EQUAL(file.search(file.begin(),"first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line",true) == file.end(), true)
	
	//Try it on the same file (but trimmed)
	file.load("data/TextFile_test_infile.txt",true);

	TEST_EQUAL(file.search(file.begin(),"first_line") == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line") == file.end(), true)

	//trim
	TEST_EQUAL(file.search(file.begin(),"first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line",true) == file.end(), true)
END_SECTION

START_SECTION((Iterator search(const String& text, bool trim=false)))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.search("first_line") == file.begin(), true)
	TEST_EQUAL(file.search("middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search("space_line") == file.end(), true)
	TEST_EQUAL(file.search("tab_line") == file.end(), true)
	TEST_EQUAL(file.search("last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search("invented_line") == file.end(), true)
	TEST_EQUAL(file.search(" ") == (file.begin()+5), true)
	TEST_EQUAL(file.search("\t") == (file.begin()+6), true)
	
	//trim
	TEST_EQUAL(file.search("first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("invented_line",true) == file.end(), true)
	
	//Try it on the same file (but trimmed)
	file.load("data/TextFile_test_infile.txt",true);

	TEST_EQUAL(file.search("first_line") == file.begin(), true)
	TEST_EQUAL(file.search("middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search("invented_line") == file.end(), true)

	//trim
	TEST_EQUAL(file.search("first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("invented_line",true) == file.end(), true)
END_SECTION

START_SECTION((Iterator searchSuffix(const Iterator& start, const String& text, bool trim=false)))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.searchSuffix(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_space_line",true) == file.begin()+7, true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_tab_line",true) == file.begin()+8, true)
	TEST_EQUAL(file.searchSuffix(file.begin()+8,"back_space_line",true) == file.end(), true)
	
	TEST_EQUAL(file.searchSuffix(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_space_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_tab_line") == file.end(), true)
END_SECTION

START_SECTION((Iterator searchSuffix(const String& text, bool trim=false)))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.searchSuffix("invented_line",true) == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_space_line",true) == file.begin()+7, true)
	TEST_EQUAL(file.searchSuffix("back_tab_line",true) == file.begin()+8, true)
	
	TEST_EQUAL(file.searchSuffix("invented_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_space_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_tab_line") == file.end(), true)
END_SECTION

START_SECTION((ConstIterator search(const ConstIterator& start, const String& text, bool trim=false) const))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.search(file.begin(),"first_line") == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search(file.begin(),"space_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin(),"last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()," ") == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"\t") == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin()+9,"\t") == (file.end()), true)
	
	//trim
	TEST_EQUAL(file.search(file.begin(),"first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line",true) == file.end(), true)
	
	//Try it on the same file (but trimmed)
	file.load("data/TextFile_test_infile.txt",true);

	TEST_EQUAL(file.search(file.begin(),"first_line") == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line") == file.end(), true)

	//trim
	TEST_EQUAL(file.search(file.begin(),"first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search(file.begin(),"space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search(file.begin(),"tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.search(file.begin()+1,"first_line",true) == file.end(), true)
END_SECTION

START_SECTION((ConstIterator search(const String& text, bool trim=false) const))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.search("first_line") == file.begin(), true)
	TEST_EQUAL(file.search("middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search("space_line") == file.end(), true)
	TEST_EQUAL(file.search("tab_line") == file.end(), true)
	TEST_EQUAL(file.search("last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search("invented_line") == file.end(), true)
	TEST_EQUAL(file.search(" ") == (file.begin()+5), true)
	TEST_EQUAL(file.search("\t") == (file.begin()+6), true)
	
	//trim
	TEST_EQUAL(file.search("first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("invented_line",true) == file.end(), true)
	
	//Try it on the same file (but trimmed)
	file.load("data/TextFile_test_infile.txt",true);

	TEST_EQUAL(file.search("first_line") == file.begin(), true)
	TEST_EQUAL(file.search("middle_line") == (file.begin()+3), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("last_line") == (file.end()-1), true)
	TEST_EQUAL(file.search("invented_line") == file.end(), true)

	//trim
	TEST_EQUAL(file.search("first_line",true) == file.begin(), true)
	TEST_EQUAL(file.search("space_line",true) == (file.begin()+5), true)
	TEST_EQUAL(file.search("tab_line",true) == (file.begin()+6), true)
	TEST_EQUAL(file.search("invented_line",true) == file.end(), true)
END_SECTION

START_SECTION((ConstIterator searchSuffix(const ConstIterator& start, const String& text, bool trim=false) const))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.searchSuffix(file.begin(),"invented_line",true) == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_space_line",true) == file.begin()+7, true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_tab_line",true) == file.begin()+8, true)
	TEST_EQUAL(file.searchSuffix(file.begin()+8,"back_space_line",true) == file.end(), true)
	
	TEST_EQUAL(file.searchSuffix(file.begin(),"invented_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_space_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix(file.begin(),"back_tab_line") == file.end(), true)
END_SECTION

START_SECTION((ConstIterator searchSuffix(const String& text, bool trim=false) const))
	TextFile file("data/TextFile_test_infile.txt");
	
	TEST_EQUAL(file.searchSuffix("invented_line",true) == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_space_line",true) == file.begin()+7, true)
	TEST_EQUAL(file.searchSuffix("back_tab_line",true) == file.begin()+8, true)
	
	TEST_EQUAL(file.searchSuffix("invented_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_space_line") == file.end(), true)
	TEST_EQUAL(file.searchSuffix("back_tab_line") == file.end(), true)
END_SECTION

START_SECTION((String asString() const))
	TextFile tmp;
	TEST_EQUAL(tmp.asString(),"")
	tmp.push_back("1\n");
	tmp.push_back("2\n");
	tmp.push_back("3\n");
	TEST_EQUAL(tmp.asString(),"1\n2\n3\n")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
