// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
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
TextFile* nullPointer = 0;

START_SECTION((TextFile()))
	ptr = new TextFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~TextFile()))
	delete ptr;
END_SECTION

START_SECTION((void load(const String& filename, bool trim_lines=false, Int first_n=-1) ))
	TextFile file;
	
	TEST_EXCEPTION(Exception::FileNotFound, file.load("FileDoesNotExist.txt"))	
	
	file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"));
	TEST_EQUAL(file.size(), 11)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	
	//trimmed
	file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true);
	TEST_EQUAL(file.size(), 11)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	TEST_EQUAL(file[5].trim() == "space_line", true)
	TEST_EQUAL(file[6].trim() == "tab_line", true)
	TEST_EQUAL(file[7].trim() == "back_space_line", true)
	TEST_EQUAL(file[8].trim() == "back_tab_line", true)
	
	//only first few
	file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,1);
	TEST_EQUAL(file.size(), 1)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	
	file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,3);
	TEST_EQUAL(file.size(), 3)
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[1].trim() == "", true)
	TEST_EQUAL(file[2].trim() == "", true)
	
	file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,4);
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
	TextFile file(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"));
	TEST_EQUAL(file[0].trim() == "first_line", true)
	TEST_EQUAL(file[3].trim() == "middle_line", true)
	TEST_EQUAL(file[10].trim() == "last_line", true)
	TEST_EQUAL(file.size(), 11)
	
	TextFile file2(OPENMS_GET_TEST_DATA_PATH("TextFile_test_empty_infile.txt"));
	TEST_EQUAL(file2.size(), 0)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
