// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/TextFile.h>
#include <iostream>
#include <vector>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(String, "$Id$")

/////////////////////////////////////////////////////////////

TextFile* ptr = nullptr;
TextFile* nullPointer = nullptr;

START_SECTION((TextFile()))
  ptr = new TextFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((TextFile(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false) ))
  // just some basic stuff, since the C'Tor calls load() directly
  TextFile file(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"));
  TextFile::ConstIterator file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  file_it += 3;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)
  file_it += 7;
  TEST_EQUAL(String(*file_it).trim() == "last_line", true)
  TEST_EQUAL((file.end() - file.begin()), 11)

  TextFile file2(OPENMS_GET_TEST_DATA_PATH("TextFile_test_empty_infile.txt"));
  TEST_EQUAL((file2.end() - file2.begin()), 0)
END_SECTION

START_SECTION((~TextFile()))
  delete ptr;
END_SECTION


START_SECTION((void load(const String& filename, bool trim_lines = false, Int first_n = -1, bool skip_empty_lines = false) ))
  TextFile file;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("FileDoesNotExist.txt"))

  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"));
  TEST_EQUAL((file.end() - file.begin()), 11)
  TextFile::ConstIterator file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  file_it += 3;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)
  file_it += 7;
  TEST_EQUAL(String(*file_it).trim() == "last_line", true)

  //trimmed
  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true);
  TEST_EQUAL((file.end() - file.begin()), 11)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  file_it += 3;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)
  file_it += 2;
  TEST_EQUAL(String(*file_it).trim() == "space_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "tab_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "back_space_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "back_tab_line", true)
  file_it += 2;
  TEST_EQUAL(String(*file_it).trim() == "last_line", true)

  //only first few
  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,1);
  TEST_EQUAL((file.end() - file.begin()), 1)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)

  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,3);
  TEST_EQUAL((file.end() - file.begin()), 3)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "", true)

  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true,4);
  TEST_EQUAL((file.end() - file.begin()), 4)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)

  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true, -1, true);
  TEST_EQUAL((file.end() - file.begin()), 7)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "space_line", true)
  file_it += 4;
  TEST_EQUAL(String(*file_it).trim() == "last_line", true)

  file.load(OPENMS_GET_TEST_DATA_PATH("TextFile_test_infile.txt"),true, 4, true);
  TEST_EQUAL((file.end() - file.begin()), 4)
  file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "first_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "middle_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "space_line", true)
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "tab_line", true)
END_SECTION

START_SECTION((void store(const String& filename) ))
  TextFile file;

  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/does/not/exist/FileDoesNotExist.txt"))

  file.addLine("line1");
  file.addLine("line2\n");
  file.addLine("line3\r\n");
  String filename;
  NEW_TMP_FILE(filename);
  file.store(filename);
  file.load(filename);

  // validate loaded content
  TextFile::ConstIterator file_it = file.begin();
  TEST_EQUAL(String(*file_it).trim() == "line1",true);
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "line2",true);
  ++file_it;
  TEST_EQUAL(String(*file_it).trim() == "line3",true);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
