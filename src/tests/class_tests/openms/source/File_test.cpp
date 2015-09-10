// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <QDir>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static String getExecutablePath()))
	TEST_NOT_EQUAL(File::getExecutablePath().size(), 0)
END_SECTION

START_SECTION((static bool exists(const String &file)))
	TEST_EQUAL(File::exists("does_not_exists.txt"), false)
	TEST_EQUAL(File::exists(""), false)
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
	TEST_EQUAL(File::readable(""), false)
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
	TEST_EXCEPTION(Exception::FileNotFound, File::find("File.h"))
  String s_obo = File::find("CV/psi-ms.obo");
  TEST_EQUAL(s_obo.empty(), false);
  TEST_EQUAL(File::find(s_obo), s_obo); // iterative finding should return the identical file
  
  TEST_EXCEPTION(Exception::FileNotFound, File::find(""))
END_SECTION

START_SECTION((static String findDoc(const String& filename)))
	TEST_EXCEPTION(Exception::FileNotFound,File::findDoc("non-existing-documentation"))
  // should exist in every valid source tree (we cannot test for Doxyfile since doxygen might not be installed)
  TEST_EQUAL(File::findDoc("doxygen/Doxyfile.in").hasSuffix("Doxyfile.in"), true)
  // a file from the build tree
  TEST_EQUAL(File::findDoc("code_examples/cmake_install.cmake").hasSuffix("cmake_install.cmake"), true)
END_SECTION

START_SECTION((static String absolutePath(const String &file)))
	NOT_TESTABLE
END_SECTION

START_SECTION((static String path(const String &file)))
	NOT_TESTABLE
END_SECTION

START_SECTION((static String basename(const String &file)))
	TEST_EQUAL(File::basename("/souce/config/bla/bluff.h"), "bluff.h");
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

START_SECTION(static bool removeDirRecursively(const String &dir_name))
  QDir d;
  String dirname = File::getTempDirectory() + "/" + File::getUniqueName() + "/" + File::getUniqueName() + "/";
  TEST_EQUAL(d.mkpath(dirname.toQString()), TRUE);
  TextFile tf;
  tf.store(dirname + "test.txt");
  TEST_EQUAL(File::removeDirRecursively(dirname), true)
END_SECTION

START_SECTION(static String getTempDirectory())
  TEST_NOT_EQUAL(File::getTempDirectory(), String())
  TEST_EQUAL(File::exists(File::getTempDirectory()), true)
END_SECTION

START_SECTION(static String getUserDirectory())
  TEST_NOT_EQUAL(File::getUserDirectory(), String())
  TEST_EQUAL(File::exists(File::getUserDirectory()), true)
END_SECTION

START_SECTION(static Param getSystemParameters())
  Param p = File::getSystemParameters();
  TEST_EQUAL(p.size()>0, true)
  TEST_EQUAL(p.getValue("version"), VersionInfo::getVersion())
END_SECTION

START_SECTION(static String findDatabase(const String &db_name))

  TEST_EXCEPTION(Exception::FileNotFound, File::findDatabase("filedoesnotexists"))
  String db = File::findDatabase("./CV/unimod.obo");
  //TEST_EQUAL(db,"wtf")
  TEST_EQUAL(db.hasSubstring("share/OpenMS"), true)

END_SECTION


START_SECTION(static String findExecutable(const OpenMS::String& toolName))
{
	TEST_EXCEPTION(Exception::FileNotFound, File::findExecutable("executable_does_not_exist"))
	TEST_EQUAL(File::path(File::findExecutable("File_test")) + "/", File::getExecutablePath())
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
