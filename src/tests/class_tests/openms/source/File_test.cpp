// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <QDir>

#include <fstream>
#include <filesystem>

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

START_SECTION((static UInt64 fileSize(const String& file)))
  TEST_EQUAL(File::fileSize("does_not_exists.txt"), -1)
  TEST_EQUAL(File::fileSize(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), 0)
  TEST_EQUAL(File::fileSize(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), 15)
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

#ifdef ENABLE_DOCS
START_SECTION((static String findDoc(const String& filename)))
  TEST_EXCEPTION(Exception::FileNotFound,File::findDoc("non-existing-documentation"))
  // should exist in every valid source tree (we cannot test for Doxyfile since doxygen might not be installed)
  TEST_EQUAL(File::findDoc("doxygen/Doxyfile.in").hasSuffix("Doxyfile.in"), true)
  // a file from the build tree
  TEST_EQUAL(File::findDoc("code_examples/cmake_install.cmake").hasSuffix("cmake_install.cmake"), true)
END_SECTION
#endif

START_SECTION((static String absolutePath(const String &file)))
  NOT_TESTABLE
END_SECTION

START_SECTION((static String path(const String &file)))
  TEST_EQUAL(File::path("/source/config/bla/bluff.h"), "/source/config/bla");
  TEST_EQUAL(File::path("c:\\config\\bla\\tuff.h"), "c:\\config\\bla");
  TEST_EQUAL(File::path("filename_only.h"), "."); // useful when you want to reassemble a full path using path() + '/' + basename(), but the input is only a filename
  TEST_EQUAL(File::path("/path/only/"), "/path/only");
END_SECTION

START_SECTION((static String basename(const String &file)))
  TEST_EQUAL(File::basename("/source/config/bla/bluff.h"), "bluff.h");
  TEST_EQUAL(File::basename("filename_only.h"), "filename_only.h");
  TEST_EQUAL(File::basename("/path/only/"), "");
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

START_SECTION((static String getUniqueName(bool include_hostname = true)))
  String unique_name = File::getUniqueName();
  String unique_name_no_host = File::getUniqueName(false);
  TEST_EQUAL(unique_name.size() > unique_name_no_host.size(), true)
  
  // test if the string consists of four parts
  StringList split;
  unique_name.split('_', split);
  TEST_EQUAL(split.size() >= 4, true) // if name of machine also contains '_' ...
END_SECTION

START_SECTION((static String getOpenMSDataPath()))
  NOT_TESTABLE
END_SECTION

START_SECTION((static bool isDirectory(const String& path)))
  TEST_EQUAL(File::isDirectory(""),false)
  TEST_EQUAL(File::isDirectory("."),true)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("")),true)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("does_not_exists.txt")),false)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")),false)
END_SECTION

// make source directory and copy it to new location
// check copy function and if file exists in target path
START_SECTION(static bool copyDirRecursively(const QString &fromDir, const QString &toDir,File::CopyOptions option = CopyOptions::OVERWRITE))
  // folder OpenMS/src/tests/class_tests/openms/data/XMassFile_test 
  String source_name = OPENMS_GET_TEST_DATA_PATH("XMassFile_test");
  String target_name = File::getTempDirectory() + "/" + File::getUniqueName() + "/"; 
  // test canonical path
  TEST_EQUAL(File::copyDirRecursively(source_name.toQString(),source_name.toQString()),false)
  // test default
  TEST_EQUAL(File::copyDirRecursively(source_name.toQString(),target_name.toQString()),true)
  TEST_EQUAL(File::exists(target_name + "/pdata/1/proc"),true);
  // overwrite file content 
  std::ofstream ow_ofs;
  ow_ofs.open(target_name + "/pdata/1/proc");
  ow_ofs << "This line can be used to test the overwrite option";
  ow_ofs.close();
  // check file size 
  std::ifstream infile;
  infile.open(target_name + "/pdata/1/proc"); 
  infile.seekg(0,infile.end);
  long file_size = infile.tellg();
  infile.close();
  TEST_EQUAL(file_size,50)
  // test option skip
  TEST_EQUAL(File::copyDirRecursively(source_name.toQString(),target_name.toQString(), File::CopyOptions::SKIP),true)
  infile.open(target_name + "/pdata/1/proc"); 
  infile.seekg(0,infile.end);
  file_size = infile.tellg();
  infile.close();
  TEST_EQUAL(file_size,50)
  // test option overwrite
  TEST_EQUAL(File::copyDirRecursively(source_name.toQString(),target_name.toQString(), File::CopyOptions::OVERWRITE),true)
  infile.open(target_name + "/pdata/1/proc"); 
  infile.seekg(0,infile.end);
  file_size = infile.tellg();
  infile.close();
  TEST_EQUAL(file_size,3558)
  // test option cancel 
  TEST_EQUAL(File::copyDirRecursively(source_name.toQString(),target_name.toQString(), File::CopyOptions::CANCEL),false)
  // remove temporary directory after testing
  File::removeDirRecursively(target_name);
END_SECTION

START_SECTION(static bool removeDirRecursively(const String &dir_name))
  QDir d;
  String dirname = File::getTempDirectory() + "/" + File::getUniqueName() + "/" + File::getUniqueName() + "/";
  TEST_TRUE(d.mkpath(dirname.toQString()));
  TextFile tf;
  tf.store(dirname + "test.txt");
  TEST_EQUAL(File::removeDirRecursively(dirname), true)
  END_SECTION

START_SECTION(static bool makeDir(const String& dir_name))
  File::TempDir tdir;
  String dirname = tdir.getPath() + "/" + File::getUniqueName() + "/" + File::getUniqueName() + "/";
  // absolute path
  TEST_FALSE(File::isDirectory(dirname))
  TEST_TRUE(File::makeDir(dirname))
  TEST_TRUE(File::isDirectory(dirname))
  // a relative path
  auto current_path = std::filesystem::current_path(); // get current path
  filesystem::current_path(std::filesystem::path(dirname.c_str())); // set current path to dirname
  TEST_TRUE(File::makeDir("subdir/333"))
  TEST_TRUE(File::isDirectory("./subdir/333/"))
  // try create something which should be forbidden
#if defined(OPENMS_WINDOWSPLATFORM)
  TEST_FALSE(File::makeDir("c:\\te:st")) // ':' is not allowed in path on Windows; Unix pretty much allows everything
#endif
  std::filesystem::current_path(current_path); // reset current path (enable deletion of dirname)
END_SECTION

START_SECTION(static String getTempDirectory())
  TEST_NOT_EQUAL(File::getTempDirectory(), String())
  TEST_EQUAL(File::exists(File::getTempDirectory()), true)
END_SECTION

START_SECTION(static String getUserDirectory())
  TEST_NOT_EQUAL(File::getUserDirectory(), String())
  TEST_EQUAL(File::exists(File::getUserDirectory()), true)

  // set user directory to a path set by environmental variable and test that
  // it is correctly set (no changes on the file system occur)
  QDir d;
  String dirname = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
  TEST_EQUAL(d.mkpath(dirname.toQString()), true);
#ifdef OPENMS_WINDOWSPLATFORM
  _putenv_s("OPENMS_HOME_PATH", dirname.c_str());  
#else
  setenv("OPENMS_HOME_PATH", dirname.c_str(), 0);  
#endif
  TEST_EQUAL(File::getUserDirectory(), dirname)
  // Note: this does not guarantee any more that the user directory or an
  // OpenMS.ini file exists at the new location.
END_SECTION

START_SECTION(static Param getSystemParameters())
  Param p = File::getSystemParameters();
  TEST_EQUAL(!p.empty(), true)
  TEST_EQUAL(p.getValue("version"), VersionInfo::getVersion())
END_SECTION

START_SECTION(static String findDatabase(const String &db_name))

  TEST_EXCEPTION(Exception::FileNotFound, File::findDatabase("filedoesnotexists"))
  String db = File::findDatabase("./CV/unimod.obo");
  //TEST_EQUAL(db,"wtf")
  TEST_EQUAL(db.hasSubstring("share/OpenMS"), true)

END_SECTION

START_SECTION(static bool findExecutable(OpenMS::String& exe_filename))
{
  //NOT_TESTABLE // since it depends on PATH

  // this test is somewhat brittle, but should work on most platforms (revert to NOT_TESTABLE if this does not work)
#ifdef OPENMS_WINDOWSPLATFORM
  String find = "cmd";
  TEST_EQUAL(File::findExecutable(find), true)
  TEST_EQUAL(find.suffix(7).toUpper(), "CMD.EXE") // should be C:\Windows\System32\cmd.exe or similar
#else
  String find = "echo";
  TEST_EQUAL(File::findExecutable(find), true)
  TEST_EQUAL(find.hasSuffix("echo"), true) // should be /usr/bin/echo or similar
#endif
}
END_SECTION

START_SECTION(static StringList getPathLocations(const String& path = std::getenv("PATH")))
{
  // set env-variables is not portable across platforms, thus we inject the PATH values
#ifdef OPENMS_WINDOWSPLATFORM
  String test_paths=R"(C:\WINDOWS\CCM;C:\WINDOWS\system32\config\systemprofile\AppData\Local\Microsoft\WindowsApps;C:\Program Files (x86)\Git\cmd)";
#else
  String test_paths="/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games";
#endif
  StringList l = File::getPathLocations(test_paths);
#ifdef OPENMS_WINDOWSPLATFORM
  TEST_EQUAL(ListUtils::contains(l, "C:/Program Files (x86)/Git/cmd/"), true)
#else
  TEST_EQUAL(ListUtils::contains(l, "/usr/bin/"), true)
#endif
}
END_SECTION


START_SECTION(static String findSiblingTOPPExecutable(const OpenMS::String& toolName))
{
  TEST_EXCEPTION(Exception::FileNotFound, File::findSiblingTOPPExecutable("executable_does_not_exist"))
  TEST_EQUAL(File::path(File::findSiblingTOPPExecutable("File_test")) + "/", File::getExecutablePath())
}
END_SECTION

START_SECTION(File::TempDir(bool keep_dir = false))
{
  File::TempDir* dir = new File::TempDir();
  File::TempDir* nullPointer = nullptr;
  TEST_NOT_EQUAL(dir, nullPointer)
  TEST_EQUAL(File::exists((*dir).getPath()),1)
  delete dir;
}
END_SECTION

START_SECTION(File::~TempDir())
{
  String path;
  {
    File::TempDir dir;
    path = dir.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 0)
  if (File::exists(path)) File::removeDir(path.toQString());
  {
    File::TempDir dir2(true);
    path = dir2.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 1)
  if (File::exists(path)) File::removeDir(path.toQString());
}
END_SECTION


START_SECTION(File::download(std::string url, std::string filename))
{
  std::string url = R"(http://raw.githubusercontent.com/OpenMS/images/master/logos/OpenMS_transparent_background.png)";
  std::string folder = File::getTempDirectory();
  File::download(url, folder);
  std::string output_file_path = folder + "/OpenMS_transparent_background.png";

  TEST_EQUAL(File::exists(output_file_path), 1);
  if (File::exists(output_file_path))
  {
    File::removeDir(QString(output_file_path.c_str()));
  }
}
END_SECTION

START_SECTION(static File::MatchingFileListsStatus validateMatchingFileNames(const StringList& sl1, const StringList& sl2, bool basename, bool ignore_extension))
{
  // Test exact match
  {
    StringList list1 = {"file1.txt", "file2.txt"};
    StringList list2 = {"file1.txt", "file2.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2) ==  File::MatchingFileListsStatus::MATCH)
  }

  // Test order mismatch
  {
    StringList list1 = {"file1.txt", "file2.txt"};
    StringList list2 = {"file2.txt", "file1.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2) == File::MatchingFileListsStatus::ORDER_MISMATCH)
  }

  // Test different sets
  {
    StringList list1 = {"file1.txt", "file2.txt"};
    StringList list2 = {"file1.txt", "file3.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2) == File::MatchingFileListsStatus::SET_MISMATCH)
  }

  // Test different counts
  {
    StringList list1 = {"file1.txt", "file2.txt"};
    StringList list2 = {"file1.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2) ==  File::MatchingFileListsStatus::SET_MISMATCH)
  }

  // Test basename comparison
  {
    StringList list1 = {"/path/to/file1.txt", "/different/path/file2.txt"};
    StringList list2 = {"/other/path/file1.txt", "/somewhere/file2.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2, true, false) ==  File::MatchingFileListsStatus::MATCH)
  }

  // Test basename with order mismatch
  {
    StringList list1 = {"/path/to/file1.txt", "/different/path/file2.txt"};
    StringList list2 = {"/somewhere/file2.txt", "/other/path/file1.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2, true, false) ==  File::MatchingFileListsStatus::ORDER_MISMATCH)
  }

  // Test ignore extension
  {
    StringList list1 = {"file1.txt", "file2.mzML"};
    StringList list2 = {"file1.mzML", "file2.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2, false, true) == File::MatchingFileListsStatus::MATCH)
  }

  // Test ignore extension with different basenames
  {
    StringList list1 = {"file1.txt", "file2.mzML"};
    StringList list2 = {"file1.mzML", "file3.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2, false, true) ==  File::MatchingFileListsStatus::SET_MISMATCH)
  }

  // Test with both basename and ignore extension
  {
    StringList list1 = {"/path/to/file1.txt", "/different/path/file2.mzML"};
    StringList list2 = {"/other/path/file1.mzML", "/somewhere/file2.txt"};
    TEST_TRUE(File::validateMatchingFileNames(list1, list2, true, true) == File::MatchingFileListsStatus::MATCH)
  }

  // Test with empty lists
  {
    StringList list1, list2;
    TEST_TRUE(File::validateMatchingFileNames(list1, list2) == File::MatchingFileListsStatus::MATCH)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
