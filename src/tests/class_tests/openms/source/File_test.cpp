// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

/////////////////////////////////////////////////////////////

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <atomic>
#include <ctime>
#include <fstream>
#include <filesystem>
#include <thread>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static std::string getExecutablePath()))
  TEST_NOT_EQUAL(File::getExecutablePath().size(), 0)
END_SECTION

START_SECTION((static bool exists(const std::string &file)))
  TEST_EQUAL(File::exists("does_not_exists.txt"), false)
  TEST_EQUAL(File::exists(""), false)
  TEST_EQUAL(File::exists(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), true)
END_SECTION

START_SECTION((static bool empty(const std::string &file)))
  TEST_EQUAL(File::empty("does_not_exists.txt"), true)
  TEST_EQUAL(File::empty(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)
  TEST_EQUAL(File::empty(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), false)
END_SECTION

START_SECTION((static UInt64 fileSize(const std::string& file)))
  TEST_EQUAL(File::fileSize("does_not_exists.txt"), -1)
  TEST_EQUAL(File::fileSize(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), 0)
  TEST_EQUAL(File::fileSize(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), 15)
END_SECTION

START_SECTION((static Int64 getModificationTime(const std::string& file)))
  TEST_EQUAL(File::getModificationTime("does_not_exists.txt"), -1)

  // Bracketed by the wall clock around a file this test writes itself. That pins both properties
  // the method promises: the number is that file's actual modification time, and it counts seconds
  // from the Unix epoch rather than from file_clock's implementation-defined one -- 2174 on
  // libstdc++, 1601 on MSVC -- so it means the same thing to a reader on another platform.
  //
  // Deliberately not read off a file in the source tree: distributions building reproducibly
  // (Debian, Nix) normalise every source timestamp to the epoch, which is a property of the build
  // environment and not of this method (see issue #10018). The few seconds of slack absorb coarse
  // filesystem timestamp granularity; the errors this guards against are off by billions.
  const Int64 t_before = static_cast<Int64>(std::time(nullptr));
  std::string tmp;
  NEW_TMP_FILE(tmp);
  { std::ofstream os(tmp.c_str()); os << "x"; }
  const Int64 t_after = static_cast<Int64>(std::time(nullptr));

  const Int64 t_tmp = File::getModificationTime(tmp);
  TEST_TRUE(t_tmp >= t_before - 5)
  TEST_TRUE(t_tmp <= t_after + 5)
END_SECTION

START_SECTION((static bool remove(const std::string &file)))
  //deleting non-existing file
  TEST_EQUAL(File::remove("does_not_exists.txt"), true)

  //deleting existing file
  std::string filename;
  NEW_TMP_FILE(filename);
  ofstream os;
  os.open (filename.c_str(), ofstream::out);
  os << "File_test dummy file to delete" << endl;
  os.close();
  TEST_EQUAL(File::remove(filename), true)
END_SECTION

START_SECTION((static bool readable(const std::string &file)))
  TEST_EQUAL(File::readable("does_not_exists.txt"), false)
  TEST_EQUAL(File::readable(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)
  TEST_EQUAL(File::readable(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")), true)
  TEST_EQUAL(File::readable(""), false)
END_SECTION

START_SECTION((static bool writable(const std::string &file)))
  TEST_EQUAL(File::writable("/this/file/cannot/be/written.txt"), false)
  TEST_EQUAL(File::writable(OPENMS_GET_TEST_DATA_PATH("File_test_empty.txt")), true)
  TEST_EQUAL(File::writable(OPENMS_GET_TEST_DATA_PATH("File_test_imaginary.txt")), true)
  TEST_EQUAL(File::writable(""), false)

  std::string filename;
  NEW_TMP_FILE(filename);
  TEST_EQUAL(File::writable(filename), true)

  // writable() is a query: asking the question must not change the answer, and must not
  // touch anything on disk that the caller did not ask about.
  {
    std::string probed;
    NEW_TMP_FILE(probed);
    { std::ofstream os(probed.c_str()); os << "payload"; }
    TEST_EQUAL(File::writable(probed), true)
    TEST_EQUAL(File::exists(probed), true)                        // still there
    TEST_EQUAL(size_t(std::filesystem::file_size(probed)), size_t(7)) // and still intact

    std::string absent;
    NEW_TMP_FILE(absent);
    std::filesystem::remove(absent);
    TEST_EQUAL(File::writable(absent), true)
    TEST_EQUAL(File::exists(absent), false)                       // probing left no litter
  }

  // A long path must not make writable() lie. writable() answers "could this be
  // created?" by creating a file of its own in the same directory, and that name is
  // longer than a typical caller's -- so close to the platform's path limit there is a
  // band where the caller's own file still fits but the probe no longer does. Reporting
  // "not writable" there is exactly the false negative this function exists to avoid.
  //
  // The limit is not the same on every platform (PATH_MAX, MAX_PATH, and whatever the
  // filesystem says), so rather than hard-code a length, grow a directory until the OS
  // refuses to go deeper and then sweep the last stretch of headroom. At every depth the
  // requirement is the same and is checked against reality rather than against an
  // expected answer: if the file can actually be created, writable() must say true.
  {
    File::TempDir long_path_dir;
    std::error_code ec;
    std::string deep = long_path_dir.getPath();
    const std::string chunk(60, 'd');

    // Cap the descent: a filesystem with no practical path limit must not loop forever.
    for (int level = 0; level < 200; ++level)
    {
      const std::string next = deep + chunk + "/";
      std::filesystem::create_directory(next, ec);
      if (ec || !std::filesystem::is_directory(next, ec)) break;
      deep = next;
    }

    int agreed = 0, checked = 0;
    for (size_t tail = 1; tail < 60; tail += 6)
    {
      const std::string leaf = deep + std::string(tail, 'e');
      std::filesystem::create_directory(leaf, ec);
      if (ec || !std::filesystem::is_directory(leaf, ec)) continue;

      // Ground truth: let the OS itself say whether this file can be created.
      const std::string target = leaf + "/o";
      std::ofstream os(target.c_str());
      const bool creatable = os.is_open() && os.good();
      os.close();
      if (creatable) std::filesystem::remove(target, ec);

      ++checked;
      if (File::writable(target) == creatable) ++agreed;
      std::filesystem::remove(leaf, ec);
    }

    TEST_NOT_EQUAL(checked, 0)      // the sweep has to have probed something
    TEST_EQUAL(agreed, checked)     // ... and writable() has to agree with the OS at every depth

    std::filesystem::remove_all(deep, ec);
  }

  // Concurrent callers must not sabotage each other. Both of the following are races, so a
  // single attempt proves nothing -- before the fix they went wrong in roughly 2% and 79% of
  // attempts respectively, which a one-shot check sails straight past. Repeat instead, and
  // use a spin barrier: without one the first thread finishes before the second even starts,
  // and the two never overlap at all.
  //
  // The spins yield rather than burn the core. That costs nothing natively (500 repeats run in
  // ~0.1 s either way) but matters by three orders of magnitude under valgrind, whose scheduler
  // hands a spinning thread a full quantum before switching: measured on the same binary,
  // 5 repeats take 67 s with a bare spin and 0.11 s with the yield.
  //
  // Note this section can only catch the bug where the threads genuinely overlap; pinned to a
  // single core they serialise and the pre-fix implementation passes too. It is a regression
  // guard on ordinary multi-core runners, not a proof on every machine.
  {
    const int repeats = 500;
    int both_writable = 0;   // two probes of the same new path must both answer 'true'
    int output_survived = 0; // a probe must not delete what a concurrent writer produced

    std::string shared;
    NEW_TMP_FILE(shared);
    std::string contended;
    NEW_TMP_FILE(contended);

    for (int i = 0; i < repeats; ++i)
    {
      // Two callers race to probe one path that does not exist yet. Probing used to create
      // and delete that very path, so one prober would delete the other's file and then
      // report it unwritable -- the intermittent 'Cannot write output file' seen on Windows
      // CI when two TOPP tools were given the same -out path.
      std::filesystem::remove(shared);
      std::atomic<int> go_probe{0};
      bool first = false, second = false;
      std::thread p1([&]() { while (go_probe == 0) { std::this_thread::yield(); } first = File::writable(shared); });
      std::thread p2([&]() { while (go_probe == 0) { std::this_thread::yield(); } second = File::writable(shared); });
      go_probe = 1;
      p1.join();
      p2.join();
      if (first && second) ++both_writable;

      // A probe running alongside a real writer must leave that writer's output alone.
      std::filesystem::remove(contended);
      std::atomic<int> go_write{0};
      std::thread writer([&]() { while (go_write == 0) { std::this_thread::yield(); } std::ofstream os(contended.c_str()); os << "important"; });
      std::thread prober([&]() { while (go_write == 0) { std::this_thread::yield(); } File::writable(contended); });
      go_write = 1;
      writer.join();
      prober.join();
      std::error_code ec;
      if (std::filesystem::exists(contended, ec) && std::filesystem::file_size(contended, ec) == 9) ++output_survived;
    }

    TEST_EQUAL(both_writable, repeats)
    TEST_EQUAL(output_survived, repeats)

    std::filesystem::remove(shared);
    std::filesystem::remove(contended);
  }
END_SECTION

START_SECTION((static std::string find(const std::string &filename, StringList directories=StringList())))
  TEST_EXCEPTION(Exception::FileNotFound, File::find("File.h"))
  std::string s_obo = File::find("CV/psi-ms.obo");
  TEST_EQUAL(s_obo.empty(), false);
  TEST_EQUAL(File::find(s_obo), s_obo); // iterative finding should return the identical file

  TEST_EXCEPTION(Exception::FileNotFound, File::find(""))

  // a missing standard data file should point the user at the resolved data path and OPENMS_DATA_PATH (issue #9636)
  {
    std::string msg;
    try
    {
      File::find("CHEMISTRY/this_file_does_not_exist_9636.xml");
    }
    catch (Exception::FileNotFound& e)
    {
      msg = e.getMessage();
    }
    TEST_EQUAL(msg.find("OPENMS_DATA_PATH") != std::string::npos, true)
    TEST_EQUAL(msg.find(File::getOpenMSDataPath()) != std::string::npos, true)
  }
END_SECTION

#ifdef ENABLE_DOCS
START_SECTION((static std::string findDoc(const std::string& filename)))
  TEST_EXCEPTION(Exception::FileNotFound,File::findDoc("non-existing-documentation"))
  // should exist in every valid source tree (we cannot test for Doxyfile since doxygen might not be installed)
  TEST_EQUAL(StringUtils::hasSuffix(File::findDoc("doxygen/Doxyfile.in"), "Doxyfile.in"), true)
  // a file from the build tree
  TEST_EQUAL(StringUtils::hasSuffix(File::findDoc("code_examples/cmake_install.cmake"), "cmake_install.cmake"), true)
END_SECTION
#endif

START_SECTION((static std::string absolutePath(const std::string &file)))
  NOT_TESTABLE
END_SECTION

START_SECTION((static std::string path(const std::string &file)))
  TEST_EQUAL(File::path("/source/config/bla/bluff.h"), "/source/config/bla");
  TEST_EQUAL(File::path("c:\\config\\bla\\tuff.h"), "c:\\config\\bla");
  TEST_EQUAL(File::path("filename_only.h"), "."); // useful when you want to reassemble a full path using path() + '/' + basename(), but the input is only a filename
  TEST_EQUAL(File::path("/path/only/"), "/path/only");
END_SECTION

START_SECTION((static std::string basename(const std::string &file)))
  TEST_EQUAL(File::basename("/source/config/bla/bluff.h"), "bluff.h");
  TEST_EQUAL(File::basename("filename_only.h"), "filename_only.h");
  TEST_EQUAL(File::basename("/path/only/"), "");
END_SECTION

START_SECTION((static std::string stemName(const std::string &file)))
  // basic: strips known extension from full path
  TEST_EQUAL(File::stemName("/path/to/sample.mzML"), "sample");
  // compound extension: .mzML.gz is a known compound extension
  TEST_EQUAL(File::stemName("/path/to/sample.mzML.gz"), "sample");
  // unknown extension: strips last dot segment
  TEST_EQUAL(File::stemName("/path/to/file.txt"), "file");
  // unknown compound: only strips known part
  TEST_EQUAL(File::stemName("/path/to/file.txt.tgz"), "file.txt");
  // no extension
  TEST_EQUAL(File::stemName("/path/to/file"), "file");
  // filename only (no path)
  TEST_EQUAL(File::stemName("experiment.featureXML"), "experiment");
  // empty string
  TEST_EQUAL(File::stemName(""), "");
  // dotted directory, no extension on file
  TEST_EQUAL(File::stemName("/home.with.dot/filename"), "filename");
  // Windows path
  TEST_EQUAL(File::stemName("c:\\data\\sample.idXML"), "sample");
  // extension-only name
  TEST_EQUAL(File::stemName(".mzML"), "");
END_SECTION

START_SECTION((static std::string extension(const std::string &file)))
  // known extension
  TEST_EQUAL(File::extension("/path/to/sample.mzML"), ".mzML");
  // compound extension
  TEST_EQUAL(File::extension("/path/to/sample.mzML.gz"), ".mzML.gz");
  // unknown extension
  TEST_EQUAL(File::extension("/path/to/file.txt"), ".txt");
  // no extension
  TEST_EQUAL(File::extension("/path/to/file"), "");
  // filename only
  TEST_EQUAL(File::extension("experiment.featureXML"), ".featureXML");
  // empty string
  TEST_EQUAL(File::extension(""), "");
  // dotted directory, no extension on file
  TEST_EQUAL(File::extension("/home.with.dot/filename"), "");
  // Windows path
  TEST_EQUAL(File::extension("c:\\data\\sample.idXML"), ".idXML");
  // extension-only name
  TEST_EQUAL(File::extension(".mzML"), ".mzML");
END_SECTION

START_SECTION((static StringList listDirectories(const std::string &dir)))
  // create temp structure with subdirectories
  File::TempDir tdir;
  std::string base = tdir.getPath();
  File::makeDir(base + "/subA");
  File::makeDir(base + "/subB");
  // also create a file (should NOT appear in results)
  {
    std::ofstream f(std::string(base + "/afile.txt"));
    f << "test";
  }

  StringList dirs = File::listDirectories(base);
  TEST_EQUAL(dirs.size(), 2);
  // results are sorted
  TEST_TRUE(StringUtils::hasSuffix(dirs[0], "subA"));
  TEST_TRUE(StringUtils::hasSuffix(dirs[1], "subB"));

  // non-existent directory returns empty list
  StringList empty = File::listDirectories("/nonexistent_path_xyz");
  TEST_EQUAL(empty.size(), 0);

  // empty string returns empty list (not a directory)
  StringList from_empty = File::listDirectories("");
  // just verify it doesn't crash - result depends on cwd
END_SECTION

START_SECTION((static bool fileList(const std::string &dir, const std::string &file_pattern, StringList &output, bool full_path=false)))
  StringList vec;
  TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "*.bliblaluff", vec), false);
  TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "File_test_text.txt", vec), true);
  TEST_EQUAL(vec[0], "File_test_text.txt");
  TEST_EQUAL(File::fileList(OPENMS_GET_TEST_DATA_PATH(""), "File_test_text.txt", vec, true), true);
  // can't use "path + separator + filename", because sep. might be "/" or "\\"
  TEST_EQUAL(StringUtils::hasPrefix(vec[0], OPENMS_GET_TEST_DATA_PATH("")), true);
  TEST_EQUAL(StringUtils::hasSuffix(vec[0], "File_test_text.txt"), true);
END_SECTION

START_SECTION((static std::string getUniqueName(bool include_hostname = true)))
  std::string unique_name = File::getUniqueName();
  std::string unique_name_no_host = File::getUniqueName(false);
  TEST_EQUAL(unique_name.size() > unique_name_no_host.size(), true)
  
  // test if the string consists of four parts
  StringList split;
  StringUtils::split(unique_name, '_', split);
  TEST_EQUAL(split.size() >= 4, true) // if name of machine also contains '_' ...
END_SECTION

START_SECTION((static std::string getOpenMSDataPath()))
  NOT_TESTABLE
END_SECTION

START_SECTION((static bool isDirectory(const std::string& path)))
  TEST_EQUAL(File::isDirectory(""),false)
  TEST_EQUAL(File::isDirectory("."),true)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("")),true)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("does_not_exists.txt")),false)
  TEST_EQUAL(File::isDirectory(OPENMS_GET_TEST_DATA_PATH("File_test_text.txt")),false)
END_SECTION

// make source directory and copy it to new location
// check copy function and if file exists in target path
START_SECTION(static bool copyDirRecursively(const std::string &fromDir, const std::string &toDir, File::CopyOptions option = CopyOptions::OVERWRITE))
  // folder OpenMS/src/tests/class_tests/openms/data/XMassFile_test 
  std::string source_name = OPENMS_GET_TEST_DATA_PATH("XMassFile_test");
  std::string target_name = File::getTempDirectory() + "/" + File::getUniqueName() + "/"; 
  // test canonical path
  // Rejecting the self-copy is reported on the error log. That report is expected here -- the
  // assertion below is what checks it -- so keep it off the console: an 'Error:' line in the
  // output of a passing test reads as a broken build to anyone packaging OpenMS (issue #10019).
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogError(), std::cerr);
    TEST_EQUAL(File::copyDirRecursively(source_name,source_name),false)
  }
  // test default
  TEST_EQUAL(File::copyDirRecursively(source_name,target_name),true)
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
  // SKIP announces every file it leaves alone on the warning log; expected, and equally noisy.
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogWarn(), std::cerr);
    TEST_EQUAL(File::copyDirRecursively(source_name,target_name, File::CopyOptions::SKIP),true)
  }
  infile.open(target_name + "/pdata/1/proc"); 
  infile.seekg(0,infile.end);
  file_size = infile.tellg();
  infile.close();
  TEST_EQUAL(file_size,50)
  // test option overwrite
  TEST_EQUAL(File::copyDirRecursively(source_name,target_name, File::CopyOptions::OVERWRITE),true)
  infile.open(target_name + "/pdata/1/proc"); 
  infile.seekg(0,infile.end);
  file_size = infile.tellg();
  infile.close();
  TEST_EQUAL(file_size,3558)
  // test option cancel 
  TEST_EQUAL(File::copyDirRecursively(source_name,target_name, File::CopyOptions::CANCEL),false)
  // remove temporary directory after testing
  File::removeDirRecursively(target_name);
END_SECTION

START_SECTION(static bool removeDirRecursively(const std::string &dir_name))
  std::string dirname = File::getTempDirectory() + "/" + File::getUniqueName() + "/" + File::getUniqueName() + "/";
  TEST_TRUE(File::makeDir(dirname));
  TextFile tf;
  tf.store(dirname + "test.txt");
  TEST_EQUAL(File::removeDirRecursively(dirname), true)
  END_SECTION

START_SECTION(static bool makeDir(const std::string& dir_name))
  File::TempDir tdir;
  std::string dirname = tdir.getPath() + "/" + File::getUniqueName() + "/" + File::getUniqueName() + "/";
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

START_SECTION(static std::string getTempDirectory())
  TEST_NOT_EQUAL(File::getTempDirectory(), std::string())
  TEST_EQUAL(File::exists(File::getTempDirectory()), true)
END_SECTION

START_SECTION(static std::string getUserDirectory())
  TEST_NOT_EQUAL(File::getUserDirectory(), std::string())
  TEST_EQUAL(File::exists(File::getUserDirectory()), true)

  // set user directory to a path set by environmental variable and test that
  // it is correctly set (no changes on the file system occur)
  std::string dirname = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
  TEST_EQUAL(File::makeDir(dirname), true);
#ifdef OPENMS_WINDOWSPLATFORM
  _putenv_s("OPENMS_HOME_PATH", dirname.c_str());  
#else
  setenv("OPENMS_HOME_PATH", dirname.c_str(), 0);  
#endif
  TEST_EQUAL(File::getUserDirectory(), dirname)
  // Note: this does not guarantee any more that the user directory or an
  // OpenMS.ini file exists at the new location.
END_SECTION

START_SECTION((static std::string getOpenMSConfigDir()))
  std::string config_dir = File::getOpenMSConfigDir();
  TEST_NOT_EQUAL(config_dir, std::string())
  // every platform branch resolves to a folder named "OpenMS" with no trailing separator
  TEST_EQUAL(StringUtils::hasSuffix(config_dir, "OpenMS"), true)
  TEST_EQUAL(StringUtils::hasSuffix(config_dir, "/"), false)
#ifdef __unix__
  // on unix-like systems, XDG_CONFIG_HOME takes precedence when set
  const char* xdg_backup = getenv("XDG_CONFIG_HOME");
  setenv("XDG_CONFIG_HOME", "/tmp/openms_xdg_test", 1);
  TEST_EQUAL(File::getOpenMSConfigDir(), "/tmp/openms_xdg_test/OpenMS")
  // restore previous environment to avoid side effects on later tests
  if (xdg_backup) { setenv("XDG_CONFIG_HOME", xdg_backup, 1); }
  else { unsetenv("XDG_CONFIG_HOME"); }
#endif
END_SECTION

START_SECTION(static Param getSystemParameters())
  Param p = File::getSystemParameters();
  TEST_EQUAL(!p.empty(), true)
  TEST_EQUAL(p.getValue("version"), VersionInfo::getVersion())
END_SECTION

START_SECTION(static std::string findDatabase(const std::string &db_name))

  // findDatabase() logs the miss before rethrowing; the exception is what this asserts on.
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogError(), std::cerr);
    TEST_EXCEPTION(Exception::FileNotFound, File::findDatabase("filedoesnotexists"))
  }
  // The success path is chatty too -- it announces the resolved path on the info log (which goes
  // to stdout, so the stderr check above does not cover it), and since nothing flushes it until
  // teardown it surfaces *after* the test's PASSED banner. Same packaging-log noise as the errors.
  std::string db;
  {
    Logger::LogSinkGuard quiet(getThreadLocalLogInfo(), std::cout);
    db = File::findDatabase("./CV/unimod.obo");
  }
  //TEST_EQUAL(db,"wtf")
  TEST_EQUAL(StringUtils::hasSubstring(db, "share/OpenMS"), true)

END_SECTION

START_SECTION(static bool findExecutable(std::string& exe_filename))
{
  //NOT_TESTABLE // since it depends on PATH

  // this test is somewhat brittle, but should work on most platforms (revert to NOT_TESTABLE if this does not work)
#ifdef OPENMS_WINDOWSPLATFORM
  std::string find = "cmd";
  TEST_EQUAL(File::findExecutable(find), true)
  TEST_EQUAL(StringUtils::toUppered(StringUtils::suffix(find, 7)), "CMD.EXE") // should be C:\Windows\System32\cmd.exe or similar
#else
  std::string find = "echo";
  TEST_EQUAL(File::findExecutable(find), true)
  TEST_EQUAL(StringUtils::hasSuffix(find, "echo"), true) // should be /usr/bin/echo or similar
#endif

  // negative case: a name that is neither an existing path nor on PATH is not found
  std::string missing = "openms_no_such_executable_4711";
  TEST_EQUAL(File::findExecutable(missing), false)
  TEST_EQUAL(missing, "openms_no_such_executable_4711") // left unchanged when not found
}
END_SECTION

START_SECTION(static StringList getPathLocations(const std::string& path))
{
  // set env-variables is not portable across platforms, thus we inject the PATH values
#ifdef OPENMS_WINDOWSPLATFORM
  std::string test_paths=R"(C:\WINDOWS\CCM;C:\WINDOWS\system32\config\systemprofile\AppData\Local\Microsoft\WindowsApps;C:\Program Files (x86)\Git\cmd)";
#else
  std::string test_paths="/usr/local/bin:/usr/bin:/bin:/usr/local/games:/usr/games";
#endif
  StringList l = File::getPathLocations(test_paths);
#ifdef OPENMS_WINDOWSPLATFORM
  TEST_EQUAL(ListUtils::contains(l, "C:/Program Files (x86)/Git/cmd/"), true)
#else
  TEST_EQUAL(ListUtils::contains(l, "/usr/bin/"), true)
#endif

  StringList empty = File::getPathLocations("");
  TEST_EQUAL(empty.empty(), true)

  // Also exercise the no-argument overload, which reads PATH safely even when the
  // environment lookup would otherwise return nullptr.
  File::getPathLocations();
}
END_SECTION


START_SECTION(static std::string findSiblingTOPPExecutable(const std::string& toolName))
{
  TEST_EXCEPTION(Exception::FileNotFound, File::findSiblingTOPPExecutable("executable_does_not_exist"))
  TEST_EQUAL(File::path(File::findSiblingTOPPExecutable("File_test")) + "/", File::getExecutablePath())

  // Sibling resolution looks next to the running executable only; it does NOT fall back to a
  // PATH lookup. A binary that exists only on PATH (e.g. the "echo"/"cmd" resolved by findExecutable
  // above) is therefore reported as not found here. This pins the current "siblings only" contract.
#ifdef OPENMS_WINDOWSPLATFORM
  TEST_EXCEPTION(Exception::FileNotFound, File::findSiblingTOPPExecutable("cmd"))
#else
  TEST_EXCEPTION(Exception::FileNotFound, File::findSiblingTOPPExecutable("echo"))
#endif
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
  std::string path;
  {
    File::TempDir dir;
    path = dir.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 0)
  if (File::exists(path)) File::removeDir(path);
  {
    File::TempDir dir2(true);
    path = dir2.getPath();
    TEST_EQUAL(File::exists(path), 1)
  }
  TEST_EQUAL(File::exists(path), 1)
  if (File::exists(path)) File::removeDir(path);
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
