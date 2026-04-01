// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/PathUtils.h>
#include <OpenMS/openms_data_path.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <algorithm>
#include <atomic>
#include <cerrno>
#include <filesystem>
#include <fstream>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h> // for GetCurrentProcessId() && GetModuleFileName() && GetComputerNameA()
#include <Shlwapi.h> // for PathMatchSpecA
#include <io.h>      // for _access_s
#pragma comment(lib, "Shlwapi.lib")
#else
#include <fnmatch.h>
#include <unistd.h> // for gethostname()
#endif

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h> // for readLink() and getpid()
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

namespace fs = std::filesystem;

using namespace std;

namespace OpenMS
{

  File::TempDir::TempDir(bool keep_dir)
    : keep_dir_(keep_dir)
  {
    temp_dir_ = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'\n";
    fs::create_directories(to_path(temp_dir_));
  };

  File::TempDir::TempDir(const String& base_dir, bool keep_dir)
    : keep_dir_(keep_dir)
  {
    // Create a unique subdirectory under the provided base_dir
    temp_dir_ = base_dir;
    if (!temp_dir_.empty() && !temp_dir_.hasSuffix("/"))
    {
      temp_dir_ += "/";
    }
    temp_dir_ += "OpenMSTempDir_" + File::getUniqueName() + "/";
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'\n";
    fs::create_directories(to_path(temp_dir_));
  };

  File::TempDir::~TempDir()
  {
    if (keep_dir_)
    {
      OPENMS_LOG_DEBUG << "Keeping temporary files in directory '" << temp_dir_ << '\n';
      return;
    }

    File::removeDirRecursively(temp_dir_);
  };

  const String& File::TempDir::getPath() const
  {
    return temp_dir_;
  }

  String File::getExecutablePath()
  {
    // see http://stackoverflow.com/questions/1023306/finding-current-executables-path-without-proc-self-exe/1024937#1024937 for more OS' (if needed)
    // Use immediately evaluated lambda to protect static variable from concurrent access.
    static const String spath = [&]() -> String {
        String rpath = "";

        char path[1024]; // maximum path length

#ifdef OPENMS_WINDOWSPLATFORM
        int size = sizeof(path);
        if (GetModuleFileNameA(NULL, path, size))
#elif  defined(__APPLE__)
        uint size = sizeof(path);
        if (_NSGetExecutablePath(path, &size) == 0)
#else // LINUX
        // note: implementation as suggested by readlink man page
        ssize_t len = ::readlink("/proc/self/exe", path, sizeof(path)-1);
        if (len != -1) //add 0 terminator at end
        {
          path[len] = '\0';
        }

        if (len != -1)
#endif
        {
          rpath = File::path(String(path));
          if (File::exists(rpath)) // check if directory exists
          {
            // ensure path ends with a "/", such that we can just write path + "ToolX", and to not worry about if its empty or a path.
            rpath.ensureLastChar('/');
          }
          else
          {
            std::cerr << "Path '" << rpath << "' extracted from Executable Path '" << path << "' does not exist! Returning empty string!\n";
            rpath = "";
          }
        } else {
          std::cerr << "Cannot get Executable Path! Not using a path prefix!\n";
        }

        return rpath;
    }();
    return spath;
  }

  bool File::exists(const String& file)
  {
    return fs::exists(to_path(file));
  }

  bool File::empty(const String& file)
  {
    std::error_code ec;
    auto p = to_path(file);
    if (!fs::exists(p, ec)) return true;
    auto sz = fs::file_size(p, ec);
    if (ec) return true; // e.g. if it is a directory
    return sz == 0;
  }

  bool File::executable(const String& file)
  {
    auto p = to_path(file);
    std::error_code ec;
    auto st = fs::status(p, ec);
    if (ec || !fs::exists(st)) return false;
#ifdef OPENMS_WINDOWSPLATFORM
    // On Windows, all files are considered executable
    return true;
#else
    auto perms = st.permissions();
    return (perms & (fs::perms::owner_exec | fs::perms::group_exec | fs::perms::others_exec)) != fs::perms::none;
#endif
  }

  UInt64 File::fileSize(const String& file)
  {
    if (!File::exists(file)) return -1;

    std::error_code ec;
    auto sz = fs::file_size(to_path(file), ec);
    if (ec) return -1;
    return sz;
  }

  bool File::rename(const String& from, const String& to, bool overwrite_existing, bool verbose)
  {
    // check for equality
    std::error_code ec;
    auto canon_from = fs::canonical(to_path(from), ec);
    if (!ec)
    {
      auto canon_to = fs::canonical(to_path(to), ec);
      if (!ec && canon_from == canon_to)
      { // same file; no need to do anything
        return true;
      }
    }

    // existing file? std::filesystem::rename overwrites on POSIX but not always on Windows.
    // To be consistent, handle overwrite explicitly:
    if (overwrite_existing && exists(to) && !remove(to))
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not overwrite existing file '" << to << "'\n";
      }
      return false;
    }
    // move the file to the actual destination:
    std::error_code rename_ec;
    fs::rename(to_path(from), to_path(to), rename_ec);
    if (rename_ec)
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not move '" << from << "' to '" << to << "'\n";
      }
      return false;
    }
    return true;
  }

  bool File::copyDirRecursively(const String& from_dir, const String& to_dir, File::CopyOptions option)
  {
    auto source_path = to_path(from_dir);
    auto target_path = to_path(to_dir);

    std::error_code ec;
    auto canonical_source = fs::canonical(source_path, ec);
    // if source doesn't exist, canonical will fail
    if (ec)
    {
      OPENMS_LOG_ERROR << "Error: Source directory '" << from_dir << "' does not exist.\n";
      return false;
    }

    // If target exists, check canonical to avoid copy onto self
    if (fs::exists(target_path, ec))
    {
      auto canonical_target = fs::canonical(target_path, ec);
      if (!ec && canonical_source == canonical_target)
      {
        OPENMS_LOG_ERROR << "Error: Could not copy  " << from_dir << " to " << to_dir << ". Same path given.\n";
        return false;
      }
    }

    // make directory if not present
    fs::create_directories(target_path, ec);

    // copy folder recursively
    for (const auto& entry : fs::directory_iterator(source_path))
    {
      auto target_file = target_path / entry.path().filename();

      if (entry.is_directory())
      {
        if (!copyDirRecursively(entry.path().generic_string(), target_file.generic_string(), option))
        {
          return false;
        }
      }
      else
      {
        if (fs::exists(target_file))
        {
          switch (option)
          {
            case CopyOptions::CANCEL:
              return false;
            case CopyOptions::SKIP:
              OPENMS_LOG_WARN << "The file " << entry.path().filename().string() << " was skipped.\n";
              continue;
            case CopyOptions::OVERWRITE:
              fs::remove(target_file, ec);
              break;
          }
        }
        std::error_code copy_ec;
        fs::copy_file(entry.path(), target_file, copy_ec);
        if (copy_ec)
        {
          return false;
        }
      }
    }
    return true;
  }

  bool File::copy(const String& from, const String& to)
  {
    std::error_code ec;
    return fs::copy_file(to_path(from), to_path(to), ec);
  }

  bool File::remove(const String& file)
  {
    if (!exists(file))
    {
      return true;
    }
    if (std::remove(file.c_str()) != 0)
    {
      return false;
    }
    return true;
  }

  bool File::removeDir(const String& dir_name)
  {
    std::error_code ec;
    fs::remove_all(to_path(dir_name), ec); // recursive, matching original Qt behavior
    if (ec)
    {
      std::cerr << "Could not remove directory " << dir_name << ": " << ec.message() << std::endl;
      return false;
    }
    return true;
  }

  bool File::makeDir(const String& dir_name)
  {
    std::error_code ec;
    fs::create_directories(to_path(dir_name), ec);
    // create_directories returns false if the directory already existed, but that is not an error for us
    return !ec;
  }

  bool File::removeDirRecursively(const String& dir_name)
  {
    std::error_code ec;
    fs::remove_all(to_path(dir_name), ec);
    if (ec)
    {
      std::cerr << "Could not remove directory " << dir_name << ": " << ec.message() << std::endl;
      return false;
    }
    return true;
  }

  String File::absolutePath(const String& file)
  {
    if (file.empty()) return fs::current_path().generic_string();
#ifdef OPENMS_WINDOWSPLATFORM
    // On Windows, paths starting with '/' are treated as absolute by Qt
    // but fs::absolute() prepends the current drive letter. Preserve Qt behavior.
    if (file[0] == '/') return file;
#endif
    return fs::absolute(to_path(file)).generic_string();
  }

  String File::basename(const String& file)
  { // using well-defined overflow of unsigned ints here if path separator is not found
    return file.substr(file.find_last_of("\\/") + 1);
  }

  String File::stemName(const String& file)
  {
    return FileHandler::stripExtension(basename(file));
  }

  String File::extension(const String& file)
  {
    String base = basename(file);
    String stem = FileHandler::stripExtension(base);
    if (stem.size() >= base.size())
    {
      return ""; // no extension (stripExtension returned the same or longer string)
    }
    return base.substr(stem.size()); // everything after the stem, including leading '.'
  }

  StringList File::listDirectories(const String& dir)
  {
    StringList result;
    auto dir_path = to_path(dir);
    std::error_code ec;
    if (!fs::is_directory(dir_path, ec)) return result;

    for (fs::directory_iterator it(dir_path, ec), end; !ec && it != end; it.increment(ec))
    {
      if (it->is_directory(ec))
      {
        result.push_back(fs::absolute(it->path(), ec).generic_string());
      }
    }
    if (ec) return {};
    std::sort(result.begin(), result.end());
    return result;
  }

  String File::path(const String& file)
  {
    size_t pos = file.find_last_of("\\/");
    // do NOT return an empty string, because this leads to issues when in generic code you do:
    // String new_path = path("a.txt") + '/' + basename("a.txt");
    // , as this would lead to "/a.txt", i.e. create a wrong absolute path from a relative name
    String no_path = ".";
    return pos == string::npos ? no_path : file.substr(0, pos);
  }

  bool File::readable(const String& file)
  {
    auto p = to_path(file);
    std::error_code ec;
    if (!fs::exists(p, ec)) return false;
#ifdef OPENMS_WINDOWSPLATFORM
    return _access_s(file.c_str(), 4) == 0; // 4 = read permission
#else
    return access(file.c_str(), R_OK) == 0;
#endif
  }

  bool File::writable(const String& file)
  {
    auto p = to_path(file);
    std::error_code ec;

    if (fs::exists(p, ec))
    {
#ifdef OPENMS_WINDOWSPLATFORM
      return _access_s(file.c_str(), 2) == 0; // 2 = write permission
#else
      return access(file.c_str(), W_OK) == 0;
#endif
    }
    else
    {
      // File does not exist: probe by trying to create it
      std::ofstream f(file.c_str());
      bool ok = f.is_open() && f.good();
      f.close();
      if (ok)
      {
        std::remove(file.c_str());
      }
      return ok;
    }
  }

  String File::find(const String& filename, StringList directories)
  {
    // maybe we do not need to do anything?!
    // This check is required since calling File::find(File::find("CHEMISTRY/unimod.xml")) will otherwise fail
    // because the outer call receives an absolute path already
    if (exists(filename))
    {
      return filename;
    }
    String filename_new = filename;

    // empty string cannot be found, so throw Exception.
    // The code below would return success on empty string, since a path is prepended and thus the location exists
    if (filename_new.trim().empty())
    {
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }
    //add data dir in OpenMS data path
    directories.push_back(getOpenMSDataPath());

    //add path suffix to all specified directories
    String path = File::path(filename);
    if (!path.empty())
    {
      for (String& str : directories)
      {
        str.ensureLastChar('/');
        str += path;
      }
      filename_new = File::basename(filename);
    }

    //look up file
    for (const String& str : directories)
    {
      String loc = str;
      loc.ensureLastChar('/');
      loc = loc + filename_new;

      if (exists(loc))
      {
        return to_path(loc).lexically_normal().generic_string();
      }
    }

    //if the file was not found, throw an exception
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }

  bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
  {
    output.clear();
    auto dir_path = to_path(dir);
    std::error_code ec;
    if (!fs::is_directory(dir_path, ec)) return false;

    for (const auto& entry : fs::directory_iterator(dir_path, ec))
    {
      if (!entry.is_regular_file()) continue;
      String fname = entry.path().filename().string();
#ifdef OPENMS_WINDOWSPLATFORM
      if (!PathMatchSpecA(fname.c_str(), file_pattern.c_str())) continue;
#else
      if (fnmatch(file_pattern.c_str(), fname.c_str(), 0) != 0) continue;
#endif
      output.push_back(full_path ? entry.path().generic_string() : fname);
    }
    std::sort(output.begin(), output.end());
    return !output.empty();
  }

  String File::findDoc(const String& filename)
  {
    StringList search_dirs;
    search_dirs.push_back(String(OPENMS_BINARY_PATH) + "/../../doc/");
    // source path is host/openms so doc is ../doc
    search_dirs.push_back(String(OPENMS_SOURCE_PATH) + "/../../doc/");
    search_dirs.push_back(getOpenMSDataPath() + "/../../doc/");
    search_dirs.push_back(OPENMS_DOC_PATH);
    search_dirs.push_back(OPENMS_INSTALL_DOC_PATH);

    // needed for OpenMS Mac OS X packages where documentation is stored in <package-root>/Documentation
#if defined(__APPLE__)
    search_dirs.push_back(String(OPENMS_BINARY_PATH) + "/Documentation/");
    search_dirs.push_back(String(OPENMS_SOURCE_PATH) + "/Documentation/");
    search_dirs.push_back(getOpenMSDataPath() + "/../../Documentation/");
#endif

    return File::find(filename, search_dirs);
  }

  String File::getUniqueName(bool include_hostname)
  {
    DateTime now = DateTime::now();
    String pid;
#ifdef OPENMS_WINDOWSPLATFORM
    pid = (String)GetCurrentProcessId();
#else
    pid = (String)getpid();
#endif
    static std::atomic_int number = 0;
    String hostname_str;
    if (include_hostname)
    {
      char hbuf[256] = {};
#ifdef OPENMS_WINDOWSPLATFORM
      DWORD sz = sizeof(hbuf);
      if (!GetComputerNameA(hbuf, &sz))
      {
        OPENMS_LOG_WARN << "Warning: GetComputerNameA failed (error " << GetLastError() << ")" << std::endl;
        hbuf[0] = '\0';
      }
#else
      if (gethostname(hbuf, sizeof(hbuf)) != 0)
      {
        OPENMS_LOG_WARN << "Warning: gethostname failed (errno " << errno << ")" << std::endl;
        hbuf[0] = '\0';
      }
#endif
      if (hbuf[0] != '\0')
      {
        hostname_str = String(hbuf) + "_";
      }
    }
    return now.getDate().remove('-') + "_" + now.getTime().remove(':') + "_" + hostname_str + pid + "_" + (++number);
  }

  String File::getOpenMSDataPath()
  {
    // Use immediately evaluated lambda to protect static variable from concurrent access.
    static const String path = [&]() -> String {
      String path;
      bool path_checked = false;

      String found_path_from;
      bool from_env(false);
      if (getenv("OPENMS_DATA_PATH") != nullptr)
      {
        path = getenv("OPENMS_DATA_PATH");
        from_env = true;
        path_checked = isOpenMSDataPath_(path);
        if (path_checked)
        {
          found_path_from = "OPENMS_DATA_PATH (environment)";
        }
      }

      // probe the install path
      if (!path_checked)
      {
        path = OPENMS_INSTALL_DATA_PATH;
        path_checked = isOpenMSDataPath_(path);
        if (path_checked)
        {
          found_path_from = "OPENMS_INSTALL_DATA_PATH (compiled)";
        }
      }

      // probe the OPENMS_DATA_PATH macro
      if (!path_checked)
      {
        path = OPENMS_DATA_PATH;
        path_checked = isOpenMSDataPath_(path);
        if (path_checked) found_path_from = "OPENMS_DATA_PATH (compiled)";
      }

  #if defined(__APPLE__)
      // try to find it relative to the executable in the bundle (e.g. TOPPView)
      if (!path_checked)
      {
        path = getExecutablePath() + "../../../share/OpenMS";
        path_checked = isOpenMSDataPath_(path);
        if (path_checked) found_path_from = "bundle path (run time)";
      }
  #endif

      // On Linux and Apple check relative from the executable
      if (!path_checked)
      {
        path = getExecutablePath() + "../share/OpenMS";
        path_checked = isOpenMSDataPath_(path);
        if (path_checked)
        {
          found_path_from = "tool path (run time)";
        }
      }

      // make its a proper path:
      path = path.substitute("\\", "/").ensureLastChar('/').chop(1);

      if (!path_checked) // - now we're in big trouble as './share' is not were its supposed to be...
      { // - do NOT use OPENMS_LOG_ERROR or similar for the messages below! (it might not even usable at this point)
        std::cerr << "OpenMS FATAL ERROR!\n  Cannot find shared data! OpenMS cannot function without it!\n";
        if (from_env)
        {
          String p = getenv("OPENMS_DATA_PATH");
          std::cerr << "  The environment variable 'OPENMS_DATA_PATH' currently points to '" << p << "', which is incorrect!\n";
        }
  #ifdef OPENMS_WINDOWSPLATFORM
        String share_dir = R"(c:\Program Files\OpenMS\share\OpenMS)";
  #else
        String share_dir = "/usr/share/OpenMS";
  #endif
        std::cerr << "  To resolve this, set the environment variable 'OPENMS_DATA_PATH' to the OpenMS share directory (e.g., '" + share_dir + "').\n";
        std::cerr << "Exiting now.\n";
        exit(1);
      }
      return path;
    }();

    return path;
  }

  bool File::isOpenMSDataPath_(const String& path)
  {
    bool found = exists(path + "/CHEMISTRY/unimod.xml");
    return found;
  }

  bool File::isDirectory(const String& path)
  {
    return fs::is_directory(to_path(path));
  }

  String File::getTempDirectory()
  {
    Param p = getSystemParameters();
    String dir;
    if (getenv("OPENMS_TMPDIR") != nullptr)
    {
      dir = getenv("OPENMS_TMPDIR");
    }
    else if (p.exists("temp_dir") && !String(p.getValue("temp_dir").toString()).trim().empty())
    {
      dir = p.getValue("temp_dir").toString();
    }
    else
    {
      dir = fs::temp_directory_path().generic_string();
    }
    return dir;
  }

  /// The current OpenMS user data path (for result files)
  String File::getUserDirectory()
  {
    Param p = getSystemParameters();
    String dir;
    if (getenv("OPENMS_HOME_PATH") != nullptr)
    {
      dir = getenv("OPENMS_HOME_PATH");
    }
    else if (p.exists("home_dir") && !String(p.getValue("home_dir").toString()).trim().empty())
    {
      dir = p.getValue("home_dir").toString();
    }
    else
    {
#ifdef OPENMS_WINDOWSPLATFORM
      const char* home = getenv("USERPROFILE");
#else
      const char* home = getenv("HOME");
#endif
      dir = home ? String(home).substitute('\\', '/') : String(".");
    }
    dir.ensureLastChar('/');
    return dir;
  }

  String File::findDatabase(const String& db_name)
  {
    Param sys_p = getSystemParameters();
    String full_db_name;
    try
    {
      full_db_name = find(db_name, ListUtils::toStringList<std::string>(sys_p.getValue("id_db_dir")));
      OPENMS_LOG_INFO << "Augmenting database name '" << db_name << "' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '" << full_db_name << "'\n";
    }
    catch (Exception::FileNotFound& e)
    {
      OPENMS_LOG_ERROR << "Input database '" + db_name + "' not found (" << e.what() << "). Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!\n";
      throw;
    }

    return full_db_name;
  }

  String File::getOpenMSHomePath()
  {
    String home_path;
    // set path where OpenMS.ini is found from environment or use default
    if (getenv("OPENMS_HOME_PATH") != nullptr)
    {
      home_path = getenv("OPENMS_HOME_PATH");
    }
    else
    {
#ifdef OPENMS_WINDOWSPLATFORM
      const char* home = getenv("USERPROFILE");
#else
      const char* home = getenv("HOME");
#endif
      home_path = home ? String(home).substitute('\\', '/') : String(".");
    }
    return home_path;
  }

  Param File::getSystemParameters()
  {
    String home_path = File::getOpenMSHomePath();
    String filename;
    //Comply with https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html on unix identifying systems
    #ifdef __unix__
      if (getenv("XDG_CONFIG_HOME"))
      {
        filename = String(getenv("XDG_CONFIG_HOME")) + "/OpenMS/OpenMS.ini";
      }
      else
      {
        filename = File::getOpenMSHomePath() + "/.config/OpenMS/OpenMS.ini";
      }
    #else
      filename = home_path + "/.OpenMS/OpenMS.ini";
    #endif

    Param p;
    if (!File::readable(filename)) // no file, lets keep it that way
    {
      p = getSystemParameterDefaults_();
    }
    else
    {
      ParamXMLFile paramFile;
      paramFile.load(filename, p);

      // check version
      if (!p.exists("version") || (p.getValue("version") != VersionInfo::getVersion()))
      {
        if (!p.exists("version"))
        {
          OPENMS_LOG_WARN << "Broken file '" << filename << "' discovered. The 'version' tag is missing.\n";
        }
        else // old version
        {
          OPENMS_LOG_WARN << "File '" << filename << "' is deprecated.\n";
        }
        OPENMS_LOG_WARN << "Updating missing/wrong entries in '" << filename << "' with defaults!\n";
        Param p_new = getSystemParameterDefaults_();
        p.setValue("version", VersionInfo::getVersion()); // update old version, such that p_new:version does not get overwritten during update()
        p_new.update(p);
        // no new version is stored
      }
    }
    return p;
  }

  Param File::getSystemParameterDefaults_()
  {
    Param p;
    p.setValue("version", VersionInfo::getVersion());
    p.setValue("home_dir", ""); // only active when user enters something in this value
    p.setValue("temp_dir", ""); // only active when user enters something in this value
    p.setValue("id_db_dir", std::vector<std::string>(),
               String("Default directory for FASTA and psq files used as databased for id engines. ") + \
               "This allows you to specify just the filename of the DB in the " + \
               "respective TOPP tool, and the database will be searched in the directories specified here " + \
               ""); // only active when user enters something in this value
    p.setValue("threads", 1);

    return p;
  }

#ifdef OPENMS_WINDOWSPLATFORM
  StringList File::executableExtensions_(const String& ext)
  {
    // check if content of env-var %PATHEXT% makes sense
    StringList exts;
    ext.split(';', exts);
    // sanity check
    if (ListUtils::contains(exts, ".exe", ListUtils::CASE::INSENSITIVE)) return exts;
    // .. use fallback otherwise
    else return {".exe", ".bat" };
  }
#endif

  StringList File::getPathLocations(const String& path)
  {
    // split by ":" or ";", depending on platform
    StringList paths;
#ifdef OPENMS_WINDOWSPLATFORM
    path.split(';', paths);
#else
    path.split(':', paths);
#endif
    // ensure it ends with '/'
    for (String& p : paths) p.substitute('\\', '/').ensureLastChar('/');
    return paths;
  }

  bool File::findExecutable(OpenMS::String& exe_filename)
  {
    if (exists(exe_filename) && !isDirectory(exe_filename))
    {
      return true;
    }
    StringList paths = getPathLocations();
    StringList exe_filenames = { exe_filename };
#ifdef OPENMS_WINDOWSPLATFORM
    // try extensions like .exe on Windows
    if (!exe_filename.has('.'))
    {
      StringList exts = executableExtensions_();
      for (String& ext : exts) ext = exe_filename + ext;
      exe_filenames = exts;
    }
#endif
    // try all filenames (on Windows its potentially more than one) in each path...
    for (const String& p : paths)
    {
      for (const String& fn : exe_filenames)
      {
        if (exists(p + fn) && !isDirectory(p + fn))
        {
          exe_filename = p + fn;
          return true;
        }
      }
    }
    return false;
  }

  String File::findSiblingTOPPExecutable(const OpenMS::String& toolName)
  {
    // we first try the executablePath
    String exec = File::getExecutablePath() + toolName;

#if OPENMS_WINDOWSPLATFORM
    if (!exec.hasSuffix(".exe")) exec += ".exe";
#endif

    if (File::exists(exec))
    {
      return exec;
    }
#if defined(__APPLE__)
    // check if we are in one of the bundles (only built, not installed)
    exec = File::getExecutablePath() + "../../../" + toolName;
    if (File::exists(exec)) return exec;

    // check if we are in one of the bundles in an installed bundle (old bundles)
    exec = File::getExecutablePath() + "../../../TOPP/" + toolName;
    if (File::exists(exec)) return exec;

    // check if we are in one of the bundles in an installed bundle (new bundles)
    exec = File::getExecutablePath() + "../../../bin/" + toolName;
    if (File::exists(exec)) return exec;
#endif
    // TODO(aiche): probe in PATH

    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, toolName);
  }

  String File::getTemporaryFile(const String& alternative_file)
  {
    // take no action
    if (!alternative_file.empty())
    {
      return alternative_file;
    }
    // create temporary (and schedule for deletion)
    return temporary_files_.newFile();
  }


  File::TemporaryFiles_::TemporaryFiles_()
    : filenames_()
  {
  }

  String File::TemporaryFiles_::newFile()
  {
    String s = getTempDirectory().ensureLastChar('/') + getUniqueName();
    std::lock_guard<std::mutex> _(mtx_);
    filenames_.push_back(s);
    // do NOT return filenames_.back() by ref, since another thread might resize the vector and invalidate the reference!
    return s; // uses RVO, so its efficient
  }

  File::TemporaryFiles_::~TemporaryFiles_()
  {
    std::lock_guard<std::mutex> _(mtx_);
    for (Size i = 0; i < filenames_.size(); ++i)
    {
      if (File::exists(filenames_[i]) && !File::remove(filenames_[i]))
      {
        std::cerr << "Warning: unable to remove temporary file '" << filenames_[i] << "'" << std::endl;
      }
    }
  }

  File::MatchingFileListsStatus File::validateMatchingFileNames(const StringList& sl1,
                                                        const StringList& sl2,
                                                        bool basename,
                                                        bool ignore_extension)
  {
      // Different counts means different sets
      if (sl1.size() != sl2.size())
      {
          return MatchingFileListsStatus::SET_MISMATCH;
      }

      set<String> sl1_set;
      set<String> sl2_set;
      bool different_name_at_index = false;

      // Process and compare each filename
      for (size_t i = 0; i != sl1.size(); ++i)
      {
          String sl1_name = sl1[i];
          String sl2_name = sl2[i];

          if (basename)
          {
              sl1_name = File::basename(sl1_name);
              sl2_name = File::basename(sl2_name);
          }

          if (ignore_extension)
          {
              sl1_name = FileHandler::stripExtension(sl1_name);
              sl2_name = FileHandler::stripExtension(sl2_name);
          }

          sl1_set.insert(sl1_name);
          sl2_set.insert(sl2_name);

          if (sl1_name != sl2_name)
          {
              different_name_at_index = true;
          }
      }

      bool same_set = (sl1_set == sl2_set);

      // Check if it's an order mismatch or complete mismatch
      if (same_set)
      {
          return different_name_at_index ?
                MatchingFileListsStatus::ORDER_MISMATCH :
                MatchingFileListsStatus::MATCH;
      }

      return MatchingFileListsStatus::SET_MISMATCH;
  }

  File::TemporaryFiles_ File::temporary_files_;

} // namespace OpenMS
