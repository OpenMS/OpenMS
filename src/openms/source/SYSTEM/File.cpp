// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <OpenMS/SYSTEM/NetworkGetRequest.h>

#include <atomic>
#include <fstream>
#include <filesystem>
#include <cstring>
#include <regex>
#include <algorithm>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h> // for GetCurrentProcessId() && GetModuleFileName()
#include <direct.h>  // for _getcwd, _mkdir
#else
#include <unistd.h>  // for readlink, getpid, gethostname, getcwd
#include <limits.h>  // for PATH_MAX, HOST_NAME_MAX
#include <sys/stat.h> // for mkdir
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace std;

namespace OpenMS
{

  File::TempDir::TempDir(bool keep_dir)
    : keep_dir_(keep_dir)
  {
    temp_dir_ = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'" << std::endl;
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(temp_dir_.c_str()), ec);
    if (ec)
    {
      String error_msg = "Failed to create temporary directory '" + temp_dir_ + "': " + ec.message();
      OPENMS_LOG_ERROR << error_msg << std::endl;
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error_msg);
    }
  };

  File::TempDir::~TempDir()
  {
    if (keep_dir_)
    {
      OPENMS_LOG_DEBUG << "Keeping temporary files in directory '" << temp_dir_ << std::endl;
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
    std::error_code ec;
    return std::filesystem::exists(std::filesystem::path(file.c_str()), ec);
  }

  bool File::empty(const String& file)
  {
    std::error_code ec;
    std::filesystem::path p(file.c_str());
    if (!std::filesystem::exists(p, ec)) return true;
    return std::filesystem::file_size(p, ec) == 0;
  }

  bool File::executable(const String& file)
  {
    std::error_code ec;
    if (!std::filesystem::exists(std::filesystem::path(file.c_str()), ec)) return false;
#ifdef OPENMS_WINDOWSPLATFORM
    // On Windows, check if file ends with .exe, .bat, .cmd, or .com
    String lower = file;
    lower.toLower();
    return lower.hasSuffix(".exe") || lower.hasSuffix(".bat") ||
           lower.hasSuffix(".cmd") || lower.hasSuffix(".com");
#else
    // On Unix, check execute permission
    struct stat st;
    if (stat(file.c_str(), &st) != 0) return false;
    return (st.st_mode & S_IXUSR) != 0;
#endif
  }

  UInt64 File::fileSize(const String& file)
  {
    std::error_code ec;
    std::filesystem::path p(file.c_str());
    if (!std::filesystem::exists(p, ec)) return -1;
    return std::filesystem::file_size(p, ec);
  }

  bool File::rename(const String& from, const String& to, bool overwrite_existing, bool verbose)
  {
    std::error_code ec;
    std::filesystem::path from_path(from.c_str());
    std::filesystem::path to_path(to.c_str());

    // check for equality using canonical paths
    std::filesystem::path from_canonical = std::filesystem::canonical(from_path, ec);
    if (ec) from_canonical = from_path; // If canonical fails, use original path
    ec.clear();
    std::filesystem::path to_canonical = std::filesystem::canonical(to_path, ec);
    if (ec) to_canonical = to_path; // If canonical fails, use original path

    if (from_canonical == to_canonical)
    { // same file; no need to do anything
      return true;
    }

    // existing file? std::filesystem won't overwrite, so try to remove it:
    if (overwrite_existing && exists(to) && !remove(to))
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not overwrite existing file '" << to << "'\n";
      }
      return false;
    }
    // move the file to the actual destination:
    ec.clear();
    std::filesystem::rename(from_path, to_path, ec);
    if (ec)
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not move '" << from << "' to '" << to << "'\n";
      }
      return false;
    }
    return true;
  }

  // https://stackoverflow.com/questions/2536524/copy-directory-using-qt
  bool File::copyDirRecursively(const String& from_dir, const String& to_dir, File::CopyOptions option)
  {
    std::error_code ec;
    std::filesystem::path source_path(from_dir.c_str());
    std::filesystem::path target_path(to_dir.c_str());

    // Get canonical paths to check for same directory
    std::filesystem::path canonical_source = std::filesystem::canonical(source_path, ec);
    if (ec) canonical_source = source_path;
    ec.clear();
    std::filesystem::path canonical_target = std::filesystem::canonical(target_path, ec);
    if (ec) canonical_target = target_path;
    ec.clear();

    // check canonical path
    if (canonical_source == canonical_target)
    {
      OPENMS_LOG_ERROR << "Error: Could not copy " << from_dir << " to " << to_dir << ". Same path given." << std::endl;
      return false;
    }

    // make directory if not present
    if (!std::filesystem::exists(target_path, ec))
    {
      std::filesystem::create_directories(target_path, ec);
      if (ec)
      {
        OPENMS_LOG_ERROR << "Error: Could not create target directory " << to_dir << std::endl;
        return false;
      }
    }

    // copy folder recursively
    for (const auto& entry : std::filesystem::directory_iterator(source_path, ec))
    {
      std::string entry_name = entry.path().filename().string();

      // Skip . and .. (though directory_iterator should skip these automatically)
      if (entry_name == "." || entry_name == "..")
      {
        continue;
      }

      std::filesystem::path target_entry = target_path / entry_name;

      if (entry.is_directory())
      {
        if (!copyDirRecursively(entry.path().string(), target_entry.string(), option))
        {
          return false;
        }
      }
      else
      {
        // Check if target file exists
        if (std::filesystem::exists(target_entry, ec))
        {
          switch (option)
          {
            case CopyOptions::CANCEL:
              return false;
            case CopyOptions::SKIP:
              OPENMS_LOG_WARN << "The file " << entry_name << " was skipped." << std::endl;
              continue;
            case CopyOptions::OVERWRITE:
              std::filesystem::remove(target_entry, ec);
              break;
          }
        }

        // Copy the file
        std::filesystem::copy_file(entry.path(), target_entry, ec);
        if (ec)
        {
          OPENMS_LOG_ERROR << "Error: Could not copy file " << entry.path().string() << " to " << target_entry.string() << std::endl;
          return false;
        }
      }
    }
    return true;
  }

  bool File::copy(const String& from, const String& to)
  {
    std::error_code ec;
    std::filesystem::copy_file(std::filesystem::path(from.c_str()),
                                std::filesystem::path(to.c_str()),
                                std::filesystem::copy_options::overwrite_existing,
                                ec);
    return !ec;
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
    std::filesystem::path dir_path(dir_name.c_str());

    if (!std::filesystem::exists(dir_path, ec))
    {
      return true; // Already doesn't exist
    }

    // Recursively remove directory contents
    for (const auto& entry : std::filesystem::directory_iterator(dir_path, ec))
    {
      if (entry.is_directory())
      {
        if (!removeDir(entry.path().string()))
        {
          return false;
        }
      }
      else
      {
        if (!std::filesystem::remove(entry.path(), ec))
        {
          return false;
        }
      }
    }

    // Remove the directory itself
    return std::filesystem::remove(dir_path, ec);
  }

  bool File::makeDir(const String& dir_name)
  {
    std::error_code ec;
    std::filesystem::create_directories(std::filesystem::path(dir_name.c_str()), ec);
    return !ec;
  }

  bool File::removeDirRecursively(const String& dir_name)
  {
    std::error_code ec;
    std::filesystem::path dir_path(dir_name.c_str());

    // Remove directory and all its contents recursively
    std::uintmax_t removed = std::filesystem::remove_all(dir_path, ec);

    if (ec)
    {
      OPENMS_LOG_WARN << "Could not remove directory " << dir_name << ": " << ec.message() << std::endl;
      return false;
    }

    return removed > 0;
  }

  String File::absolutePath(const String& file)
  {
    std::error_code ec;
    std::filesystem::path abs_path = std::filesystem::absolute(std::filesystem::path(file.c_str()), ec);
    if (ec) return file;
    return abs_path.string();
  }

  String File::basename(const String& file)
  { // using well-defined overflow of unsigned ints here if path separator is not found
    return file.substr(file.find_last_of("\\/") + 1);
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
    if (!exists(file)) return false;
#ifdef OPENMS_WINDOWSPLATFORM
    return _access(file.c_str(), 4) == 0; // 4 = read permission
#else
    return access(file.c_str(), R_OK) == 0;
#endif
  }

  bool File::writable(const String& file)
  {
    if (exists(file))
    {
#ifdef OPENMS_WINDOWSPLATFORM
      return _access(file.c_str(), 2) == 0; // 2 = write permission
#else
      return access(file.c_str(), W_OK) == 0;
#endif
    }
    else
    {
      // Try to create and remove the file to test writability
      std::ofstream f(file.c_str());
      bool writable = f.good();
      f.close();
      if (writable)
      {
        std::remove(file.c_str());
      }
      return writable;
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
        // Normalize path (remove . and .., convert separators)
        std::filesystem::path p(loc.c_str());
        return p.lexically_normal().string();
      }
    }

    //if the file was not found, throw an exception
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }

  bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
  {
    output.clear();

    std::error_code ec;
    std::filesystem::path dir_path(dir.c_str());

    if (!std::filesystem::is_directory(dir_path, ec))
    {
      return false;
    }

    // Convert wildcard pattern to regex
    // Simple pattern matching: * -> .*, ? -> .
    String pattern_regex = file_pattern;
    pattern_regex.substitute(".", "\\.");
    pattern_regex.substitute("*", ".*");
    pattern_regex.substitute("?", ".");
    std::regex pattern("^" + pattern_regex + "$");

    std::vector<std::filesystem::path> matches;
    for (const auto& entry : std::filesystem::directory_iterator(dir_path, ec))
    {
      if (!entry.is_regular_file()) continue;

      std::string filename = entry.path().filename().string();
      if (std::regex_match(filename, pattern))
      {
        matches.push_back(entry.path());
      }
    }

    if (matches.empty())
    {
      return false;
    }

    // Sort by name
    std::sort(matches.begin(), matches.end());

    // Fill output
    for (const auto& match : matches)
    {
      output.push_back(full_path ? match.string() : match.filename().string());
    }

    return true;
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

    String hostname;
    if (include_hostname)
    {
#ifdef OPENMS_WINDOWSPLATFORM
      char buffer[256];
      DWORD size = sizeof(buffer);
      if (GetComputerNameA(buffer, &size))
      {
        hostname = String(buffer) + "_";
      }
#else
      char buffer[256];
      if (gethostname(buffer, sizeof(buffer)) == 0)
      {
        hostname = String(buffer) + "_";
      }
#endif
    }

    return now.getDate().remove('-') + "_" + now.getTime().remove(':') + "_" + hostname + pid + "_" + (++number);
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
    std::error_code ec;
    return std::filesystem::is_directory(std::filesystem::path(path.c_str()), ec);
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
      std::error_code ec;
      dir = std::filesystem::temp_directory_path(ec).string();
      if (ec) dir = "/tmp"; // fallback
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
      if (!home) home = getenv("HOME");
#else
      const char* home = getenv("HOME");
#endif
      if (home)
      {
        dir = String(home);
      }
      else
      {
        dir = "."; // fallback to current directory
      }
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
      OPENMS_LOG_INFO << "Augmenting database name '" << db_name << "' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '" << full_db_name << "'" << std::endl;
    }
    catch (Exception::FileNotFound& e)
    {
      OPENMS_LOG_ERROR << "Input database '" + db_name + "' not found (" << e.what() << "). Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!" << std::endl;
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
      if (!home) home = getenv("HOME");
#else
      const char* home = getenv("HOME");
#endif
      if (home)
      {
        home_path = String(home);
      }
      else
      {
        home_path = "."; // fallback to current directory
      }
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
          OPENMS_LOG_WARN << "Broken file '" << filename << "' discovered. The 'version' tag is missing." << std::endl;
        }
        else // old version
        {
          OPENMS_LOG_WARN << "File '" << filename << "' is deprecated." << std::endl;
        }
        OPENMS_LOG_WARN << "Updating missing/wrong entries in '" << filename << "' with defaults!" << std::endl;
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

  // construct a filename. Add number if already exists.
  String saveFileName_(const String& url)
  {
    // Extract filename from URL path
    String basename = File::basename(url);

    // Remove query parameters if present
    if (basename.has('?'))
    {
      basename = basename.prefix('?');
    }

    if (basename.empty())
    {
      basename = "download";
    }

    // If file exists, add number suffix
    if (File::exists(basename))
    {
      String original_basename = basename;
      int i = 0;
      basename = original_basename + "." + String(i);
      while (File::exists(basename))
      {
        ++i;
        basename = original_basename + "." + String(i);
      }
    }

    return basename;
  }

// static
void File::download(const std::string& url, const std::string& download_folder)
{
  NetworkGetRequest query;
  query.setUrl(url);
  query.setTimeout(600000); // 10 minutes timeout
  query.run();

  if (!query.hasError())
  {
    String folder = download_folder.empty() ? "./" : download_folder;
    folder.ensureLastChar('/');
    String filename = folder + saveFileName_(url);

    // Write file
    std::ofstream file(filename.c_str(), std::ios::binary);
    if (!file.is_open())
    {
      throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                   "Cannot open file for writing: " + filename);
    }

    const std::vector<char>& data = query.getResponseBinary();
    file.write(data.data(), data.size());
    file.close();

    OPENMS_LOG_INFO << "Download of '" << url << "' successful." << endl;
    OPENMS_LOG_INFO << "Stored as '" << filename << "'." << endl;
  }
  else
  {
    String error = "Download of '" + url + "' failed!. Error: " + query.getErrorString() + '\n';
    throw Exception::IOException(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, error);
  }
}


} // namespace OpenMS
