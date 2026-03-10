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

#include <atomic>
#include <cerrno>
#include <filesystem>
#include <chrono>
#include <cstdlib>
#include <algorithm>
#include <fstream>
#include <regex>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h> // for GetCurrentProcessId() && GetModuleFileName() && GetComputerNameA()
#include <Shlobj.h> // for SHGetFolderPath
#include <io.h> // for _access_s
#endif

#ifdef OPENMS_HAS_UNISTD_H
#include <unistd.h> // for readLink() and getpid()
#endif

#ifdef __APPLE__
#include <mach-o/dyld.h>
#endif

using namespace std;

namespace OpenMS
{
  namespace
  {
    // Helper function to get home directory cross-platform
    String getHomePath()
    {
#ifdef OPENMS_WINDOWSPLATFORM
      // Use SHGetFolderPath for Windows
      char path[MAX_PATH];
      if (SUCCEEDED(SHGetFolderPathA(NULL, CSIDL_PROFILE, NULL, 0, path)))
      {
        return String(path);
      }
      // Fallback to environment variables
      const char* homeDrive = getenv("HOMEDRIVE");
      const char* homePath = getenv("HOMEPATH");
      if (homeDrive && homePath)
      {
        return String(homeDrive) + String(homePath);
      }
      const char* userProfile = getenv("USERPROFILE");
      if (userProfile)
      {
        return String(userProfile);
      }
      return String("C:\\");
#else
      // Unix-like systems (Linux, macOS)
      const char* home = getenv("HOME");
      if (home)
      {
        return String(home);
      }
      return String("/tmp");
#endif
    }
  }

  File::TempDir::TempDir(bool keep_dir)
    : keep_dir_(keep_dir)
  {
    temp_dir_ = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'" << std::endl;
    std::filesystem::create_directories(static_cast<std::string>(temp_dir_));
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
    std::error_code ec;
    return std::filesystem::exists(static_cast<std::string>(file), ec);
  }

  bool File::empty(const String& file)
  {
    std::error_code ec;
    if (!std::filesystem::exists(static_cast<std::string>(file), ec))
    {
      return true;
    }
    return std::filesystem::file_size(static_cast<std::string>(file), ec) == 0;
  }

  bool File::executable(const String& file)
  {
    std::error_code ec;
    if (!std::filesystem::exists(static_cast<std::string>(file), ec))
    {
      return false;
    }

    [[maybe_unused]] auto perms = std::filesystem::status(static_cast<std::string>(file), ec).permissions();
    if (ec) return false;
    
    // Check if executable bit is set (on Unix-like systems)
#ifdef OPENMS_WINDOWSPLATFORM
    // On Windows, check if it's a .exe, .bat, .cmd file or similar
    std::filesystem::path p(static_cast<std::string>(file));
    String ext = p.extension().string();
    std::transform(ext.begin(), ext.end(), ext.begin(), ::tolower);
    return (ext == ".exe" || ext == ".bat" || ext == ".cmd" || ext == ".com");
#else
    return (perms & std::filesystem::perms::owner_exec) != std::filesystem::perms::none ||
           (perms & std::filesystem::perms::group_exec) != std::filesystem::perms::none ||
           (perms & std::filesystem::perms::others_exec) != std::filesystem::perms::none;
#endif
  }

  UInt64 File::fileSize(const String& file)
  {
    if (!File::exists(file)) return -1;

    std::error_code ec;
    auto size = std::filesystem::file_size(static_cast<std::string>(file), ec);
    return ec ? 0 : size;
  }

  bool File::rename(const String& from, const String& to, bool overwrite_existing, bool verbose)
  {
    // check for equality
    std::error_code ec1, ec2;
    auto from_canonical = std::filesystem::canonical(static_cast<std::string>(from), ec1);
    auto to_canonical = std::filesystem::canonical(static_cast<std::string>(to), ec2);
    
    // If both files exist and are the same, do nothing
    if (!ec1 && !ec2 && from_canonical == to_canonical)
    { // same file; no need to to anything
      return true;
    }

    // existing file? Qt won't overwrite, so try to remove it:
    if (overwrite_existing && exists(to) && !remove(to))
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not overwrite existing file '" << to << "'\n";
      }
      return false;
    }
    // move the file to the actual destination:
    std::error_code ec;
    std::filesystem::rename(std::filesystem::path(static_cast<std::string>(from)), std::filesystem::path(static_cast<std::string>(to)), ec);
    if (ec)
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "Error: Could not move '" << from << "' to '" << to << "': " << ec.message() << "\n";
      }
      return false;
    }
    return true;
  }

  // https://stackoverflow.com/questions/2536524/copy-directory-using-qt
  bool File::copyDirRecursively(const String& from_dir, const String& to_dir, File::CopyOptions option)
  {
    namespace fs = std::filesystem;
    fs::path source_path(static_cast<std::string>(from_dir));
    fs::path target_path(static_cast<std::string>(to_dir));

    try
    {
      // Check source exists and is a directory
      if (!fs::exists(source_path) || !fs::is_directory(source_path))
      {
        OPENMS_LOG_ERROR << "Error: Source directory '" << from_dir << "' does not exist or is not a directory." << std::endl;
        return false;
      }

      // If target exists, ensure it is not the same as source
      std::error_code ec_equiv;
      if (fs::exists(target_path) && fs::equivalent(source_path, target_path, ec_equiv))
      {
        OPENMS_LOG_ERROR << "Error: Could not copy  " << from_dir << " to " << to_dir << ". Same path given." << std::endl;
        return false;
      }

      // Create target if not present (parents included)
      if (!fs::exists(target_path))
      {
        fs::create_directories(target_path);
      }

      // copy folder recursively
      for (const auto& entry : fs::directory_iterator(source_path))
      {
        const fs::path& entry_path = entry.path();
        fs::path target_entry = target_path / entry_path.filename();

        if (entry.is_directory())
        {
          if (!copyDirRecursively(String(entry_path.string()), String(target_entry.string()), option))
          {
            return false;
          }
        }
        else
        {
          if (fs::exists(target_entry))
          {
            switch (option)
            {
              case CopyOptions::CANCEL:
                return false;
              case CopyOptions::SKIP:
                OPENMS_LOG_WARN << "The file " << entry_path.filename().string() << " was skipped." << std::endl;
                continue;
              case CopyOptions::OVERWRITE:
              {
                std::error_code ec_rm;
                fs::remove(target_entry, ec_rm); // ignore error; best-effort
                break;
              }
            }
          }

          std::error_code ec_cp;
          fs::copy_file(entry_path, target_entry, ec_cp);
          if (ec_cp)
          {
            OPENMS_LOG_ERROR << "Error: Could not copy file '" << entry_path.string() << "' to '" << target_entry.string() << "': " << ec_cp.message() << std::endl;
            return false;
          }
        }
      }
      return true;
    }
    catch (const std::filesystem::filesystem_error& /*e*/)
    {
      return false;
    }
  }

  bool File::copy(const String& from, const String& to)
  {
    std::error_code ec;
    bool success = std::filesystem::copy_file(std::filesystem::path(static_cast<std::string>(from)), std::filesystem::path(static_cast<std::string>(to)), ec);
    return success && !ec;
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
    try 
    {
      std::filesystem::path dir_path(static_cast<std::string>(dir_name));
      if (std::filesystem::exists(dir_path) && std::filesystem::is_directory(dir_path))
      {
        return std::filesystem::remove_all(dir_path) > 0;
      }
      return true; // Directory doesn't exist, consider it removed
    }
    catch (const std::filesystem::filesystem_error& /*e*/)
    {
      return false;
    }
  }

  bool File::makeDir(const String& dir_name)
  {
    try
    {
      const std::string s = static_cast<std::string>(dir_name);
      std::error_code ec;
      // If it already exists and is a directory, that's success (Qt QDir::mkpath behavior)
      if (std::filesystem::exists(s, ec) && std::filesystem::is_directory(s, ec))
      {
        return true;
      }
      // If a non-directory entry exists at the path, fail
      if (std::filesystem::exists(s, ec) && !std::filesystem::is_directory(s, ec))
      {
        return false;
      }
      // Otherwise, attempt to create the directory (and parents)
      return std::filesystem::create_directories(s);
    }
    catch (const std::filesystem::filesystem_error& /*e*/)
    {
      return false;
    }
  }

  bool File::removeDirRecursively(const String& dir_name)
  {
    try
    {
      std::filesystem::path dir_path(static_cast<std::string>(dir_name));
      if (std::filesystem::exists(dir_path) && std::filesystem::is_directory(dir_path))
      {
        return std::filesystem::remove_all(dir_path) > 0;
      }
      return true; // Directory doesn't exist, consider it removed
    }
    catch (const std::filesystem::filesystem_error& /*e*/)
    {
      return false;
    }
  }

  String File::absolutePath(const String& file)
  {
    std::error_code ec;
    auto abs_path = std::filesystem::absolute(static_cast<std::string>(file), ec);
    return ec ? static_cast<std::string>(file) : abs_path.string();
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
#ifdef OPENMS_WINDOWSPLATFORM
    return _access_s(file.c_str(), 4) == 0; // 4 = read permission
#else
    return access(file.c_str(), R_OK) == 0;
#endif
  }

  bool File::writable(const String& file)
  {
#ifdef OPENMS_WINDOWSPLATFORM
    if (_access_s(file.c_str(), 0) == 0) // file exists
    {
      return _access_s(file.c_str(), 2) == 0; // 2 = write permission
    }
#else
    if (access(file.c_str(), F_OK) == 0) // file exists
    {
      return access(file.c_str(), W_OK) == 0;
    }
#endif
    // File doesn't exist, try to create it temporarily to test if we can write
    std::ofstream test_file(std::string(file));
    bool is_writable = test_file.good();
    test_file.close();
    if (is_writable)
    {
      std::error_code ec;
      std::filesystem::remove(std::string(file), ec);
    }
    return is_writable;
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
        std::filesystem::path clean_path = std::filesystem::path(static_cast<std::string>(loc)).lexically_normal();
        return String(clean_path.string());
      }
    }

    //if the file was not found, throw an exception
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }

  bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
  {
    output.clear();
    
    try
    {
      std::filesystem::path dir_path(static_cast<std::string>(dir));
      
      if (!std::filesystem::exists(dir_path) || !std::filesystem::is_directory(dir_path))
      {
        return false;
      }

      std::vector<std::filesystem::path> matching_files;
      
      // Convert glob pattern to regex: escape metacharacters, then convert * and ?
      std::string pattern_str = static_cast<std::string>(file_pattern);
      std::string regex_pattern;
      regex_pattern.reserve(pattern_str.size() * 2);
      for (char c : pattern_str)
      {
        switch (c)
        {
          case '*': regex_pattern += ".*"; break;
          case '?': regex_pattern += '.'; break;
          // Escape regex metacharacters
          case '.': case '(': case ')': case '[': case ']':
          case '{': case '}': case '+': case '^': case '$':
          case '|': case '\\':
            regex_pattern += '\\';
            regex_pattern += c;
            break;
          default:
            regex_pattern += c;
            break;
        }
      }

      std::regex pattern_regex(regex_pattern);
      
      // Iterate through directory entries
      for (const auto& entry : std::filesystem::directory_iterator(dir_path))
      {
        if (entry.is_regular_file())
        {
          std::string filename = entry.path().filename().string();
          if (std::regex_match(filename, pattern_regex))
          {
            matching_files.push_back(entry.path());
          }
        }
      }
      
      if (matching_files.empty())
      {
        return false;
      }
      
      // Sort files by name (similar to QDir::Name)
      std::sort(matching_files.begin(), matching_files.end(),
                [](const std::filesystem::path& a, const std::filesystem::path& b) {
                  return a.filename().string() < b.filename().string();
                });
      
      output.resize(matching_files.size());
      for (size_t i = 0; i < matching_files.size(); ++i)
      {
        output[i] = full_path ? String(matching_files[i].string()) : String(matching_files[i].filename().string());
      }
      
      return true;
    }
    catch (const std::exception& /*e*/)
    {
      return false;
    }
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
    std::error_code ec;
    return std::filesystem::is_directory(static_cast<std::string>(path), ec);
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
      dir = String(std::filesystem::temp_directory_path().string());
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
      dir = getHomePath();
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
      home_path = getHomePath();
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
