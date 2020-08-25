// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/openms_data_path.h>

#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtNetwork/QHostInfo>

#ifdef OPENMS_WINDOWSPLATFORM
#include <Windows.h> // for GetCurrentProcessId() && GetModuleFileName()
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

  File::TempDir::TempDir(bool keep_dir)
    : keep_dir_(keep_dir)
  {
    temp_dir_ = File::getTempDirectory() + "/" + File::getUniqueName() + "/";
    OPENMS_LOG_DEBUG << "Creating temporary directory '" << temp_dir_ << "'" << std::endl;
    QDir d;
    d.mkpath(temp_dir_.toQString());
  };

  File::TempDir::~TempDir()
  {
    if (keep_dir_) {
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
    static String spath = "";
    static bool path_checked = false;

    // short route. Only inquire the path once. The result will be the same every time.
    if (path_checked) return spath;

    char path[1024];

#ifdef OPENMS_WINDOWSPLATFORM
    int size = sizeof(path);
    if (GetModuleFileName(NULL, path, size))
#elif  defined(__APPLE__)
    uint size = sizeof(path);
    if (_NSGetExecutablePath(path, &size) == 0)
#else // LINUX
    int size = sizeof(path);
    int ch = readlink("/proc/self/exe", path, size);
    if (ch != -1)
#endif
    {
      spath = File::path(String(path));
      if (File::exists(spath)) // check if directory exists
      {
        // ensure path ends with a "/", such that we can just write path + "ToolX", and to not worry about if its empty or a path.
        spath.ensureLastChar('/');
      }
      else
      {
        std::cerr << "Path extracted from Executable Path does not exist! Returning empty string!\n";
        spath = "";
      }
    }
    else
    {
      std::cerr << "Cannot get Executable Path! Not using a path prefix!\n";
    }

    path_checked = true; // enable short route for next run
    return spath;
  }

  bool File::exists(const String& file)
  {
    QFileInfo fi(file.toQString());
    return fi.exists();
  }

  bool File::empty(const String& file)
  {
    QFileInfo fi(file.toQString());
    return !fi.exists() || fi.size() == 0;
  }

  bool File::rename(const String& from, const String& to, bool overwrite_existing, bool verbose)
  {
    // check for equality
    if (QFileInfo(from.c_str()).canonicalFilePath() == QFileInfo(to.c_str()).canonicalFilePath())
    { // same file; no need to to anything
      return true;
    }

    // existing file? Qt won't overwrite, so try to remove it:
    if (overwrite_existing && exists(to) && !remove(to))
    {
      if (verbose) OPENMS_LOG_ERROR << "Error: Could not overwrite existing file '" << to << "'\n";
      return false;
    }
    // move the file to the actual destination:
    if (!QFile::rename(from.toQString(), to.toQString()))
    {
      if (verbose) OPENMS_LOG_ERROR << "Error: Could not move '" << from << "' to '" << to << "'\n";
      return false;
    }
    return true;
  }

  // https://stackoverflow.com/questions/2536524/copy-directory-using-qt
  bool File::copyDirRecursively(const QString &from_dir, const QString &to_dir, File::CopyOptions option)
  {
    QDir source_dir(from_dir);
    QDir target_dir(to_dir);

    QString canonical_source_dir = source_dir.canonicalPath();
    QString canonical_target_dir = target_dir.canonicalPath();
   
    // check canonical path  
    if (canonical_source_dir == canonical_target_dir)
    {
      OPENMS_LOG_ERROR << "Error: Could not copy  " << from_dir.toStdString() << " to " << to_dir.toStdString() << ". Same path given." << std::endl;;  
      return false;
    }

    // make directory if not present 
    if (!target_dir.exists())
    {
      target_dir.mkpath(to_dir);
    }
  
    // copy folder recurively
    QFileInfoList file_list = source_dir.entryInfoList();
    for (const QFileInfo& entry : file_list)   
    {
      if (entry.fileName() == "." || entry.fileName() == "..")
      {
        continue;
      }
      if (entry.isDir())
      {
        if (!copyDirRecursively(entry.filePath(), target_dir.filePath(entry.fileName()), option))
        {
          return false;
        }
      }
      else
      {
        if (target_dir.exists(entry.fileName()))
        {
          switch (option)
            {
              case CopyOptions::CANCEL: 
                return false;
              case CopyOptions::SKIP: 
                OPENMS_LOG_WARN << "The file " << entry.fileName().toStdString() << " was skipped." << std::endl; 
                continue;
              case CopyOptions::OVERWRITE:
                target_dir.remove(entry.fileName());
            }
        }
        if (!QFile::copy(entry.filePath(), target_dir.filePath(entry.fileName())))
        {
          return false;
        }
      }
    }
    return true;
  }

  bool File::remove(const String& file)
  {
    if (!exists(file))
      return true;

    if (std::remove(file.c_str()) != 0)
      return false;

    return true;
  }

  bool File::removeDir(const QString& dir_name)
  {
    bool result = true;
    QDir dir(dir_name);

    if (dir.exists(dir_name))
    {
      Q_FOREACH(QFileInfo info, dir.entryInfoList(QDir::NoDotAndDotDot | QDir::System | QDir::Hidden  | QDir::AllDirs | QDir::Files, QDir::DirsFirst))
        {
          if (info.isDir())
          {
            result = removeDir(info.absoluteFilePath());
          }
          else
          {
            result = QFile::remove(info.absoluteFilePath());
          }
          if (!result)
          {
            return result;
          }
        }
      result = dir.rmdir(dir_name);
    }
    return result;
  }

  bool File::removeDirRecursively(const String& dir_name)
  {
    bool fail = false;
    QString path = dir_name.toQString();
    QDir dir(path);
    QStringList files = dir.entryList(QDir::Files | QDir::NoDotAndDotDot);
    foreach(const QString &file_name, files)
    {
      if (!dir.remove(file_name))
      {
        OPENMS_LOG_WARN << "Could not remove file " << String(file_name) << "!" << std::endl;
        fail = true;
      }
    }
    QStringList contained_dirs = dir.entryList(QDir::Dirs | QDir::NoDotAndDotDot);
    foreach(const QString &contained_dir, contained_dirs)
    {
      if (!removeDirRecursively(path + QDir::separator() + contained_dir))
      {
        fail = true;
      }
    }

    QDir parent_dir(path);
    if (parent_dir.cdUp())
    {
      if (!parent_dir.rmdir(path))
      {
        std::cerr << "Could not remove directory " << String(dir.dirName()) << "!" << std::endl;
        fail = true;
      }
    }

    return !fail;
  }

  String File::absolutePath(const String& file)
  {
    QFileInfo fi(file.toQString());
    return fi.absoluteFilePath();
  }

  String File::basename(const String& file)
  {
    QFileInfo fi(file.toQString());
    return fi.fileName();
  }

  String File::path(const String& file)
  {
    QFileInfo fi(file.toQString());
    return fi.path();
  }

  bool File::readable(const String& file)
  {
    QFileInfo fi(file.toQString());
    return fi.exists() && fi.isReadable();
  }

  bool File::writable(const String& file)
  {
    QFileInfo fi(file.toQString());

    bool tmp = false;
    if (fi.exists())
    {
      tmp = fi.isWritable();
    }
    else
    {
      QFile f;
      f.setFileName(file.toQString());
      f.open(QIODevice::WriteOnly);
      tmp = f.isWritable();
      f.remove();
    }

    return tmp;
  }

  String File::find(const String& filename, StringList directories)
  {
    // maybe we do not need to do anything?!
    // This check is required since calling File::find(File::find("CHEMISTRY/Elements.xml")) will otherwise fail
    // because the outer call receives an absolute path already
    if (exists(filename)) return filename;

    String filename_new = filename;

    // empty string cannot be found, so throw Exception.
    // The code below would return success on empty string, since a path is prepended and thus the location exists
    if (filename_new.trim().empty()) throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);

    //add data dir in OpenMS data path
    directories.push_back(getOpenMSDataPath());

    //add path suffix to all specified directories
    String path = File::path(filename);
    if (path != "")
    {
      for (StringList::iterator it = directories.begin(); it != directories.end(); ++it)
      {
        it->ensureLastChar('/');
        *it += path;
      }
      filename_new = File::basename(filename);
    }

    //look up file
    for (StringList::const_iterator it = directories.begin(); it != directories.end(); ++it)
    {
      String loc = *it;
      loc.ensureLastChar('/');
      loc = loc + filename_new;

      if (exists(loc))
      {
        return String(QDir::cleanPath(loc.toQString()));
      }
    }

    //if the file was not found, throw an exception
    throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }

  bool File::fileList(const String& dir, const String& file_pattern, StringList& output, bool full_path)
  {
    QDir d(dir.toQString(), file_pattern.toQString(), QDir::Name, QDir::Files);
    QFileInfoList list = d.entryInfoList();

    //clear and check if empty
    output.clear();
    if (list.empty())
    {
      return false;
    }

    //resize output
    output.resize(list.size());

    //fill output
    UInt i = 0;
    for (QFileInfoList::const_iterator it = list.constBegin(); it != list.constEnd(); ++it)
    {
      output[i++] = full_path ? it->filePath() : it->fileName();
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
    static int number = 0;
    return now.getDate().remove('-') + "_" + now.getTime().remove(':') + "_" + (include_hostname ? String(QHostInfo::localHostName()) + "_" : "")  + pid + "_" + (++number);
  }

  String File::getOpenMSDataPath()
  {
    static String path;
    static bool path_checked = false;

    // we already checked the path, just return it
    // we do not support moving the path while OpenMS is running
    if (path_checked) return path;

    String found_path_from;
    bool from_env(false);
    if (getenv("OPENMS_DATA_PATH") != nullptr)
    {
      path = getenv("OPENMS_DATA_PATH");
      from_env = true;
      path_checked = isOpenMSDataPath_(path);
      if (path_checked) found_path_from = "OPENMS_DATA_PATH (environment)";
    }

    // probe the install path
    if (!path_checked)
    {
      path = OPENMS_INSTALL_DATA_PATH;
      path_checked = isOpenMSDataPath_(path);
      if (path_checked) found_path_from = "OPENMS_INSTALL_DATA_PATH (compiled)";
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
      if (path_checked) found_path_from = "tool path (run time)";
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
      String share_dir = "c:\\Program Files\\OpenMS\\share\\OpenMS";
#else
      String share_dir = "/usr/share/OpenMS";
#endif
      std::cerr << "  To resolve this, set the environment variable 'OPENMS_DATA_PATH' to the OpenMS share directory (e.g., '" + share_dir + "').\n";
      std::cerr << "Exiting now.\n";
      exit(1);
    }

    return path;
  }

  bool File::isOpenMSDataPath_(const String& path)
  {
    bool found = exists(path + "/CHEMISTRY/Elements.xml");
    return found;
  }

  bool File::isDirectory(const String& path)
  {
    QFileInfo fi(path.toQString());
    return fi.isDir();
  }

  String File::getTempDirectory()
  {
    Param p = getSystemParameters();
    String dir;
    if (getenv("OPENMS_TMPDIR") != 0)
    {
      dir = getenv("OPENMS_TMPDIR");
    }
    else if (p.exists("temp_dir") && String(p.getValue("temp_dir")).trim() != "")
    {
      dir = p.getValue("temp_dir");
    }
    else
    {
      dir = String(QDir::tempPath());
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
    else if (p.exists("home_dir") && String(p.getValue("home_dir")).trim() != "")
    {
      dir = p.getValue("home_dir");
    }
    else
    {
      dir = String(QDir::homePath());
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
      full_db_name = find(db_name, sys_p.getValue("id_db_dir"));
      OPENMS_LOG_INFO << "Augmenting database name '" << db_name << "' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '" << full_db_name << "'" << std::endl;
    }
    catch (Exception::FileNotFound& e)
    {
      OPENMS_LOG_ERROR << "Input database '" + db_name + "' not found (" << e.getMessage() << "). Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!" << std::endl;
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
      home_path = String(QDir::homePath());
    }
    return home_path;
  }

  Param File::getSystemParameters()
  {
    String home_path = File::getOpenMSHomePath();

    String filename = home_path + "/.OpenMS/OpenMS.ini";

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
    p.setValue("id_db_dir", ListUtils::create<String>(""),
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
    if (exists(exe_filename) && !isDirectory(exe_filename)) return true;

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

    if (File::exists(exec)) return exec;

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

  const String& File::getTemporaryFile(const String& alternative_file)
  {
    // take no action
    if (!alternative_file.empty()) return alternative_file;

    // create temporary (and schedule for deletion)
    return temporary_files_.newFile();
  }


  File::TemporaryFiles_::TemporaryFiles_()
    : filenames_()
  {
  }

  const String& File::TemporaryFiles_::newFile()
  {
    String s = getTempDirectory().ensureLastChar('/') + getUniqueName();
    filenames_.push_back(s);
    return filenames_.back();
  }

  File::TemporaryFiles_::~TemporaryFiles_()
  {
    for (Size i = 0; i < filenames_.size(); ++i)
    {
      if (File::exists(filenames_[i]) && !File::remove(filenames_[i])) 
      {
        std::cerr << "Warning: unable to remove temporary file '" << filenames_[i] << "'" << std::endl;
      }
    }
  }

  bool File::validateMatchingFileNames(const StringList& sl1, const StringList& sl2, bool basename, bool ignore_extension, bool strict)
  {
    // same number of filenames?
    if (sl1.size() != sl2.size()) return false;

    set<String> sl1_set;
    set<String> sl2_set;
    bool different_name_at_index = false;
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

      if (sl1_name != sl2_name) different_name_at_index = true;      
    }

    // Check for common mistake that order of input files have been switched.
    // This is the case if names (or basenames) are identical but the order does not match.
    bool same_set = (sl1_set == sl2_set);
    if (same_set && different_name_at_index) return false;

    // If we enforce a strict check then the sets of filenames must be identical.
    // Note that this can lead to problems if a workflow engine assigns random names to intermediate results.
    if (strict && !same_set) return false;

    return true;
  }

  File::TemporaryFiles_ File::temporary_files_;

} // namespace OpenMS
