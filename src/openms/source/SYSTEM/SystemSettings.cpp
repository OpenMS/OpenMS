// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/SystemSettings.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <cstdlib>
#include <filesystem>
#include <string>
#include <vector>

namespace fs = std::filesystem;
using std::getenv;

namespace OpenMS
{
  std::string SystemSettings::getTempDirectory()
  {
    Param p = getSystemParameters();
    std::string dir;
    if (getenv("OPENMS_TMPDIR") != nullptr)
    {
      dir = getenv("OPENMS_TMPDIR");
    }
    else if (p.exists("temp_dir") && !StringUtils::trimmed(p.getValue("temp_dir").toString()).empty())
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
  std::string SystemSettings::getUserDirectory()
  {
    Param p = getSystemParameters();
    std::string dir;
    if (getenv("OPENMS_HOME_PATH") != nullptr)
    {
      dir = getenv("OPENMS_HOME_PATH");
    }
    else if (p.exists("home_dir") && !StringUtils::trimmed(p.getValue("home_dir").toString()).empty())
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
      dir = home ? std::string(home) : std::string(".");
      StringUtils::substitute(dir, '\\', '/');
    }
    StringUtils::ensureLastChar(dir, '/');
    return dir;
  }

  std::string SystemSettings::findDatabase(const std::string& db_name)
  {
    Param sys_p = getSystemParameters();
    std::string full_db_name;
    try
    {
      full_db_name = File::find(db_name, ListUtils::toStringList<std::string>(sys_p.getValue("id_db_dir")));
      OPENMS_LOG_INFO << "Augmenting database name '" << db_name << "' with path given in 'OpenMS.ini:id_db_dir'. Full name is now: '" << full_db_name << "'\n";
    }
    catch (Exception::FileNotFound& e)
    {
      OPENMS_LOG_ERROR << "Input database '" + db_name + "' not found (" << e.what() << "). Make sure it exists (and check 'OpenMS.ini:id_db_dir' if you used relative paths. Aborting!\n";
      throw;
    }

    return full_db_name;
  }

  std::string SystemSettings::getOpenMSHomePath()
  {
    std::string home_path;
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
      home_path = home ? std::string(home) : std::string(".");
      StringUtils::substitute(home_path, '\\', '/');
    }
    return home_path;
  }

  std::string SystemSettings::getOpenMSConfigDir()
  {
    // Comply with https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html on unix identifying systems.
    // This is the single source of truth for the per-user config dir (OpenMS.ini, update-check .ver files, ...).
    #ifdef __unix__
      if (getenv("XDG_CONFIG_HOME"))
      {
        return std::string(getenv("XDG_CONFIG_HOME")) + "/OpenMS";
      }
      return SystemSettings::getOpenMSHomePath() + "/.config/OpenMS";
    #else
      return SystemSettings::getOpenMSHomePath() + "/.OpenMS";
    #endif
  }

  Param SystemSettings::getSystemParameters()
  {
    std::string filename = SystemSettings::getOpenMSConfigDir() + "/OpenMS.ini";

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

  Param SystemSettings::getSystemParameterDefaults_()
  {
    Param p;
    p.setValue("version", VersionInfo::getVersion());
    p.setValue("home_dir", ""); // only active when user enters something in this value
    p.setValue("temp_dir", ""); // only active when user enters something in this value
    p.setValue("id_db_dir", std::vector<std::string>(),
               std::string("Default directory for FASTA and psq files used as databased for id engines. ") + \
               "This allows you to specify just the filename of the DB in the " + \
               "respective TOPP tool, and the database will be searched in the directories specified here " + \
               ""); // only active when user enters something in this value
    p.setValue("threads", 1);

    return p;
  }

  // --------------------------------------------------------------------------
  // Deprecated File forwarders (kept for one release; see SystemSettings)
  // --------------------------------------------------------------------------
  std::string File::getOpenMSHomePath() { return SystemSettings::getOpenMSHomePath(); }
  std::string File::getOpenMSConfigDir() { return SystemSettings::getOpenMSConfigDir(); }
  std::string File::getTempDirectory() { return SystemSettings::getTempDirectory(); }
  std::string File::getUserDirectory() { return SystemSettings::getUserDirectory(); }
  Param File::getSystemParameters() { return SystemSettings::getSystemParameters(); }
  std::string File::findDatabase(const std::string& db_name) { return SystemSettings::findDatabase(db_name); }

} // namespace OpenMS
