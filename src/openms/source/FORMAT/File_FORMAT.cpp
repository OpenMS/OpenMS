// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow, Marc Sturm $
// --------------------------------------------------------------------------

// This file contains the definitions of those OpenMS::File member functions
// that require knowledge from the FORMAT layer (file-type recognition via
// FileHandler::stripExtension and INI-file loading via ParamXMLFile).
//
// File itself lives in the SYSTEM layer (a low-level utility), so SYSTEM/File.cpp
// must not include FORMAT headers (that would be an upward/back-edge in the
// include layering). These few convenience methods are therefore defined here,
// in the FORMAT layer, where including FileHandler.h / ParamXMLFile.h is a normal
// forward edge. The methods stay as File:: members to keep File's public API
// (and all existing callers) unchanged.

#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>

#include <OpenMS/DATASTRUCTURES/Param.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <set>

using namespace std;

namespace OpenMS
{

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

} // namespace OpenMS
