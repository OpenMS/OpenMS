// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/UpdateCheck.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#ifdef OPENMS_WINDOWSPLATFORM
#include <sys/utime.h>
#elif __APPLE__
#include <utime.h>
#else
#include <utime.h>
#endif

#include <sys/stat.h>

#include <OpenMS/SYSTEM/NetworkGetRequest.h>
#include <QtCore/QDir>
#include <QtCore/QCoreApplication>
#include <QtCore/QDateTime>
#include <QtCore/QTimer>

#include <OpenMS/CONCEPT/VersionInfo.h>

using namespace std;
  
namespace OpenMS
{
  void UpdateCheck::run(const String& tool_name, const String& version, int debug_level)
  {
    String architecture = QSysInfo::WordSize == 32 ? "32" : "64";

    // if the revision info is meaningful, show it as well
    String revision("UNKNOWN");
    if (!VersionInfo::getRevision().empty() && VersionInfo::getRevision() != "exported")
    {
      revision = VersionInfo::getRevision();
    }
    String platform;

#ifdef OPENMS_WINDOWSPLATFORM
    platform = "Win";
#elif __APPLE__
    platform = "Mac";
#elif __linux__
    platform = "Linux";
#elif __unix__
    platform = "Unix";
#else
    platform = "unknown";
#endif

    // write to tmp + userid folder

    // e.g.: OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0
    String tool_version_string;
    String config_path;
    //Comply with https://specifications.freedesktop.org/basedir-spec/basedir-spec-latest.html on unix identifying systems
    #ifdef __unix__
    if (getenv("XDG_CONFIG_HOME"))
    {
      config_path = String(getenv("XDG_CONFIG_HOME")) + "/OpenMS";
    }
    else
    {
      config_path = File::getOpenMSHomePath() + "/.config/OpenMS";
    }
    #else
    config_path =  File::getOpenMSHomePath() + "/.OpenMS";
    #endif
    tool_version_string = String("OpenMS") + "_" + "Default_" + platform + "_" + architecture + "_" + tool_name + "_" + version;

    String version_file_name = config_path + "/" + tool_name + ".ver";

    // create version file if it doesn't exist yet
    bool first_run(false);
    if (!File::exists(version_file_name) || !File::readable(version_file_name))
    {
      // create OpenMS folder for .ver files
      QDir dir(config_path.toQString());

      if (!dir.exists())
      {
        dir.mkpath(".");
      }

      // touch file to create it and set initial modification time stamp
      QFile f;
      f.setFileName(version_file_name.toQString());
      f.open(QIODevice::WriteOnly);
      f.close();
      first_run = true;
    }

    if (File::readable(version_file_name))
    {
      QDateTime last_modified_dt = QFileInfo(version_file_name.toQString()).lastModified();
      QDateTime current_dt = QDateTime::currentDateTime();

      // check if at least one day passed since last request
      if (first_run || current_dt > last_modified_dt.addDays(1))
      {
        // update modification time stamp
        struct stat old_stat;
        struct utimbuf new_times;
        stat(version_file_name.c_str(), &old_stat);
        new_times.actime = old_stat.st_atime; // keep accession time unchanged 
        new_times.modtime = time(nullptr);  // mod time to current time
        utime(version_file_name.c_str(), &new_times);          

        if (debug_level > 0)
        {
          OPENMS_LOG_INFO << "The OpenMS team is collecting usage statistics for quality control and funding purposes." << endl;
          OPENMS_LOG_INFO << "We will never give out your personal data, but you may disable this functionality by " << endl;
          OPENMS_LOG_INFO << "setting the environmental variable OPENMS_DISABLE_UPDATE_CHECK to ON." << endl;
        }
      
        // We need to use a QCoreApplication to fire up the  QEventLoop to process the signals and slots.
        char const * argv2[] = { "dummyname", nullptr };
        int argc = 1;
        QCoreApplication event_loop(argc, const_cast<char**>(argv2));
        NetworkGetRequest* query = new NetworkGetRequest(&event_loop);
        query->setUrl(QUrl(QString("http://openms-update.cs.uni-tuebingen.de/check/") + tool_version_string.toQString()));
        QObject::connect(query, SIGNAL(done()), &event_loop, SLOT(quit()));
        QTimer::singleShot(1000, query, SLOT(run()));          
        QTimer::singleShot(5000, query, SLOT(timeOut()));
        event_loop.exec();

        if (!query->hasError())
        {
          if (debug_level > 0)
          {
            OPENMS_LOG_INFO << "Connecting to REST server successful. " << endl;
          }

          QString response = query->getResponse();
          VersionInfo::VersionDetails server_version = VersionInfo::VersionDetails::create(response);
          if (server_version != VersionInfo::VersionDetails::EMPTY)
          {
            if (VersionInfo::getVersionStruct() < server_version)
            {
              OPENMS_LOG_INFO << "Version " + version + " of " + tool_name + " is available at www.OpenMS.de" << endl;
            }
          }
        }
        else
        {
          if (debug_level > 0)
          {
            OPENMS_LOG_INFO << "Connecting to REST server failed. Skipping update check." << endl;
            OPENMS_LOG_INFO << "Error: " << String(query->getErrorString()) << endl;
          }
        }
        delete query;
        event_loop.quit();
      }
    }
  }

}

