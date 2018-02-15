// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#include <QDir>
#include <QFile>
#include <QCoreApplication>
#include <QtCore/QTimer>
#include <QtCore/QDateTime>

#include <OpenMS/DATASTRUCTURES/String.h>

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
#else
    platform = "Linux";
#endif

    // write to tmp + userid folder

    // e.g.: OpenMS_Default_Win_64_FeatureFinderCentroided_2.0.0
    String tool_version_string;
    tool_version_string = String("OpenMS") + "_" + "Default_" + platform + "_" + architecture + "_" + tool_name + "_" + version;

    String version_file_name = File::getOpenMSHomePath() + "/.OpenMS/" + tool_name + ".ver";

    // create version file if it doesn't exist yet
    bool first_run(false);
    if (!File::exists(version_file_name) || !File::readable(version_file_name))
    {
      // create OpenMS folder for .ver files
      String home_path = File::getOpenMSHomePath();
      String dirname = home_path + "/.OpenMS";
      QDir dir(dirname.toQString());

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
          LOG_INFO << "The OpenMS team is collecting usage statistics for quality control and funding purposes." << endl;
          LOG_INFO << "We will never give out your personal data, but you may disable this functionality by " << endl;
          LOG_INFO << "setting the environmental variable OPENMS_DISABLE_UPDATE_CHECK to ON." << endl;
        }
      
        // We need to use a QCoreApplication to fire up the  QEventLoop to process the signals and slots.
        char const * argv2[] = { "dummyname", nullptr };
        int argc = 1;
        QCoreApplication event_loop(argc, const_cast<char**>(argv2));
        NetworkGetRequest* query = new NetworkGetRequest(&event_loop);
        query->setUrl(QUrl(QString("http://openms-update.informatik.uni-tuebingen.de/check/") + tool_version_string.toQString()));
        QObject::connect(query, SIGNAL(done()), &event_loop, SLOT(quit()));
        QTimer::singleShot(1000, query, SLOT(run()));          
        QTimer::singleShot(5000, query, SLOT(timeOut()));
        event_loop.exec();

        if (!query->hasError())
        {
          if (debug_level > 0)
          {
            LOG_INFO << "Connecting to REST server successful. " << endl;
          }

          QString response = query->getResponse();
          VersionInfo::VersionDetails server_version = VersionInfo::VersionDetails::create(response);
          if (server_version != VersionInfo::VersionDetails::EMPTY)
          {
            if (VersionInfo::getVersionStruct() < server_version)
            {
              LOG_INFO << "Version " + version + " of " + tool_name + " is available at www.OpenMS.de" << endl;
            }
          }
        }
        else
        {
          if (debug_level > 0)
          {
            LOG_INFO << "Connecting to REST server failed. Skipping update check." << endl;
            LOG_INFO << "Error: " << String(query->getErrorString()) << endl;
          }
        }
        delete query;
      }
    }
  }

}

