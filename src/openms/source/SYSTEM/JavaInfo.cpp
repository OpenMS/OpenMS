// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/JavaInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QtCore/QDir>

namespace OpenMS
{

  bool JavaInfo::canRun(const String& java_executable, bool verbose_on_error)
  {
    QProcess qp;
    qp.start(java_executable.toQString(), QStringList() << "-version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (!success && verbose_on_error)
    {
        OPENMS_LOG_ERROR << "Java-Check:\n";
        if (qp.error() == QProcess::Timedout)
        {
          OPENMS_LOG_ERROR
            << "  Java was found at '" << java_executable << "' but the process timed out (can happen on very busy systems).\n"
            << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check." << std::endl;
        }
        else if (qp.error() == QProcess::FailedToStart)
        {
          OPENMS_LOG_ERROR
            << "  Java not found at '" << java_executable << "'!\n"
            << "  Make sure Java is installed and this location is correct.\n";
          if (QDir::isRelativePath(java_executable.toQString()))
          {
            static String path;
            if (path.empty())
            {
              path = getenv("PATH");
            }
            OPENMS_LOG_ERROR << "  You might need to add the Java binary to your PATH variable\n"
              << "  or use an absolute path+filename pointing to Java.\n" 
              << "  The current SYSTEM PATH is: '" << path << "'.\n\n"
  #ifdef __APPLE__
              << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Java!\n"
  #endif
              << std::endl;
          }
          else 
          {
            OPENMS_LOG_ERROR << "  You gave an absolute path to Java. Please check if it's correct.\n"
              << "  You can also try 'java' if your system path is correctly configured.\n"
              << std::endl;
          }
        }
        else
        {
          OPENMS_LOG_ERROR << "  Error executing '" << java_executable << "'!\n"
                    << "  Error description: '" << qp.errorString().toStdString() << "'.\n";
        }
    }
    return success;
  }

} // namespace OpenMS
