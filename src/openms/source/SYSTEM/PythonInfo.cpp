// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/PythonInfo.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>
#include <QtCore/QDir>

#include <sstream>

using namespace std;

namespace OpenMS
{
  bool PythonInfo::canRun(String& python_executable, String& error_msg)
  {
    stringstream ss;
    String py_original = python_executable;
    if (!File::findExecutable(python_executable))
    {
      ss << "  Python not found at '" << python_executable << "'!\n"
         << "  Make sure Python is installed and this location is correct.\n";
      if (QDir::isRelativePath(python_executable.toQString()))
      {
        static String path;
        if (path.empty())
        {
          path = getenv("PATH");
        }
        ss << "  You might need to add the Python binary to your PATH variable\n"
           << "  or use an absolute path+filename pointing to Python.\n"
           << "  The current SYSTEM PATH is: '" << path << "'.\n\n";
#ifdef __APPLE__
        ss << "  On MacOSX, application bundles change the system PATH; Open your executable (e.g. KNIME/TOPPAS/TOPPView) from within the bundle (e.g. ./TOPPAS.app/Contents/MacOS/TOPPAS) to preserve the system PATH or use an absolute path to Python!\n";
#endif
      }
      error_msg = ss.str();
      return false;
    }

    if (python_executable != py_original)
    {
      ss << "Python executable ('" << py_original << "') resolved to '" << python_executable << "'\n";
    }

    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "--version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (!success)
    {
      if (qp.error() == QProcess::Timedout)
      {
        ss << "  Python was found at '" << python_executable << "' but the process timed out (can happen on very busy systems).\n"
           << "  Please free some resources or if you want to run the TOPP tool nevertheless set the TOPP tools 'force' flag in order to avoid this check.\n";
      }
      else if (qp.error() == QProcess::FailedToStart)
      {
        ss << "  Python found at '" << python_executable << "' but failed to run!\n"
           << "  Make sure you have the rights to execute this binary file.\n";
      }
      else
      {
        ss << "  Error executing '" << python_executable << "'!\n"
           << "  Error description: '" << qp.errorString().toStdString() << "'.\n";
      }
    }

    error_msg = ss.str();
    return success;
  }

  bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name)
  {
    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "-c" << (String("import ") + package_name).c_str(), QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    return (success && qp.exitStatus() == QProcess::ExitStatus::NormalExit && qp.exitCode() == 0);
  }

  String PythonInfo::getVersion(const String& python_executable)
  {
    String v;
    QProcess qp;
    qp.start(python_executable.toQString(), QStringList() << "--version", QIODevice::ReadOnly);
    bool success = qp.waitForFinished();
    if (success && qp.exitStatus() == QProcess::ExitStatus::NormalExit && qp.exitCode() == 0)
    {
      v = qp.readAllStandardOutput().toStdString(); // some pythons report is on stdout
      v += qp.readAllStandardError().toStdString();  // ... some on stderr
      v.trim(); // remove '\n'
    }
    return v;
  }

} // namespace OpenMS
