// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/RWrapper.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/StringListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QProcess>

namespace OpenMS
{

  bool RWrapper::runScript( const String& script_file, const QStringList& cmd_args, const QString& executable /*= "Rscript"*/, bool find_R /*= false */, bool verbose /*= true */)
  {
    if (find_R && !findR(executable, verbose))
    {
      return false;
    }
    
    String fullscript;
    try
    {
      fullscript = findScript(script_file, verbose);
    }
    catch (...)
    {
      return false;
    }

    if (verbose)
    {
      OPENMS_LOG_INFO << "Running R script '" << fullscript << "' ...";
    }
    QStringList args;
    args << "--vanilla" << "--quiet" << fullscript.toQString();
    args.append(cmd_args);

    QProcess p;
    p.start(executable, args);
    p.waitForFinished(-1);

    if (p.error() == QProcess::FailedToStart || p.exitStatus() == QProcess::CrashExit || p.exitCode() != 0)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "\n--- ERROR MESSAGES ---\n";
        OPENMS_LOG_ERROR << QString(p.readAllStandardError()).toStdString();
        OPENMS_LOG_ERROR << "\n--- OTHER MESSAGES ---\n";
        OPENMS_LOG_ERROR << QString(p.readAllStandardOutput()).toStdString();
        OPENMS_LOG_ERROR << "\n\nScript failed. See above for an error description. " << std::endl;
      }
      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  bool RWrapper::findR( const QString& executable /*= "Rscript"*/, bool verbose /*= true*/ )
  {
    if (verbose) OPENMS_LOG_INFO << "Finding R interpreter 'Rscript' ...";

    QStringList args(QStringList() << "--vanilla" << "-e" << "sessionInfo()");
    QProcess p;
    p.setProcessChannelMode(QProcess::MergedChannels); // stdout receives all messages (stderr is empty)
    p.start(executable, args);
    p.waitForFinished(-1);

    if (p.error() == QProcess::FailedToStart)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        String out = QString(p.readAllStandardOutput()).toStdString();
        OPENMS_LOG_ERROR << "Error: Could not find or run '" << executable.toStdString() << "' executable (FailedToStart).\n";
        if (!out.empty())
        {
          OPENMS_LOG_ERROR << "Output was:\n------>\n"
                    << out
                    << "\n<------\n";
        }
        OPENMS_LOG_ERROR << "Please install 'Rscript', make sure it's in PATH and is flagged as executable." << std::endl;
      }

      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << "Trying to invoke 'Rscript' ...";
    }
    if (p.exitStatus() != QProcess::NormalExit || p.exitCode() != 0)
    {
      if (verbose)
      {
        OPENMS_LOG_INFO << " failed" << std::endl;
        OPENMS_LOG_ERROR << "Error: 'Rscript' executable returned with error (command: 'Rscript " << args.join(" ").toStdString() << "')\n"
                  << "Output was:\n------>\n"
                  << QString(p.readAllStandardOutput()).toStdString()
                  << "\n<------\n"
                  << "Make sure 'Rscript' is installed properly." << std::endl;
      }
      return false;
    }
    if (verbose)
    {
      OPENMS_LOG_INFO << " success" << std::endl;
    }
    return true;
  }

  OpenMS::String RWrapper::findScript( const String& script_file, bool verbose /*= true*/ )
  {
    String s;
    try
    {
      s = File::find(script_file, StringList(1, File::getOpenMSDataPath().ensureLastChar('/') + "SCRIPTS"));
    }
    catch (...)
    {
      if (verbose)
      {
        OPENMS_LOG_ERROR << "\n\nCould not find R script '" << script_file << "'!\n" << std::endl;
      }
      throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, script_file);
    }
    return s;
  }

} // namespace OpenMS
