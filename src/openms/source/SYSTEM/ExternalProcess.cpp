// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <OpenMS/DATASTRUCTURES/String.h>


#include <QtCore/QCoreApplication>
#include <QtCore/QProcess>
#include <QtCore/QStringList>
#include <utility>


namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcess::ExternalProcess()
    : ExternalProcess([&](const String& /*out*/) {}, [&](const String& /*out*/) {}) // call other Ctor to connect signals!
  {
  }

  ExternalProcess::ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
    : qp_(new QProcess),
    callbackStdOut_(std::move(callbackStdOut)),
    callbackStdErr_(std::move(callbackStdErr))
  {
    connect(qp_, &QProcess::readyReadStandardOutput, this, &ExternalProcess::processStdOut_);
    connect(qp_, &QProcess::readyReadStandardError, this, &ExternalProcess::processStdErr_);
  }

  ExternalProcess::~ExternalProcess()
  {
    delete qp_;
  }

  /// re-wire the callbacks used using run()
  void ExternalProcess::setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
  {
    callbackStdOut_ = std::move(callbackStdOut);
    callbackStdErr_ = std::move(callbackStdErr);
  }


  ExternalProcess::RETURNSTATE ExternalProcess::run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, IO_MODE io_mode)
  {
    String error_msg;
    return run(exe, args, working_dir, verbose, error_msg, io_mode);
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg, IO_MODE io_mode)
  {
    error_msg.clear();
    if (!working_dir.isEmpty())
    {
      qp_->setWorkingDirectory(working_dir);
    }

    if (verbose)
    {
      callbackStdOut_("Running: " + (QStringList() << exe << args).join(' ') + '\n');
    }
    // Map IO_MODE enum value to QIODevice value
    QIODevice::OpenModeFlag mode;
    switch (io_mode)
    {
      case IO_MODE::NO_IO:
        mode = QIODevice::NotOpen;
        break;
      case IO_MODE::READ_ONLY:
        mode = QIODevice::ReadOnly;
        break;
      case IO_MODE::WRITE_ONLY:
        mode = QIODevice::WriteOnly;
        break;
      default:
        mode = QIODevice::ReadWrite;
    }

    qp_->start(exe, args, mode);
    if (!(qp_->waitForStarted()))
    {
      error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
      if (verbose)
      {
        callbackStdErr_(error_msg + '\n');
      }
      return RETURNSTATE::FAILED_TO_START;
    }
    while (qp_->state() == QProcess::Running)
    {
      QCoreApplication::processEvents();
      if (qp_->waitForReadyRead(50)) // wait 50ms. Small enough to have the GUI repaint when switching windows
      {
        processStdOut_();
        processStdErr_();
      }
    }

    if (qp_->exitStatus() != QProcess::NormalExit)
    {
      error_msg = "Process '" + exe + "' crashed hard (segfault-like). Please check the log.";
      if (verbose)
      {
        callbackStdErr_(error_msg + '\n');
      }
      return RETURNSTATE::CRASH;
    }
    else if (qp_->exitCode() != 0)
    {
      error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + String(int(qp_->exitCode())).toQString() + "). Please check the log.";
      if (verbose)
      {
        callbackStdErr_(error_msg + '\n');
      }
      return RETURNSTATE::NONZERO_EXIT;
    }
    
    if (verbose)
    {
      callbackStdOut_("Executed '" + String(exe) + "' successfully!\n");
    }
    return RETURNSTATE::SUCCESS;
  }

  void ExternalProcess::processStdOut_()
  {
    String s(QString(qp_->readAllStandardOutput()));
    //std::cout << s << "\n";
    callbackStdOut_(s);
  }

  void ExternalProcess::processStdErr_()
  {
    String s(QString(qp_->readAllStandardError()));
    //std::cout << s << "\n";
    callbackStdErr_(s);
  }

} // namespace OpenMS
