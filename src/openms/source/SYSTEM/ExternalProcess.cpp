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
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/SYSTEM/ExternalProcess.h>

#include <OpenMS/DATASTRUCTURES/String.h>


#include <QProcess>
#include <QStringList>
#include <QCoreApplication>


namespace OpenMS
{

  /// default Ctor; callbacks for stdout/stderr are empty
  ExternalProcess::ExternalProcess()
    : ExternalProcess([&](const String& /*out*/) {}, [&](const String& /*out*/) {}) // call other Ctor to connect signals!
  {
  }

  ExternalProcess::ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
    : qp_(new QProcess),
    callbackStdOut_(callbackStdOut),
    callbackStdErr_(callbackStdErr)
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
    callbackStdOut_ = callbackStdOut;
    callbackStdErr_ = callbackStdErr;
  }


  ExternalProcess::RETURNSTATE ExternalProcess::run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose)
  {
    String error_msg;
    return run(exe, args, working_dir, verbose, error_msg);
  }

  ExternalProcess::RETURNSTATE ExternalProcess::run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg)
  {
    error_msg.clear();
    if (!working_dir.isEmpty())
    {
      qp_->setWorkingDirectory(working_dir);
    }

    if (verbose)  callbackStdOut_("Running: " + (QStringList() << exe << args).join(' ') + '\n');

    qp_->start(exe, args);
    if (!(qp_->waitForStarted()))
    {
      error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::FAILED_TO_START;
    }
    while (qp_->state() == QProcess::Running)
    {
      QCoreApplication::processEvents();
      if (qp_->waitForReadyRead(50)) // wait 50msecs. Small enough to have the GUI repaint when switching windows
      {
        processStdOut_();
        processStdErr_();
      }
    }

    if (qp_->exitStatus() != QProcess::NormalExit)
    {
      error_msg = "Process '" + exe + "' crashed hard (segfault-like). Please check the log.";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::CRASH;
    }
    else if (qp_->exitCode() != 0)
    {
      error_msg = "Process '" + exe + "' did not finish successfully (exit code: " + int(qp_->exitCode()) + "). Please check the log.";
      if (verbose) callbackStdErr_(error_msg + '\n');
      return RETURNSTATE::NONZERO_EXIT;
    }
    
    if (verbose) callbackStdOut_("Executed '" + String(exe) + "' successfully!\n");
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
} // ns OpenMS
