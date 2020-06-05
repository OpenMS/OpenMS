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

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/String.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <QTabWidget>

#include <vector>
#include <utility> // for std::pair
#include <functional> // for std::function

#include <QMessageBox>
#include <QProcess>

#include <QCoreApplication>

namespace Ui
{
  class SwathTabWidget;
}

namespace OpenMS
{
  class InputFile;
  class OutputDirectory;
  class ParamEditor;

  namespace Internal
  {

    /**
      @brief A wrapper around QProcess to conveniently start an external program and forward its outputs

      Use the custom Ctor to provide callback functions for stdout/stderr output or set them via setCallbacks().

      Running an external program blocks the caller, so do not use this in a main GUI thread
      (unless you have some other means to tell the user that no interaction is possible at the moment).

    */
    // TODO: this will replace TOPPBase::runExternalProgram_(), but in a separate PR
    class ExternalProcess
      : public QObject
    {
      Q_OBJECT

    public:
      enum class RETURNSTATE
      {
        SUCCESS,  ///< everything went smoothly (exit-code = 0)
        NONZERO_EXIT, /// ran, but returned with an exit-code other than 0
        CRASH, ///< ran, but crashed (segfault etc)
        FAILED_TO_START ///< executable not found or not enough access rights for user
      };

      /// default Ctor; callbacks for stdout/stderr are empty
      ExternalProcess()
      : ExternalProcess([&](const String& /*out*/) {}, [&](const String& /*out*/) {}) // call other Ctor to connect signals!
      {
      }

      ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
        : qp_(new QProcess),
          callbackStdOut_(callbackStdOut),
          callbackStdErr_(callbackStdErr)
      {
        connect(qp_, &QProcess::readyReadStandardOutput, this, &ExternalProcess::processStdOut_);
        connect(qp_, &QProcess::readyReadStandardError, this, &ExternalProcess::processStdErr_);
      }

      ~ExternalProcess()
      {
        delete qp_;
      }
      
      /// re-wire the callbacks used using run()
      void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr)
      {
        callbackStdOut_ = callbackStdOut;
        callbackStdErr_ = callbackStdErr;
      }

      /**
        @brief Runs a program and calls the callback functions from time to time if output from the external program is available.

        @param parent Optional parent widget, used to show QMesssageBoxes (see @p verbose_GUI)
        @param exe The program to call (can contain spaces in path, no problem)
        @param args A list of extra arguments (can be empty)
        @param verbose Report the call issued and errors to the callbacks you provided
        @param verbose_GUI Show QMessageBoxes with errors should they occur
      */
      RETURNSTATE run(QWidget* parent, const QString& exe, const QStringList& args, bool verbose, bool verbose_GUI)
      {
        String error_msg;
        if (verbose)  callbackStdOut_("Running: " + (QStringList() << exe << args).join(' ') + '\n');
        
        qp_->start(exe, args);
        if (!(qp_->waitForStarted()))
        {
          error_msg = "Process '" + exe + "' failed to start. Does it exist? Is it executable?";
          if (verbose) callbackStdErr_(error_msg + '\n');
          if (verbose_GUI) QMessageBox::critical(parent, "Error", error_msg.toQString());
          return RETURNSTATE::FAILED_TO_START;
        }
        while (qp_->state() == QProcess::Running)
        {
          QCoreApplication::processEvents();
          qp_->waitForReadyRead();
          processStdOut_();
          processStdErr_();
        }
        if (verbose) callbackStdOut_("Exit code: " + String(qp_->exitCode()));

        bool any_failure = qp_->exitStatus() != QProcess::NormalExit || qp_->exitCode() != 0;
        if (any_failure)
        {
          error_msg = "Process '" + exe + "' did not finish successfully. Please check the log.";
          if (verbose) callbackStdErr_(error_msg + '\n');
          if (verbose_GUI) QMessageBox::critical(parent, "Error", error_msg.toQString());
          if (qp_->exitCode() !=0) return RETURNSTATE::NONZERO_EXIT;
          else return RETURNSTATE::CRASH;
        }
        return RETURNSTATE::SUCCESS;
      }

    private slots:
      void processStdOut_()
      {
        String s(QString(qp_->readAllStandardOutput()));
        //std::cout << s << "\n";
        callbackStdOut_(s);
      }
      void processStdErr_()
      {
        String s(QString(qp_->readAllStandardError()));
        //std::cout << s << "\n";
        callbackStdErr_(s);
      }

    private:
      QProcess* qp_; ///< pointer to avoid including the QProcess header here (it's huge)
      std::function<void(const String&)> callbackStdOut_;
      std::function<void(const String&)> callbackStdErr_;
    };

    /// A multi-tabbed widget for the SwathWizard offering setting of parameters, input-file specification and running Swath and more
    class OPENMS_GUI_DLLAPI SwathTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      explicit SwathTabWidget(QWidget *parent = nullptr);
      ~SwathTabWidget();

      StringList getMzMLInputFiles() const;
    
    private slots:
      void on_run_swath_clicked();
      void on_edit_advanced_parameters_clicked();
      /// update the current working directory for all file input fields
      void broadcastNewCWD_(const QString& new_cwd);

      void on_btn_runPyProphet_clicked();

    private:
      /// collect all parameters throughout the Wizard's controls and update 'swath_param_'
      void updateSwathParamFromWidgets_();

      /// update Widgets given a param object
      void updateWidgetsfromSwathParam_();

      QString getCurrentOutDir_() const;

      /// translate the current list of input mzMLs and the current output directory of OSW to a list of expected OSW output files == pyProphet input files
      /// The bool indicates if the file is already present
      std::vector<std::pair<String, bool>> getPyProphetInputFiles() const;

      /// check if input to pyProphet is already present in the output directory of OSW
      void checkPyProphetInput_();

      /// fill osw_result_files_ according to the the currently specified input mzMLs


      /// append text to the log tab
      /// @param text The text to write
      /// @param new_section Start a new block with a date and time
      void writeLog_(const QString& text, bool new_section = false);

      /// Ensure all input widgets are filled with data by the user to run OpenSwathWorkflow
      /// If anything is missing: show a Messagebox and return false.
      bool checkOSWInputReady_();

      Ui::SwathTabWidget *ui;
      Param swath_param_; ///< the global Swath parameters which will be passed to OpenSwathWorkflow.exe, once updated with parameters the Wizard holds separately
      Param swath_param_wizard_; ///< small selection of important parameters which the user can directly change in the Wizard

      StringList osw_result_files_; ///< list of .osw files produced by OSW which are currently available
      ExternalProcess ep_; ///< to run external programs and pipe their output into our log
    };

  }
} // ns OpenMS

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
