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

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MISC/ExternalProcessMBox.h>

#include <QTabWidget>

#include <vector>
#include <utility> // for std::pair

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
    class SwathTabWidget;

    /**
      @brief RAII class to switch to certain TabWidget, disable the GUI and go back to the orignal Tab when this class is destroyed
    */
    class GUILock
    {
      public:
      GUILock(SwathTabWidget* stw);

      ~GUILock();
      private:
        SwathTabWidget* stw_;
        QWidget* old_;
        bool was_enabled_;
    };

    /// A multi-tabbed widget for the SwathWizard offering setting of parameters, input-file specification and running Swath and more
    class OPENMS_GUI_DLLAPI SwathTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      friend class GUILock;

      explicit SwathTabWidget(QWidget *parent = nullptr);
      ~SwathTabWidget();

      StringList getMzMLInputFiles() const;
    
    private slots:
      void on_run_swath_clicked();
      void on_edit_advanced_parameters_clicked();
      /// update the current working directory for all file input fields
      void broadcastNewCWD_(const QString& new_cwd);


      void on_btn_runPyProphet_clicked();

      void on_btn_pyresults_clicked();

    private:
      /// find the path of a Script, given the location of python(.exe). E.g. pyprophet.exe or feature_alignment.py
      /// Returns true on success, with the full path in @p script_name
      bool findPythonScript_(const String& path_to_python_exe, String& script_name);

      /// collect all parameters throughout the Wizard's controls and update 'swath_param_'
      void updateSwathParamFromWidgets_();

      /// update Widgets given a param object
      void updateWidgetsfromSwathParam_();

      /// where to write OSW output and pyProphet output
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
      void writeLog_(const QString& text, const QColor& color = "#000000", bool new_section = false);
      /// @brief convenient overload for String
      void writeLog_(const String& text, const QColor& color = "#000000", bool new_section = false);

      /// Ensure all input widgets are filled with data by the user to run OpenSwathWorkflow
      /// If anything is missing: show a Messagebox and return false.
      bool checkOSWInputReady_();

      Ui::SwathTabWidget *ui;
      Param swath_param_; ///< the global Swath parameters which will be passed to OpenSwathWorkflow.exe, once updated with parameters the Wizard holds separately
      Param swath_param_wizard_; ///< small selection of important parameters which the user can directly change in the Wizard

      StringList osw_result_files_; ///< list of .osw files produced by OSW which are currently available
      ExternalProcessMBox ep_; ///< to run external programs and pipe their output into our log
    };

  }
} // ns OpenMS

// this is required to allow Ui_SwathTabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
