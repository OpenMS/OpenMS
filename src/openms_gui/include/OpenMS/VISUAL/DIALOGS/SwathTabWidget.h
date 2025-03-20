// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/VISUAL/MISC/ExternalProcessMBox.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/TableView.h>

#include <QTabWidget> // our base class

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

    /// A multi-tabbed widget for the SwathWizard offering setting of parameters, input-file specification and running Swath and more
    class OPENMS_GUI_DLLAPI SwathTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      template <typename> friend class WizardGUILock;

      explicit SwathTabWidget(QWidget *parent = nullptr);
      ~SwathTabWidget();

      StringList getMzMLInputFiles() const;

      QStringList getPyProphetOutputFileNames() const;

    private slots:
      void on_run_swath_clicked();
      void on_edit_advanced_parameters_clicked();
      /// update the current working directory for all file input fields
      void broadcastNewCWD_(const QString& new_cwd);


      void on_btn_runPyProphet_clicked();

      void on_btn_pyresults_clicked();

      void on_pushButton_clicked();

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
      /// @param color Color of the text
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
using TableView = OpenMS::TableView;
