// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/MISC/ExternalProcessMBox.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>
#include <OpenMS/VISUAL/TableView.h>
#include <QTabWidget> // our base class
#include <utility>    // for std::pair
#include <vector>

namespace Ui
{
  class FLASHDeconvTabWidget;
}

namespace OpenMS
{
  class InputFile;
  class OutputDirectory;
  class ParamEditor;


  namespace Internal
  {
    class FLASHDeconvTabWidget;

    /// A multi-tabbed widget for the FLASHDeconvWizard offering setting of parameters, input-file specification and running FLASHDeconv and more
    class OPENMS_GUI_DLLAPI FLASHDeconvTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      template <typename> friend class WizardGUILock;

      /// constructor
      explicit FLASHDeconvTabWidget(QWidget* parent = nullptr);
      /// Destructor
      ~FLASHDeconvTabWidget();

      /// get all the input mzML files as a string list
      StringList getMzMLInputFiles() const;

    private slots:
      void on_run_fd_clicked();
      void on_edit_advanced_parameters_clicked();
      void on_open_output_directory_clicked();
      /// update the current working directory for all file input fields
      void broadcastNewCWD_(const QString& new_cwd);

    private:
      /// collect all parameters throughout the Wizard's controls and update 'flashdeconv_param_'
      void updateFLASHDeconvParamFromWidgets_();

      /// collect output format parameters from the Wizard's control and update 'flashdeconv_output_tags_'
      void updateOutputParamFromWidgets_();

      /// update 'flashdeconv_param_outputs' with given input file name
      void updateOutputParamFromPerInputFile(const QString& input_file_name);

      /// update Widgets given a param object
      void setWidgetsfromFDDefaultParam_();

      /// where to write output
      QString getCurrentOutDir_() const;

      /// append text to the log tab
      /// @param text The text to write
      /// @param color Color for the text
      /// @param new_section Start a new block with a date and time
      void writeLog_(const QString& text, const QColor& color = "#000000", bool new_section = false);
      /// @brief convenient overload for String
      void writeLog_(const String& text, const QColor& color = "#000000", bool new_section = false);

      /// Ensure all input widgets are filled with data by the user to run FLASHDeconv
      /// If anything is missing: show a Messagebox and return false.
      bool checkFDInputReady_();

      Ui::FLASHDeconvTabWidget* ui;
      Param flashdeconv_param_;            ///< the global FLASHDeconv parameters which will be passed to FLASHDeconv.exe, once updated with parameters the Wizard holds separately
      Param flashdeconv_param_outputs_;    ///< Parameter set for different output formats
      StringList flashdeconv_output_tags_; ///< list of output parameter names checked by the user

      ExternalProcessMBox ep_; ///< to run external programs and pipe their output into our log
    };

  } // namespace Internal
} // namespace OpenMS

// this is required to allow Ui_FLASHDeconvTabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
