// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
  class FLASHQuantTabWidget;
}

namespace OpenMS
{
  class InputFile;
  class OutputDirectory;
  class ParamEditor;


  namespace Internal
  {
    class FLASHQuantTabWidget;

    /// A multi-tabbed widget for the FLASHQuantWizard offering setting of parameters, input-file specification and running FLASHQuant and more
    class OPENMS_GUI_DLLAPI FLASHQuantTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      template <typename> friend class WizardGUILock;

      /// constructor
      explicit FLASHQuantTabWidget(QWidget *parent = nullptr);
      /// Destructor
      ~FLASHQuantTabWidget();

      /// get all the input mzML files as a string list
      StringList getMzMLInputFiles() const;

    private slots:
      void on_run_fq_clicked();
      void on_open_output_directory_clicked();
      void on_consensus_names_check_clicked();
      /// update the current working directory for all file input fields
      void broadcastNewCWD_(const QString& new_cwd);

    private:
      /// collect all parameters throughout the Wizard's controls and update 'FLASHQuant_param_'
      void updateFLASHQuantParamFromWidgets_();

      /// collect output format parameters from the Wizard's control and update 'FLASHQuant_output_tags_'
      void updateOutputParamFromWidgets_();

      /// update Widgets given a param object
      void setWidgetsfromFQDefaultParam_();

      /// where to write output
      QString getCurrentOutDir_() const;
      QString infileToFQoutput(const String& infile, const String& extension) const;

      /// append text to the log tab
      /// @param text The text to write
      /// @param new_section Start a new block with a date and time
      void writeLog_(const QString& text, const QColor& color = "#000000", bool new_section = false);
      /// @brief convenient overload for String
      void writeLog_(const String& text, const QColor& color = "#000000", bool new_section = false);

      /// Ensure all input widgets are filled with data by the user to run FLASHQuant
      /// If anything is missing: show a Messagebox and return false.
      bool checkFQInputReady_();

      /// Group replicate LC-MS files with prefix of tags
      std::vector<std::vector<String>> groupReplicateFiles_(String prefix);

      /// run TopDownConsensusFeatureGroup
      void runTopDownConsensusFeatureGroup_();

      Ui::FLASHQuantTabWidget*ui;
      Param flashquant_param_; ///< the global FLASHQuant parameters which will be passed to FLASHQuant.exe, once updated with parameters the Wizard holds separately
      bool featurexml_output_; ///< true if it is requested

      ExternalProcessMBox ep_; ///< to run external programs and pipe their output into our log
    };

  } // namespace Internal
} // namespace OpenMS

// this is required to allow Ui_FLASHQuantTabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
