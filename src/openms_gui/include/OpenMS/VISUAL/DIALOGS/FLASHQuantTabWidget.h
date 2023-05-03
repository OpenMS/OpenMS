// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Jihyung Kim $
// $Authors: Jihyung Kim $
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

    /**
      @brief RAII class to switch to certain TabWidget, disable the GUI and go back to the original Tab when this class is destroyed
    */
    class FLASHQuantGUILock
    {
      public:
      FLASHQuantGUILock(FLASHQuantTabWidget* ftw);

      ~FLASHQuantGUILock();
      private:
        FLASHQuantTabWidget* ftw_;
        QWidget* old_;
        GUIHelpers::GUILock glock_;
    };

    /// A multi-tabbed widget for the FLASHQuantWizard offering setting of parameters, input-file specification and running FLASHQuant and more
    class OPENMS_GUI_DLLAPI FLASHQuantTabWidget : public QTabWidget
    {
      Q_OBJECT

    public:
      friend class FLASHQuantGUILock;

      explicit FLASHQuantTabWidget(QWidget *parent = nullptr);
      ~FLASHQuantTabWidget();

      StringList getMzMLInputFiles() const;

    private slots:
      void on_run_fdq_clicked();
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
      void setWidgetsfromFDDefaultParam_();

      /// where to write output
      QString getCurrentOutDir_() const;
      QString infileToFDQoutput(const String& infile, const String& extension) const;

      /// append text to the log tab
      /// @param text The text to write
      /// @param new_section Start a new block with a date and time
      void writeLog_(const QString& text, const QColor& color = "#000000", bool new_section = false);
      /// @brief convenient overload for String
      void writeLog_(const String& text, const QColor& color = "#000000", bool new_section = false);

      /// Ensure all input widgets are filled with data by the user to run FLASHQuant
      /// If anything is missing: show a Messagebox and return false.
      bool checkFDQInputReady_();

      /// Group replicate LC-MS files with prefix of tags
      std::vector<std::vector<String>> groupReplicateFiles_(String prefix);

      /// run TopDownConsensusFeatureGroup
      void runTopDownConsensusFeatureGroup_();

      Ui::FLASHQuantTabWidget*ui;
      Param flashquant_param_; ///< the global FLASHQuant parameters which will be passed to FLASHQuant.exe, once updated with parameters the Wizard holds separately
      bool featurexml_output_; ///< true if it is requested

      ExternalProcessMBox ep_; ///< to run external programs and pipe their output into our log
    };

  }
} // ns OpenMS

// this is required to allow Ui_FLASHQuantTabWidget (auto UIC'd from .ui) to have a InputFile member
using InputFile = OpenMS::InputFile;
using OutputDirectory = OpenMS::OutputDirectory;
using ParamEditor = OpenMS::ParamEditor;
