// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

class QLabel;
class QComboBox;
class QPushButton;
class QRadioButton;
class QString;

#include <QtGui/QDialog>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/LayerData.h>

namespace OpenMS
{
  class ParamEditor;

  /**
      @brief TOPP tool selection dialog

      In the dialog, the user can
    - select a TOPP tool
    - select the options used for the input and output file
      - and set the parameters for the tool

      This information can then be used to execute the tool.

      The offered tools depend on the data type set in the constructor.

      @ingroup Dialogs
  */
  class OPENMS_GUI_DLLAPI ToolsDialog :
    public QDialog
  {
    Q_OBJECT

public:
    /**
        @brief Constructor

        @param parent Qt parent widget
        @param ini_file The file name of the temporary INI file created by this dialog
        @param default_dir The default directory for loading and storing INI files
        @param layertype The type of data (determines the applicable tools)
    */
    ToolsDialog(QWidget * parent, String ini_file, String default_dir, LayerData::DataType layertype);
    ///Destructor
    ~ToolsDialog() override;

    /// to get the parameter name for output. Empty if no output was selected.
    String getOutput();
    /// to get the parameter name for input
    String getInput();
    /// to get the currently selected tool-name
    String getTool();

private:
    /// ParamEditor for reading ini-files
    ParamEditor * editor_;
    /// tools description label
    QLabel * tool_desc_;
    /// ComboBox for choosing a TOPP-tool
    QComboBox * tools_combo_;
    /// for choosing an input parameter
    QComboBox * input_combo_;
    /// for choosing an output parameter
    QComboBox * output_combo_;
    /// Param for loading the ini-file
    Param arg_param_;
    /// Param for loading configuration information in the ParamEditor
    Param vis_param_;
    /// ok-button connected with slot ok_()
    QPushButton * ok_button_;
    /// map for getting the parameter name from the full path in arg_param
    std::map<String, String> arg_map_;
    /// Location of the temporary INI file this dialog works on
    String ini_file_;
    /// default-dir of ini-file to open
    String default_dir_;
    /// name of ini-file
    QString filename_;

    ///Disables the ok button and input/output comboboxes
    void disable_();
    ///Enables the ok button and input/output comboboxes
    void enable_();

protected slots:

    /// if ok button pressed show the tool output in a new layer, a new window or standard output as messagebox
    void ok_();
    /// Slot that handles changing of the tool
    void setTool_(int i);
    /// Slot that retrieves and displays the defaults
    void createINI_();
    /// loads an ini-file into the editor_
    void loadINI_();
    /// stores an ini-file from the editor_
    void storeINI_();
  };

}
#endif // OPENMS_VISUAL_DIALOGS_TOOLSDIALOG_H
