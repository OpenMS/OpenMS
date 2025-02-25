// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/VISUAL/LayerDataBase.h>

class QLabel;
class QComboBox;
class QPushButton;
class QString;

#include <QtWidgets/QDialog>

namespace OpenMS
{
  class ParamEditor;
  class TVToolDiscovery;

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
      @param params Containing all TOPP tool/util params
      @param ini_file The file name of the temporary INI file created by this dialog
      @param default_dir The default directory for loading and storing INI files
      @param layer_type The type of data (determines the applicable tools)
      @param layer_name The name of the selected layer
      @param tool_scanner Pointer to the tool scanner for access to the plugins and to rerun the plugins detection
    */
    ToolsDialog(QWidget * parent, const Param& params, String ini_file, String default_dir, LayerDataBase::DataType layer_type, const String& layer_name, TVToolDiscovery* tool_scanner);
    ///Destructor
    ~ToolsDialog() override;

    /// to get the parameter name for output. Empty if no output was selected.
    String getOutput();
    /// to get the parameter name for input
    String getInput();
    /// to get the currently selected tool-name
    String getTool();
    /// get the default extension for the output file
    String getExtension();


private:
    /// ParamEditor for reading ini-files
    ParamEditor * editor_;
    /// tools description label
    QLabel * tool_desc_;
    /// ComboBox for choosing a TOPP-tool
    QComboBox * tools_combo_;
    /// Button to rerun the automatic plugin detection
    QPushButton* reload_plugins_button_;
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
    /// Location of the temporary INI file this dialog works on
    String ini_file_;
    /// default-dir of ini-file to open
    String default_dir_;
    /// name of ini-file
    QString filename_;
    /// Mapping of file extension to layer type to determine the type of a tool
    std::map<String, LayerDataBase::DataType> tool_map_;
    /// Param object containing all TOPP tool/util params
    Param tool_params_;
    /// Param object containing all plugin params
    Param plugin_params_;
    /// Pointer to the tool scanner for access to the plugins and to rerun the plugins detection
    TVToolDiscovery* tool_scanner_;
    /// The layer type of the current layer to determine all usable plugins
    LayerDataBase::DataType layer_type_;

    /// Disables the ok button and input/output comboboxes
    void disable_();
    /// Enables the ok button and input/output comboboxes
    void enable_();
    /// Determine all types a tool is compatible with by mapping each file extensions in a tools param
    std::vector<LayerDataBase::DataType> getTypesFromParam_(const Param& p) const;
    /// Fill input_combo_ and output_combo_ box with the appropriate entries from the specified param object.
    void setInputOutputCombo_(const Param& p);
    /// Create a list of all TOPP tool/util/plugins that are compatible with the active layer type
    QStringList createToolsList_();

protected slots:

    /// if ok button pressed show the tool output in a new layer, a new window or standard output as messagebox
    void ok_();
    /// Slot that handles changing of the tool
    void setTool_(int i);
    /// Slot that retrieves and displays the defaults
    void createINI_();
    /// loads an ini-file into the editor
    void loadINI_();
    /// stores an ini-file from the editor
    void storeINI_();
    /// rerun the automatic plugin detection
    void reloadPlugins_();
  };

}
