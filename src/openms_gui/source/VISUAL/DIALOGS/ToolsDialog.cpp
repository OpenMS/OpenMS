// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/ToolsDialog.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/APPLICATIONS/ToolHandler.h>
#include <OpenMS/DATASTRUCTURES/Param.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/VISUAL/TVToolDiscovery.h>
#include <OpenMS/VISUAL/MISC/CommonDefs.h>

#include <QProcess>
#include <QtCore/QStringList>
#include <QtWidgets/QCheckBox>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QRadioButton>
#include <utility>


using namespace std;

namespace OpenMS
{

  ToolsDialog::ToolsDialog(
          QWidget* parent,
          const Param& params,
          String ini_file,
          String default_dir,
          LayerDataBase::DataType layer_type,
          const String& layer_name,
          TVToolDiscovery* tool_scanner) :
            QDialog(parent),
            ini_file_(std::move(ini_file)),
            default_dir_(std::move(default_dir)),
            tool_params_(params.copy("tool_params:", true)),
            plugin_params_(),
            tool_scanner_(tool_scanner),
            layer_type_(layer_type)
  {
    auto main_grid = new QGridLayout(this);

    // Layer label
    auto layer_label = new QLabel("Selected Layer:");
    main_grid->addWidget(layer_label, 0, 0);
    auto layer_label_name = new QLabel(layer_name.toQString());
    main_grid->addWidget(layer_label_name, 0, 1);

    auto label = new QLabel("TOPP tool:");
    main_grid->addWidget(label, 1, 0);

    // Determine all available tools compatible with the layer_type
    tool_map_ = {
            {FileTypes::Type::MZML, LayerDataBase::DataType::DT_PEAK},
            {FileTypes::Type::MZXML, LayerDataBase::DataType::DT_PEAK},
            {FileTypes::Type::FEATUREXML, LayerDataBase::DataType::DT_FEATURE},
            {FileTypes::Type::CONSENSUSXML, LayerDataBase::DataType::DT_CONSENSUS},
            {FileTypes::Type::IDXML, LayerDataBase::DataType::DT_IDENT}
    };

    QStringList list = createToolsList_();

    tools_combo_ = new QComboBox;
    tools_combo_->setMinimumWidth(150);
    tools_combo_->addItems(list);
    connect(tools_combo_, CONNECTCAST(QComboBox, activated, (int)), this, &ToolsDialog::setTool_);

    main_grid->addWidget(tools_combo_, 1, 1);

    reload_plugins_button_ = new QPushButton("Reload Plugins");
    connect(reload_plugins_button_, &QPushButton::clicked, this, &ToolsDialog::reloadPlugins_);
    main_grid->addWidget(reload_plugins_button_, 0, 2);

    label = new QLabel("input argument:");
    main_grid->addWidget(label, 2, 0);
    input_combo_ = new QComboBox;
    main_grid->addWidget(input_combo_, 2, 1);

    label = new QLabel("output argument:");
    main_grid->addWidget(label, 3, 0);
    output_combo_ = new QComboBox;
    main_grid->addWidget(output_combo_, 3, 1);

    // tools description label
    tool_desc_ = new QLabel;
    tool_desc_->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    tool_desc_->setWordWrap(true);
    main_grid->addWidget(tool_desc_, 1, 2, 3, 1);

    //Add advanced mode check box
    editor_ = new ParamEditor(this);
    main_grid->addWidget(editor_, 4, 0, 1, 5);

    auto hbox = new QHBoxLayout;
    auto load_button = new QPushButton(tr("&Load"));
    connect(load_button, &QPushButton::clicked, this, &ToolsDialog::loadINI_);
    hbox->addWidget(load_button);
    auto store_button = new QPushButton(tr("&Store"));
    connect(store_button, &QPushButton::clicked, this, &ToolsDialog::storeINI_);
    hbox->addWidget(store_button);
    hbox->addStretch();

    ok_button_ = new QPushButton(tr("&Ok"));
    connect(ok_button_, &QPushButton::clicked, this, &ToolsDialog::ok_);
    hbox->addWidget(ok_button_);

    auto cancel_button = new QPushButton(tr("&Cancel"));
    connect(cancel_button, &QPushButton::clicked, this, &ToolsDialog::reject);
    hbox->addWidget(cancel_button);
    main_grid->addLayout(hbox, 5, 0, 1, 5);

    setLayout(main_grid);

    setWindowTitle(tr("Apply TOPP tool to layer"));
    disable_();
  }

  ToolsDialog::~ToolsDialog() = default;

  std::vector<LayerDataBase::DataType> ToolsDialog::getTypesFromParam_(const Param& p) const
  {
    // Containing all types a tool is compatible with
    std::vector<LayerDataBase::DataType> types;
    for (const auto& entry : p)
    {
      if (entry.name == "in")
      {
        // Map all file extension to a LayerDataBase::DataType
        for (auto& file_extension : entry.valid_strings)
        {
          // a file extension in valid_strings is of form "*.TYPE" -> convert to substr "TYPE".
          const String& file_type = file_extension.substr(2, file_extension.size());
          const auto& iter = tool_map_.find(FileTypes::nameToType(file_type));
          // If mapping was found
          if (iter != tool_map_.end())
          {
            types.push_back(iter->second);
          }
        }
        break;
      }
    }
    return types;
  }

  void ToolsDialog::setInputOutputCombo_(const Param &p)
  {
    String str;
    QStringList input_list("<select>");
    QStringList output_list("<select>");
    bool outRequired = false;
    for (Param::ParamIterator iter = p.begin(); iter != p.end(); ++iter)
    {
      // iter.getName() is either of form "ToolName:1:ItemName" or "ToolName:1:NodeName:[...]:ItemName".
      // Cut off "ToolName:1:"
      str = iter.getName().substr(iter.getName().rfind("1:") + 2, iter.getName().size());
      // Only add items and no nodes
      if (!str.empty() && str.find(":") == String::npos)
      {
        // Only add to input list if item has "input file" tag.
        if (iter->tags.find("input file") != iter->tags.end())
        {
          input_list << QStringList(str.c_str());
        }
          // Only add to output list if item has "output file" tag.
        else if (iter->tags.find("output file") != iter->tags.end())
        {
          output_list << QStringList(str.c_str());
          // Check whether the item has a required tag i.e. is mandatory.
          outRequired = (outRequired) || (iter->tags.find("required") != iter->tags.end());
        }
      }
    }
    // Clear and set input combo box
    input_combo_->clear();
    output_combo_->clear();
    input_combo_->addItems(input_list);
    Int pos = input_list.indexOf("in");
    if (pos != -1)
    {
      input_combo_->setCurrentIndex(pos);
    }
    // Clear and set output combo box
    output_combo_->addItems(output_list);
    pos = output_list.indexOf("out");
    if (pos != -1 && getTool() != "FileInfo" && outRequired)
    {
      output_combo_->setCurrentIndex(pos);
    }
  }

  QStringList ToolsDialog::createToolsList_()
  {
    //Make sure the list is empty
    QStringList list;

    const auto& tools = ToolHandler::getTOPPToolList();
    plugin_params_ = tool_scanner_->getPluginParams();

    for (auto& pair : tools)
    {
      std::vector<LayerDataBase::DataType> tool_types = getTypesFromParam_(tool_params_.copy(pair.first + ':'));
      if (std::find(tool_types.begin(), tool_types.end(), layer_type_) != tool_types.end())
      {
        list << pair.first.toQString();
      }
    }
    //TODO: Plugins get added to the list just like tools and can't be differentiated in the GUI
    for (const auto& name : tool_scanner_->getPlugins())
    {
      std::vector<LayerDataBase::DataType> tool_types = getTypesFromParam_(plugin_params_.copy(name + ":"));
      if (std::find(tool_types.begin(), tool_types.end(), layer_type_) != tool_types.end())
      {
        list << String(name).toQString();
      }
    }

    //sort list alphabetically
    list.sort();
    list.push_front("<select tool>");
    return list;
  }

  void ToolsDialog::createINI_()
  {
    enable_();
    if (!arg_param_.empty())
    {
       tool_desc_->clear();
       arg_param_.clear();
       vis_param_.clear();
       editor_->clear();
    }
    auto tool_name = getTool();
    arg_param_ = tool_params_.copy(tool_name + ":");
    if (arg_param_.empty())
    {
      arg_param_ = plugin_params_.copy(tool_name + ":");
    }

    tool_desc_->setText(String(arg_param_.getSectionDescription(tool_name)).toQString());
    vis_param_ = arg_param_.copy(tool_name + ":1:", true);
    vis_param_.remove("log");
    vis_param_.remove("no_progress");
    vis_param_.remove("debug");

    editor_->load(vis_param_);

    setInputOutputCombo_(arg_param_);

    editor_->setFocus(Qt::MouseFocusReason);
  }

  void ToolsDialog::setTool_(int i)
  {
    editor_->clear();

    // no tool selected
    if (i == 0)
    {
      disable_();
      return;
    }

    createINI_();
  }

  void ToolsDialog::disable_()
  {
    ok_button_->setEnabled(false);
    input_combo_->setCurrentIndex(0);
    input_combo_->setEnabled(false);
    output_combo_->setCurrentIndex(0);
    output_combo_->setEnabled(false);
  }

  void ToolsDialog::enable_()
  {
    ok_button_->setEnabled(true);
    input_combo_->setEnabled(true);
    output_combo_->setEnabled(true);
  }

  void ToolsDialog::ok_()
  {
    if (input_combo_->currentText() == "<select>" || tools_combo_->currentText() == "<select>")
    {
      QMessageBox::critical(this, "Error", "You have to select a tool and an input argument!");
    }
    else
    {
      editor_->store();
      arg_param_.insert(getTool() + ":1:", vis_param_);
      if (!File::writable(ini_file_))
      {
        QMessageBox::critical(this, "Error", (String("Could not write to '") + ini_file_ + "'!").c_str());
      }
      ParamXMLFile paramFile;
      paramFile.store(ini_file_, arg_param_);
      accept();
    }
  }

  void ToolsDialog::loadINI_()
  {
    QString string;
    filename_ = QFileDialog::getOpenFileName(this, tr("Open ini file"), default_dir_.c_str(), tr("ini files (*.ini);; all files (*.*)"));
    //not file selected
    if (filename_.isEmpty())
    {
      return;
    }
    enable_();
    if (!arg_param_.empty())
    {
      arg_param_.clear();
      vis_param_.clear();
      editor_->clear();
    }
    try
    {
      ParamXMLFile paramFile;
      paramFile.load(filename_.toStdString(), arg_param_);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this, "Error", QString("Error loading INI file: ") + e.what());
      arg_param_.clear();
      return;
    }
    //set tool combo
    Param::ParamIterator iter = arg_param_.begin();
    String str;
    string = iter.getName().substr(0, iter.getName().find(":")).c_str();
    Int pos = tools_combo_->findText(string);
    if (pos == -1)
    {
      QMessageBox::critical(this, "Error", (String("Cannot apply '") + string + "' tool to this layer type. Aborting!").c_str());
      arg_param_.clear();
      return;
    }
    tools_combo_->setCurrentIndex(pos);
    //Extract the required parameters
    vis_param_ = arg_param_.copy(getTool() + ":1:", true);
    vis_param_.remove("log");
    vis_param_.remove("no_progress");
    vis_param_.remove("debug");
    //load data into editor
    editor_->load(vis_param_);

    setInputOutputCombo_(arg_param_);
  }

  void ToolsDialog::storeINI_()
  {
    //nothing to save
    if (arg_param_.empty())
      return;

    filename_ = QFileDialog::getSaveFileName(this, tr("Save ini file"), default_dir_.c_str(), tr("ini files (*.ini)"));
    //not file selected
    if (filename_.isEmpty())
    {
      return;
    }

    if (!filename_.endsWith(".ini"))
    {
      filename_.append(".ini");
    }
    editor_->store();
    arg_param_.insert(getTool() + ":1:", vis_param_);
    try
    {
      ParamXMLFile paramFile;
      paramFile.store(filename_.toStdString(), arg_param_);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this, "Error", QString("Error storing INI file: ") + e.what());
      return;
    }
  }

  void ToolsDialog::reloadPlugins_()
  {
    QStringList list = createToolsList_();

    int32_t selected_index = list.indexOf(tools_combo_->currentText());

    if (selected_index < 1)
    {
      tool_desc_->clear();
      arg_param_.clear();
      vis_param_.clear();
      editor_->clear();
      input_combo_->clear();
      output_combo_->clear();
      disable_();
    }
    tools_combo_->clear();
    tools_combo_->addItems(list);
    if (selected_index > 0)
    {
      tools_combo_->setCurrentIndex(selected_index);
      createINI_();
    }
  }

  String ToolsDialog::getOutput()
  {
    if (output_combo_->currentText() == "<select>")
      return "";

    return output_combo_->currentText();
  }

  String ToolsDialog::getInput()
  {
    return input_combo_->currentText();
  }

  String ToolsDialog::getTool()
  {
    return tools_combo_->currentText();
  }

  String ToolsDialog::getExtension()
  {
    // Try to Return the first valid string for the extension on the output parameter
    // If we can't get any valid strings show an error.
    String extension = FileTypes::typeToName(FileTypes::UNKNOWN);
    auto validStrings = vis_param_.getValidStrings(output_combo_->currentText().toStdString()); 
    // If we have only one valid output type use that
    if (validStrings.size() == 1)
    {
      extension = validStrings[0];
      // Remove the leading .*
      extension = extension.suffix(extension.size() - 2);
    }
    // Otherwise the type is unknown
    else 
    {
      // If we have no valid types, produce an error
      if (validStrings.empty())
        {
          QMessageBox::critical(this, "Error", QString("Error determining output type from tool. Tool is not compatible with TOPPView"));
        }
          // If we have multiple valid output types, we don't know what the file actually contains, so use UNKNOWN

    }
    return extension;
  }

}
