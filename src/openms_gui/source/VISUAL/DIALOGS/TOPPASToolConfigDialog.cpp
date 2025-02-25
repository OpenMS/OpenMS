// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASToolConfigDialog.h>

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QtCore/QStringList>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QHBoxLayout>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QLabel>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QRadioButton>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCheckBox>
#include <QProcess>

using namespace std;

namespace OpenMS
{
  TOPPASToolConfigDialog::TOPPASToolConfigDialog(QWidget* parent, Param& param, const String& default_dir, const String& tool_name, const String& tool_type, const String& tool_desc, const QVector<String>& hidden_entries) :
    QDialog(parent),
    param_(&param),
    default_dir_(default_dir),
    tool_name_(tool_name),
    tool_type_(tool_type),
    hidden_entries_(hidden_entries)
  {
    QGridLayout* main_grid = new QGridLayout(this);

    QLabel* description = new QLabel;
    description->setAlignment(Qt::AlignTop | Qt::AlignLeft);
    description->setWordWrap(true);
    description->setText(tool_desc.toQString());
    main_grid->addWidget(description, 0, 0, 1, 1);

    //Add advanced mode check box
    editor_ = new ParamEditor(this);
    editor_->setMinimumSize(500, 500);
    main_grid->addWidget(editor_, 1, 0, 1, 1);

    QHBoxLayout* hbox = new QHBoxLayout;
    QPushButton* load_button = new QPushButton(tr("&Load config from .INI file"));
    connect(load_button, SIGNAL(clicked()), this, SLOT(loadINI_()));
    hbox->addWidget(load_button);
    QPushButton* store_button = new QPushButton(tr("&Store config to .INI file"));
    connect(store_button, SIGNAL(clicked()), this, SLOT(storeINI_()));
    hbox->addWidget(store_button);
    hbox->addStretch();

    // cancel button
    QPushButton* cancel_button = new QPushButton(tr("&Cancel"));
    connect(cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    hbox->addWidget(cancel_button);

    // ok button
    QPushButton* ok_button_ = new QPushButton(tr("&Ok"));
    connect(ok_button_, SIGNAL(clicked()), this, SLOT(ok_()));
    hbox->addWidget(ok_button_);

    main_grid->addLayout(hbox, 2, 0, 1, 1);

    setLayout(main_grid);

    editor_->load(*param_);
    editor_->setFocus(Qt::MouseFocusReason);

    setWindowTitle(tool_name.toQString() + " " + tr("configuration"));
  }

  TOPPASToolConfigDialog::~TOPPASToolConfigDialog() = default;

  void TOPPASToolConfigDialog::ok_()
  {
    if (editor_->isModified())
    {
      editor_->store();
      accept();
    }
    else
    {
      reject();
    }
  }

  void TOPPASToolConfigDialog::loadINI_()
  {
    QString string;
    filename_ = QFileDialog::getOpenFileName(this, tr("Open ini file"), default_dir_.c_str(), tr("ini files (*.ini);; all files (*.*)"));
    //no file selected
    if (filename_.isEmpty())
    {
      return;
    }
    if (!arg_param_.empty())
    {
      arg_param_.clear();
      param_->clear();
      editor_->clear();
    }
    try
    {
      ParamXMLFile paramFile;
      paramFile.load(filename_.toStdString(), arg_param_);
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this, "Error", (String("Error loading INI file: ") + e.what()).c_str());
      arg_param_.clear();
      return;
    }
    //Extract the required parameters
    *param_ = arg_param_.copy(tool_name_ + ":1:", true);
    //param_->remove("log");
    //param_->remove("no_progress");
    //param_->remove("debug");

    //remove parameters already explained by edges and the "type" parameter
    foreach(const String &name, hidden_entries_)
    {
      param_->remove(name);
    }

    //load data into editor
    editor_->load(*param_);
    editor_->setModified(true);
  }

  void TOPPASToolConfigDialog::storeINI_()
  {
    //nothing to save
    if (param_->empty())
      return;

    filename_ = QFileDialog::getSaveFileName(this, tr("Save ini file"), default_dir_.c_str(), tr("ini files (*.ini)"));
    //no file selected
    if (filename_.isEmpty())
      return;

    if (!filename_.endsWith(".ini"))
      filename_.append(".ini");

    bool was_modified = editor_->isModified();
    editor_->store();
    if (was_modified) editor_->setModified(true);

    arg_param_.insert(tool_name_ + ":1:", *param_);
    try
    {
      QString tmp_ini_file = File::getTempDirectory().toQString() + QDir::separator() + "TOPPAS_" + tool_name_.toQString() + "_";
      if (!tool_type_.empty())
      {
        tmp_ini_file += tool_type_.toQString() + "_";
      }
      tmp_ini_file += File::getUniqueName().toQString() + "_tmp.ini";
      //store current parameters
      ParamXMLFile paramFile;
      paramFile.store(tmp_ini_file.toStdString(), arg_param_);
      //restore other parameters that might be missing
      QString executable = File::findSiblingTOPPExecutable(tool_name_).toQString();
      QStringList args;
      args << "-write_ini" << filename_ << "-ini" << tmp_ini_file;
      if (!tool_type_.empty())
      {
        args << "-type" << tool_type_.toQString();
      }

      if (QProcess::execute(executable, args) != 0)
      {
        QMessageBox::critical(nullptr, "Error", (String("Could not execute '\"")  + executable + "\" \"" + args.join("\" \"") + "\"'!\n\nMake sure the TOPP tools are present in '" + File::getExecutablePath() + "', that you have permission to write to the temporary file path, and that there is space left in the temporary file path.").c_str());
        return;
      }
    }
    catch (Exception::BaseException& e)
    {
      QMessageBox::critical(this, "Error", (String("Error storing INI file: ") + e.what()).c_str());
      return;
    }
  }

}
