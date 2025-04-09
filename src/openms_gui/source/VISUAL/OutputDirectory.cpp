// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/OutputDirectory.h>
#include <ui_OutputDirectory.h>


#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QMessageBox>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QCompleter>
#include <QFileSystemModel>


namespace OpenMS
{
  OutputDirectory::OutputDirectory(QWidget* parent)
    : QWidget(parent),
      ui_(new Ui::OutputDirectoryTemplate)
  {
    ui_->setupUi(this);
    QCompleter* completer = new QCompleter(this);
    QFileSystemModel* dir_model = new QFileSystemModel(completer);
    dir_model->setFilter(QDir::AllDirs);
    completer->setModel(dir_model);
    ui_->line_edit->setCompleter(completer);

    connect(ui_->browse_button, &QPushButton::clicked, this, &OutputDirectory::showFileDialog);
    connect(ui_->line_edit, &QLineEdit::textChanged, this, &OutputDirectory::textEditChanged_);
  }

  OutputDirectory::~OutputDirectory()
  {
    delete ui_;
  }

  void OutputDirectory::setDirectory(const QString& dir)
  {
    ui_->line_edit->setText(dir);
    emit directoryChanged(dir);
  }

  QString OutputDirectory::getDirectory() const
  {
    return ui_->line_edit->text();
  }

  void OutputDirectory::showFileDialog()
  {
    QString dir = File::exists(File::path(getDirectory())) ? File::path(getDirectory()).toQString() : "";
    QString selected_dir = QFileDialog::getExistingDirectory(this, tr("Select output directory"), dir);
    if (!selected_dir.isEmpty())
    {
      setDirectory(selected_dir); // emits directoryChanged()
    }
  }
  
  void OutputDirectory::textEditChanged_(const QString& /*new_text*/)
  {
    emit directoryChanged(getDirectory());
  }
  

  bool OutputDirectory::dirNameValid() const
  {
    if (!QFileInfo(getDirectory()).isDir())
    {
      return false;
    }
    QString file_name = getDirectory();
    if (!file_name.endsWith(QDir::separator()))
    {
      file_name += QDir::separator();
    }
    file_name += "test_file";
    return File::writable(file_name);
  }


} // namespace
