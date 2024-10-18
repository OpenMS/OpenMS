// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <ui_TOPPASOutputFilesDialog.h>

#include <OpenMS/SYSTEM/File.h>

#include <QtWidgets/QMessageBox>
#include <QtCore/QDir>


namespace OpenMS
{
  TOPPASOutputFilesDialog::TOPPASOutputFilesDialog(const QString& dir_name, int num_jobs)
    : ui_(new Ui::TOPPASOutputFilesDialogTemplate)
  {
    ui_->setupUi(this);
    if (dir_name != "")
    {
      ui_->out_dir->setDirectory(dir_name);
    }
    else
    {
      ui_->out_dir->setDirectory(QDir::currentPath());
    }
    if (num_jobs >= 1)
    {
      ui_->num_jobs_box->setValue(num_jobs);
    }
    
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    
    // make Ok the default (just pressing Enter will run the workflow)
    ui_->ok_button->setFocus();
  }

  TOPPASOutputFilesDialog::~TOPPASOutputFilesDialog()
  {
    delete ui_;
  }

  void TOPPASOutputFilesDialog::showFileDialog()
  {
    ui_->out_dir->showFileDialog();
  }

  QString TOPPASOutputFilesDialog::getDirectory() const
  {
    return ui_->out_dir->getDirectory();
  }

  int TOPPASOutputFilesDialog::getNumJobs() const
  {
    return ui_->num_jobs_box->value();
  }

  void TOPPASOutputFilesDialog::checkValidity_()
  {
    if (!ui_->out_dir->dirNameValid())
    {
      QMessageBox::warning(nullptr, "Invalid directory", "Either the specified path is no directory, or you have no permission to write there.");
      return;
    }

    accept();
  }


} // namespace
