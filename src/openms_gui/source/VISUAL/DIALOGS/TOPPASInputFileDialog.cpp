// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFileDialog.h>
#include <ui_TOPPASInputFileDialog.h>

#include <QtWidgets/QMessageBox>
#include <QtCore/QFileInfo>

namespace OpenMS
{
  TOPPASInputFileDialog::TOPPASInputFileDialog(const QString& file_name)
    : QDialog(),
      ui_(new Ui::TOPPASInputFileDialogTemplate)
  {
    ui_->setupUi(this);

    ui_->input_file->setFilename(file_name);
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(checkValidity_()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
  }

  TOPPASInputFileDialog::~TOPPASInputFileDialog()
  {
    delete ui_;
  }

  void TOPPASInputFileDialog::setFileFormatFilter(const QString& fff)
  {
    ui_->input_file->setFileFormatFilter(fff);
  }

  QString TOPPASInputFileDialog::getFilename() const
  {
    return ui_->input_file->getFilename();
  }

  void TOPPASInputFileDialog::checkValidity_()
  {
    QFileInfo fi(getFilename());
    if (!(fi.exists() && fi.isReadable() && (!fi.isDir())))
    {
      QMessageBox::warning(nullptr, "Invalid file name", "Filename does not exist!");
      return; // do not close the dialog
    }

    accept();
  }


} // namespace
