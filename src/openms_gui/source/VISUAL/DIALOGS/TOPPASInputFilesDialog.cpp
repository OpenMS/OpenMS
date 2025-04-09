// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <ui_TOPPASInputFilesDialog.h>

#include <OpenMS/VISUAL/InputFileList.h>

namespace OpenMS
{
  TOPPASInputFilesDialog::TOPPASInputFilesDialog(const QStringList& list, const QString& cwd, QWidget* parent)
    : QDialog(parent),
      ui_(new Ui::TOPPASInputFilesDialogTemplate)
  {
    ui_->setupUi(this);
    ifl_ = (InputFileList*)ui_->input_file_list;
    ifl_->setCWD(cwd);
    ifl_->setFilenames(list);

    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
    setAcceptDrops(true);
  }

  TOPPASInputFilesDialog::~TOPPASInputFilesDialog()
  {
    delete ui_;
  }

  void TOPPASInputFilesDialog::getFilenames(QStringList& files) const
  {
    ifl_->getFilenames(files);
    if (ui_->flag_sort_list->isChecked())
      files.sort();
  }

  const QString& TOPPASInputFilesDialog::getCWD() const
  {
    return ifl_->getCWD();
  }



} // namespace
