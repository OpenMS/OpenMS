// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/TOPPASVertexNameDialog.h>
#include <ui_TOPPASVertexNameDialog.h>

#include <QRegExpValidator> 

#include <iostream>

namespace OpenMS
{
  TOPPASVertexNameDialog::TOPPASVertexNameDialog(const QString& name, const QString& input_regex)
    : ui_(new Ui::TOPPASVertexNameDialogTemplate)
  {
    ui_->setupUi(this);
    
    if (!input_regex.isEmpty())
    {
      QRegExp rx(input_regex);
      QRegExpValidator* v = new QRegExpValidator(rx, ui_->line_edit);
      ui_->line_edit->setValidator(v);
    }

    ui_->line_edit->setText(name);
    connect(ui_->ok_button, SIGNAL(clicked()), this, SLOT(accept()));
    connect(ui_->cancel_button, SIGNAL(clicked()), this, SLOT(reject()));
  }

  TOPPASVertexNameDialog::~TOPPASVertexNameDialog()
  {
    delete ui_;
  }

  QString TOPPASVertexNameDialog::getName()
  {
    
    return ui_->line_edit->text();
  }

} // namespace
