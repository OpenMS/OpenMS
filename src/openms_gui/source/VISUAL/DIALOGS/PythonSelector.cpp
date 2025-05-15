// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/DIALOGS/PythonSelector.h>
#include <ui_PythonSelector.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/PythonInfo.h>

#include <QString>
#include <QtWidgets/QFileDialog>
#include <QMessageBox>

using namespace std;

namespace OpenMS
{
  namespace Internal
  {
    PythonSelector::PythonSelector(QWidget* parent) :
    QWidget(parent),
    ui_(new Ui::PythonSelector)
    {
      ui_->setupUi(this);
      
      connect(ui_->btn_browse, SIGNAL(clicked()), this, SLOT(showFileDialog_()));
      connect(ui_->line_edit, SIGNAL(editingFinished()), this, SLOT(validate_()));

      // load/update UI
      ui_->line_edit->setText(last_known_python_exe_.toQString());

      // internally check
      validate_();
    }

    PythonSelector::~PythonSelector()
    {
      delete ui_;
      // TODO: store UI to INI?
    }

    
    void PythonSelector::showFileDialog_()
    {
      QString file_name = QFileDialog::getOpenFileName(this, tr("Specify Python executable"), tr(""), tr(/*valid formats*/ ""));
      if (!file_name.isEmpty())
      {
        ui_->line_edit->setText(file_name); // will not trigger the validator
        emit ui_->line_edit->editingFinished(); // simulate loosing focus or pressing return (to trigger validate_())
      }
    }

    void PythonSelector::validate_()
    {
      String exe = ui_->line_edit->text();
      
      String error;
      bool success = PythonInfo::canRun(exe, error);
      if (success)
      {
        last_known_python_exe_ = exe;
        ui_->label->setText(PythonInfo::getVersion(exe).toQString());
        currently_valid_ = true;
      }
      else
      {
        QMessageBox::warning(nullptr, QString("Python not found"), error.toQString());
        // no need to currently_valid_=false, since we will revert to 'last_known_python_exe_'
      }

      // reset to last known
      ui_->line_edit->setText(last_known_python_exe_.toQString());

      emit valueChanged(last_known_python_exe_.toQString(), currently_valid_);
    }


  }   //namespace Internal
} //namspace OpenMS
