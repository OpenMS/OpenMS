// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Plot1DGoToDialog.h>
#include <ui_Plot1DGoToDialog.h>

#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  Plot1DGoToDialog::Plot1DGoToDialog(QWidget * parent) :
    QDialog(parent),
    ui_(new Ui::Plot1DGoToDialogTemplate)
  {
    ui_->setupUi(this);
  }

  Plot1DGoToDialog::~Plot1DGoToDialog()
  {
    delete ui_;
  }


  void Plot1DGoToDialog::setRange(float min, float max)
  {
    ui_->min_->setText(QString::number(min));
    ui_->max_->setText(QString::number(max));
  }

  void Plot1DGoToDialog::setMinMaxOfRange(float min, float max)
  {
    ui_->min_const_->setText(QString("min: ") + QString::number(min));
    ui_->max_const_->setText(QString("max: ") + QString::number(max));
  }

  bool Plot1DGoToDialog::checked()
  {
    return ui_->clip_checkbox->checkState() == Qt::Checked;
  }

  void Plot1DGoToDialog::fixRange()
  {
    // load from GUI
    float min_mz = ui_->min_->text().toFloat();
    float max_mz = ui_->max_->text().toFloat();

    // ensure correct order of min and max
    if (min_mz > max_mz) swap(min_mz, max_mz);

    // do not allow range of 0 --> extend to 1
    if (min_mz == max_mz)
    {
      min_mz -= 0.5;
      max_mz += 0.5;
    }

    // store in GUI
    ui_->min_->setText(QString::number(min_mz));
    ui_->max_->setText(QString::number(max_mz));
  }

  float Plot1DGoToDialog::getMin() const
  {
    return ui_->min_->text().toFloat();
  }

  float Plot1DGoToDialog::getMax() const
  {
    return ui_->max_->text().toFloat();
  }

} //namespace OpenMS
