// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/Plot2DGoToDialog.h>
#include <ui_Plot2DGoToDialog.h>


#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  Plot2DGoToDialog::Plot2DGoToDialog(QWidget * parent, std::string_view x_name, std::string_view y_name) :
    QDialog(parent),
    ui_(new Ui::Plot2DGoToDialogTemplate)
  {
    ui_->setupUi(this);
    ui_->dimx_->setText(x_name.data());
    ui_->dimy_->setText(y_name.data());
  }

  Plot2DGoToDialog::~Plot2DGoToDialog()
  {
    delete ui_;
  }

  void Plot2DGoToDialog::setRange(const AreaXYType& range)
  {
    ui_->min_x_->setText(QString::number(range.minX()));
    ui_->max_x_->setText(QString::number(range.maxX()));
    ui_->min_y_->setText(QString::number(range.minY()));
    ui_->max_y_->setText(QString::number(range.maxY()));
  }

  void Plot2DGoToDialog::setMinMaxOfRange(const AreaXYType& max_range)
  {
    ui_->min_x_const_->setText("min: " + QString::number(max_range.minX()));
    ui_->max_x_const_->setText("max: " + QString::number(max_range.maxX()));
    ui_->min_y_const_->setText("min: " + QString::number(max_range.minY()));
    ui_->max_y_const_->setText("max: " + QString::number(max_range.maxY()));
  }

  Plot2DGoToDialog::AreaXYType Plot2DGoToDialog::getRange()
  {
    AreaXYType r{
      ui_->min_x_->text().toFloat(),
      ui_->min_y_->text().toFloat(),
      ui_->max_x_->text().toFloat(),
      ui_->max_y_->text().toFloat(),
    };
    r.ensureMinSpan({1,1});
    return r;
  }

  bool Plot2DGoToDialog::checked()
  {
    return ui_->clip_checkbox->checkState() == Qt::Checked;
  }

  void Plot2DGoToDialog::enableFeatureNumber(bool enabled)
  {
    ui_->feature_label_->setEnabled(enabled);
    ui_->nr_->setEnabled(enabled);
    ui_->feature_number_->setEnabled(enabled);
    //Reorder tab order
    if (enabled)
    {
      setTabOrder(ui_->feature_number_, ui_->ok_button_);
      setTabOrder(ui_->ok_button_, ui_->cancel_button_);
      setTabOrder(ui_->cancel_button_, ui_->min_x_);
      setTabOrder(ui_->min_x_, ui_->max_x_);
      setTabOrder(ui_->max_x_, ui_->min_y_);
      setTabOrder(ui_->min_y_, ui_->max_y_);
    }
    else
    {
      setTabOrder(ui_->min_x_, ui_->max_x_);
      setTabOrder(ui_->max_x_, ui_->min_y_);
      setTabOrder(ui_->min_y_, ui_->max_y_);
      setTabOrder(ui_->max_y_, ui_->ok_button_);
      setTabOrder(ui_->ok_button_, ui_->cancel_button_);
    }
  }

  String Plot2DGoToDialog::getFeatureNumber() const
  {
    return ui_->feature_number_->text();
  }

  bool Plot2DGoToDialog::showRange() const
  {
    return getFeatureNumber().trim().empty();
  }

} //namespace OpenMS
