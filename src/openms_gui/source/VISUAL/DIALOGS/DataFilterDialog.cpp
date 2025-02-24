// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/DIALOGS/DataFilterDialog.h>
#include <ui_DataFilterDialog.h>

//Qt includes
#include <QDoubleValidator>
#include <QIntValidator>
#include <QtWidgets/QMessageBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  DataFilterDialog::DataFilterDialog(DataFilters::DataFilter & filter, QWidget * parent) :
    QDialog(parent),
    filter_(filter),
    ui_(new Ui::DataFilterDialogTemplate)
  {
    ui_->setupUi(this);
    connect(ui_->ok_button_, SIGNAL(clicked()), this, SLOT(check_()));
    connect(ui_->field_, SIGNAL(activated(const QString &)), this, SLOT(field_changed_(const QString &)));
    connect(ui_->op_, SIGNAL(activated(const QString &)), this, SLOT(op_changed_(const QString &)));

    //set values for edit mode
    ui_->field_->setCurrentIndex((UInt)filter.field);
    ui_->op_->setCurrentIndex((UInt)filter.op);
    if (filter.field == DataFilters::META_DATA)
    {
      ui_->meta_name_field_->setText(filter.meta_name.toQString());
      // if the value stored in filter is numerical, get value from filter.value (a double)
      if (filter.value_is_numerical)
      {
        ui_->value_->setText(QString::number(filter.value));
      }
      else       // get value from filter.value_string (a String)
      {
        ui_->value_->setText(filter.value_string.toQString());
      }
      ui_->meta_name_field_->setEnabled(true);
      ui_->meta_name_label_->setEnabled(true);
      if (filter.op == DataFilters::EXISTS)
      {
        ui_->value_->setEnabled(false);
        ui_->value_label_->setEnabled(false);
      }
    }
    else     // for non meta data, the value is always numerical
    {
      ui_->value_->setText(QString::number(filter.value));
    }

    //focus the value if this is an edit operation
    if (filter != DataFilters::DataFilter())
    {
      ui_->value_->selectAll();
      setTabOrder(ui_->value_, ui_->cancel_button_);
      setTabOrder(ui_->cancel_button_, ui_->ok_button_);
      setTabOrder(ui_->ok_button_, ui_->field_);
      setTabOrder(ui_->field_, ui_->meta_name_field_);
      setTabOrder(ui_->meta_name_field_, ui_->op_);
    }
  }

  DataFilterDialog::~DataFilterDialog()
  {
    delete ui_;
  }

  void DataFilterDialog::field_changed_(const QString & field)
  {
    QString op(ui_->op_->currentText());
    if (field == "Meta data")
    {
      ui_->meta_name_field_->setEnabled(true);
      ui_->meta_name_label_->setEnabled(true);
    }
    else
    {
      ui_->meta_name_field_->setEnabled(false);
      ui_->meta_name_label_->setEnabled(false);
    }
  }

  void DataFilterDialog::op_changed_(const QString & op)
  {
    QString field(ui_->field_->currentText());
    if (op != "exists")
    {
      ui_->value_->setEnabled(true);
      ui_->value_label_->setEnabled(true);
    }
    else
    {
      ui_->value_->setEnabled(false);
      ui_->value_label_->setEnabled(false);
    }
  }

  void DataFilterDialog::check_()
  {
    QString field = ui_->field_->currentText();
    QString op = ui_->op_->currentText();
    QString value = ui_->value_->text();
    QString meta_name_field = ui_->meta_name_field_->text();
    bool not_numerical = true;
    int tmp;

    //meta data
    if (field == "Meta data")
    {
      QDoubleValidator dv(this);
      not_numerical = dv.validate(value, tmp) == QValidator::Invalid;

      if (meta_name_field.isEmpty())
      {
        QMessageBox::warning(this, "Insufficient arguments", "You must specify a meta name!");
        return;
      }
      if (op == "<=" || op == ">=")
      {
        if (not_numerical)
        {
          QMessageBox::warning(this, "Invalid value", "<= and >= are defined for numerical values only!");
          return;
        }
      }
    }
    //intensity, quality, charge:
    else
    {
      if (op == "exists")
      {
        QMessageBox::warning(this, "Invalid operation", "Operation \"exists\" is defined for meta data only!");
        return;
      }
      //double
      if (field == "Intensity" || field == "Quality")
      {
        QDoubleValidator v(this);
        if (v.validate(value, tmp) == QValidator::Invalid)
        {
          QMessageBox::warning(this, "Invalid value", "A real value is required!");
          return;
        }
      }
      //int
      if (field == "Charge" || field == "Size")
      {
        QIntValidator v(this);
        if (v.validate(value, tmp) == QValidator::Invalid)
        {
          QMessageBox::warning(this, "Invalid value", "An integer value is required!");
          return;
        }
      }
    }

    //write result
    if (field == "Intensity")
      filter_.field = DataFilters::INTENSITY;
    else if (field == "Quality")
      filter_.field = DataFilters::QUALITY;
    else if (field == "Charge")
      filter_.field = DataFilters::CHARGE;
    else if (field == "Size")
      filter_.field = DataFilters::SIZE;
    else if (field == "Meta data")
    {
      filter_.field = DataFilters::META_DATA;
      filter_.meta_name = meta_name_field;
      if (not_numerical)       // entered value is not numerical, store it in value_string (as String)
      {
        filter_.value_string = String(value);
        filter_.value_is_numerical = false;
      }
      else       // value is numerical, store it in value (as double)
      {
        filter_.value = value.toDouble();
        filter_.value_is_numerical = true;
      }
    }

    if (op == ">=")
      filter_.op = DataFilters::GREATER_EQUAL;
    else if (op == "=")
      filter_.op = DataFilters::EQUAL;
    else if (op == "<=")
      filter_.op = DataFilters::LESS_EQUAL;
    else if (op == "exists")
      filter_.op = DataFilters::EXISTS;

    if (field == "Intensity" || field == "Quality")
      filter_.value = value.toDouble();
    else if (field == "Charge" || field == "Size")
      filter_.value = value.toInt();

    accept();
  }

}
