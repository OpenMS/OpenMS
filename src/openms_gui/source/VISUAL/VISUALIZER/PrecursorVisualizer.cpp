// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QListWidget>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  PrecursorVisualizer::PrecursorVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Precursor>()
  {
    addLabel_("Modify processing method information.");

    addSeparator_();

    addDoubleLineEdit_(mz_, "m/z");
    addDoubleLineEdit_(int_, "intensity");
    addIntLineEdit_(charge_, "charge");

    addDoubleLineEdit_(window_low_, "Lower offset from target m/z");
    addDoubleLineEdit_(window_up_, "Upper offset from target m/z");

    addListView_(activation_methods_, "Activation methods");
    addDoubleLineEdit_(activation_energy_, "Activation energy");

    finishAdding_();
  }

  void PrecursorVisualizer::update_()
  {
    mz_->setText(String(temp_.getMZ()).c_str());
    int_->setText(String(temp_.getIntensity()).c_str());
    charge_->setText(String(temp_.getCharge()).c_str());

    window_low_->setText(String(temp_.getIsolationWindowLowerOffset()).c_str());
    window_up_->setText(String(temp_.getIsolationWindowUpperOffset()).c_str());

    //actions
    activation_methods_->clear();
    for (Size i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
    {
      QListWidgetItem * item = new QListWidgetItem(activation_methods_);
      item->setText(QString::fromStdString(Precursor::NamesOfActivationMethod[i]));
      if (temp_.getActivationMethods().count(Precursor::ActivationMethod(i)) == 1)
      {
        item->setCheckState(Qt::Checked);
      }
      else
      {
        item->setCheckState(Qt::Unchecked);
      }
      if (isEditable())
      {
        item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
      }
      else
      {
        item->setFlags(Qt::ItemIsEnabled);
      }
      activation_methods_->addItem(item);
    }

    activation_energy_->setText(String(temp_.getActivationEnergy()).c_str());
  }

  void PrecursorVisualizer::store()
  {
    ptr_->setMZ(mz_->text().toFloat());
    ptr_->setIntensity(int_->text().toFloat());
    ptr_->setCharge(charge_->text().toInt());

    ptr_->setIsolationWindowLowerOffset(window_low_->text().toFloat());
    ptr_->setIsolationWindowUpperOffset(window_up_->text().toFloat());

    ptr_->getActivationMethods().clear();
    for (UInt i = 0; i < Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
    {
      if (activation_methods_->item(i)->checkState() == Qt::Checked)
      {
        ptr_->getActivationMethods().insert(Precursor::ActivationMethod(i));
      }
    }
    ptr_->setActivationEnergy(activation_energy_->text().toFloat());

    temp_ = (*ptr_);
  }

  void PrecursorVisualizer::undo_()
  {
    update_();
  }

}
