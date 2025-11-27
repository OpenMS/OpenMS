// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/IonDetectorVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>


//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  IonDetectorVisualizer::IonDetectorVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IonDetector>()
  {
    addLabel_("Modify iondetector information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(type_, "Type");
    addComboBox_(ac_mode_, "Acquisition mode");
    addDoubleLineEdit_(res_, "Resolution (in ns)");
    addDoubleLineEdit_(freq_, "ADC sampling frequency (in Hz)");

    finishAdding_();
  }

  void IonDetectorVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfType[temp_.getType()], 1);
      fillComboBox_(ac_mode_, &temp_.NamesOfAcquisitionMode[temp_.getAcquisitionMode()], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfType, IonDetector::SIZE_OF_TYPE);
      fillComboBox_(ac_mode_, temp_.NamesOfAcquisitionMode, IonDetector::SIZE_OF_ACQUISITIONMODE);
      type_->setCurrentIndex(temp_.getType());
      ac_mode_->setCurrentIndex(temp_.getAcquisitionMode());
    }

    order_->setText(String(temp_.getOrder()).c_str());
    res_->setText(String(temp_.getResolution()).c_str());
    freq_->setText(String(temp_.getADCSamplingFrequency()).c_str());
  }

  void IonDetectorVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setResolution(res_->text().toDouble());
    ptr_->setADCSamplingFrequency(freq_->text().toDouble());
    ptr_->setType((IonDetector::Type)type_->currentIndex());
    ptr_->setAcquisitionMode((IonDetector::AcquisitionMode)ac_mode_->currentIndex());

    temp_ = (*ptr_);
  }

  void IonDetectorVisualizer::undo_()
  {
    update_();
  }

}
