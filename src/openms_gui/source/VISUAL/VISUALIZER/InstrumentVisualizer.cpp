// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QComboBox>

//STL
#include <iostream>
#include <string>

using namespace std;

namespace OpenMS
{

  InstrumentVisualizer::InstrumentVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Instrument>()
  {
    addLabel_("Modify instrument information.");
    addSeparator_();
    addLineEdit_(name_, "Name");
    addLineEdit_(vendor_, "Vendor");
    addLineEdit_(model_, "Model");
    addTextEdit_(customizations_, "Customizations");
    addComboBox_(ion_optics_, "Ion optics");

    finishAdding_();
  }

  void InstrumentVisualizer::update_()
  {
    name_->setText(temp_.getName().c_str());
    vendor_->setText(temp_.getVendor().c_str());
    model_->setText(temp_.getModel().c_str());
    customizations_->setText(temp_.getCustomizations().c_str());

    if (!isEditable())
    {
      fillComboBox_(ion_optics_, &temp_.NamesOfIonOpticsType[temp_.getIonOptics()], 1);
    }
    else
    {
      fillComboBox_(ion_optics_, temp_.NamesOfIonOpticsType, Instrument::SIZE_OF_IONOPTICSTYPE);
    }
  }

  void InstrumentVisualizer::store()
  {
    ptr_->setName(name_->text());
    ptr_->setVendor(vendor_->text());
    ptr_->setModel(model_->text());
    ptr_->setCustomizations(customizations_->toPlainText());
    ptr_->setIonOptics((Instrument::IonOpticsType)ion_optics_->currentIndex());

    temp_ = (*ptr_);
  }

  void InstrumentVisualizer::undo_()
  {
    update_();
  }

}
