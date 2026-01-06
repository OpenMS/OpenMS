// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>

//QT
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  IonSourceVisualizer::IonSourceVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IonSource>()
  {
    addLabel_("Modify ionsource information.");
    addSeparator_();

    addIntLineEdit_(order_, "Order");
    addComboBox_(inlet_type_, "Inlet type");
    addComboBox_(ionization_method_, "Ionization method");
    addComboBox_(polarity_, "Polarity");

    finishAdding_();
  }

  void IonSourceVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(inlet_type_, &temp_.NamesOfInletType[temp_.getInletType()], 1);
      fillComboBox_(ionization_method_, &temp_.NamesOfIonizationMethod[temp_.getIonizationMethod()], 1);
      fillComboBox_(polarity_, &temp_.NamesOfPolarity[temp_.getPolarity()], 1);
    }
    else
    {
      fillComboBox_(inlet_type_, temp_.NamesOfInletType, static_cast<int>(IonSource::InletType::SIZE_OF_INLETTYPE));
      fillComboBox_(ionization_method_, temp_.NamesOfIonizationMethod, static_cast<int>(IonSource::IonizationMethod::SIZE_OF_IONIZATIONMETHOD));
      fillComboBox_(polarity_, temp_.NamesOfPolarity, static_cast<int>(IonSource::Polarity::SIZE_OF_POLARITY));

      inlet_type_->setCurrentIndex(static_cast<int>(temp_.getInletType()));
      ionization_method_->setCurrentIndex(static_cast<int>(temp_.getIonizationMethod()));
      polarity_->setCurrentIndex(static_cast<int>(temp_.getPolarity()));
    }

    order_->setText(String(temp_.getOrder()).c_str());
  }

  void IonSourceVisualizer::store()
  {
    ptr_->setOrder(order_->text().toInt());
    ptr_->setInletType((IonSource::InletType)inlet_type_->currentIndex());
    ptr_->setIonizationMethod((IonSource::IonizationMethod)ionization_method_->currentIndex());
    ptr_->setPolarity((IonSource::Polarity)polarity_->currentIndex());

    temp_ = (*ptr_);
  }

  void IonSourceVisualizer::undo_()
  {
    update_();
  }

}
