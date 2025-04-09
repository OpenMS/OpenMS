// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/ProductVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  ProductVisualizer::ProductVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Product>()
  {
    addLabel_("Modify processing method information.");

    addSeparator_();

    addDoubleLineEdit_(product_mz_, "m/z");
    addDoubleLineEdit_(product_window_low_, "Lower offset from target m/z");
    addDoubleLineEdit_(product_window_up_, "Upper offset from target m/z");

    finishAdding_();
  }

  void ProductVisualizer::update_()
  {
    product_mz_->setText(String(temp_.getMZ()).c_str());
    product_window_low_->setText(String(temp_.getIsolationWindowLowerOffset()).c_str());
    product_window_up_->setText(String(temp_.getIsolationWindowUpperOffset()).c_str());
  }

  void ProductVisualizer::store()
  {
    ptr_->setMZ(product_mz_->text().toFloat());
    ptr_->setIsolationWindowLowerOffset(product_window_low_->text().toFloat());
    ptr_->setIsolationWindowUpperOffset(product_window_up_->text().toFloat());

    temp_ = (*ptr_);
  }

  void ProductVisualizer::undo_()
  {
    update_();
  }

}
