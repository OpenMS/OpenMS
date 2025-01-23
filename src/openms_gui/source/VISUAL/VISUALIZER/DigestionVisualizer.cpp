// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>

//QT
#include <QValidator>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  DigestionVisualizer::DigestionVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Digestion>()
  {
    addLabel_("Modify Digestion information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(digestionenzyme_, "Enzyme");
    addDoubleLineEdit_(digestiontime_, "Digestion time (in min)");
    addDoubleLineEdit_(digestiontemperature_, "Temperature (in deg. C)");
    addDoubleLineEdit_(digestionPH_, "pH");

    finishAdding_();
  }

  void DigestionVisualizer::update_()
  {
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    digestionenzyme_->setText(temp_.getEnzyme().c_str());
    digestiontime_->setText(String(temp_.getDigestionTime()).c_str());
    digestiontemperature_->setText(String(temp_.getTemperature()).c_str());
    digestionPH_->setText(String(temp_.getPh()).c_str());
  }

  void DigestionVisualizer::store()
  {
    ptr_->setComment(treatmentcomment_->toPlainText());
    ptr_->setEnzyme(digestionenzyme_->text());
    ptr_->setDigestionTime(digestiontime_->text().toFloat());
    ptr_->setTemperature(digestiontime_->text().toFloat());
    ptr_->setPh(digestiontime_->text().toFloat());

    temp_ = (*ptr_);
  }

  void DigestionVisualizer::undo_()
  {
    update_();
  }

}
