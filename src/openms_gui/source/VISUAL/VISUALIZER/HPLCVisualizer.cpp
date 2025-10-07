// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>
#include <QValidator>
#include <iostream>

using namespace std;

namespace OpenMS
{

  HPLCVisualizer::HPLCVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<HPLC>()
  {
    addLabel_("Modify HPLC information");
    addSeparator_();
    addLineEdit_(hplcinstrument_, "Instrument");
    addLineEdit_(hplccolumn_, "Column");
    addIntLineEdit_(hplctemperature_, "Temperature (in deg. C)");
    addIntLineEdit_(hplcpressure_, "Pressure (in bar)");
    addIntLineEdit_(hplcflux_, "Flux (in ul/sec)");
    addTextEdit_(hplccomment_, "Comment");

    finishAdding_();
  }

  void HPLCVisualizer::update_()
  {
    hplcinstrument_->setText(temp_.getInstrument().c_str());
    hplccolumn_->setText(temp_.getColumn().c_str());
    hplctemperature_->setText(String(temp_.getTemperature()).c_str());
    hplcpressure_->setText(String(temp_.getPressure()).c_str());
    hplcflux_->setText(String(temp_.getFlux()).c_str());
    hplccomment_->setText(temp_.getComment().c_str());
  }

  void HPLCVisualizer::store()
  {
    ptr_->setInstrument(hplcinstrument_->text());
    ptr_->setColumn(hplccolumn_->text());
    ptr_->setTemperature(hplctemperature_->text().toInt());
    ptr_->setPressure(hplcpressure_->text().toInt());
    ptr_->setFlux(hplcflux_->text().toInt());
    ptr_->setComment(hplccomment_->toPlainText());
    temp_ = (*ptr_);
  }

  void HPLCVisualizer::undo_()
  {
    update_();
  }

}
