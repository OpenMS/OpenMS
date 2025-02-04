// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s


#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QValidator>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SampleVisualizer::SampleVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Sample>()
  {
    addLabel_("Modify Sample information");
    addSeparator_();
    addLineEdit_(samplename_, "Name");
    addLineEdit_(samplenumber_, "Number");
    addLineEdit_(sampleorganism_, "Organism");
    addTextEdit_(samplecomment_, "Comment");
    addComboBox_(samplestate_, "State");
    addDoubleLineEdit_(samplemass_, "Mass (in gram)");
    addDoubleLineEdit_(samplevolume_, "Volume (in ml)");
    addDoubleLineEdit_(sampleconcentration_, "Concentration (in g/l)");

    finishAdding_();
  }

  void SampleVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(samplestate_, &temp_.NamesOfSampleState[temp_.getState()], 1);
    }
    else
    {
      fillComboBox_(samplestate_, temp_.NamesOfSampleState, Sample::SIZE_OF_SAMPLESTATE);
      samplestate_->setCurrentIndex(temp_.getState());
    }

    samplename_->setText(temp_.getName().c_str());
    samplenumber_->setText(temp_.getNumber().c_str());
    sampleorganism_->setText(temp_.getOrganism().c_str());
    samplecomment_->setText(temp_.getComment().c_str());

    samplemass_->setText(String(temp_.getMass()).c_str());
    samplevolume_->setText(String(temp_.getVolume()).c_str());
    sampleconcentration_->setText(String(temp_.getConcentration()).c_str());
  }

  void SampleVisualizer::store()
  {
    ptr_->setName(samplename_->text());
    ptr_->setNumber(samplenumber_->text());
    ptr_->setOrganism(sampleorganism_->text());
    ptr_->setComment(samplecomment_->toPlainText());
    ptr_->setState((Sample::SampleState)samplestate_->currentIndex());
    ptr_->setMass(samplemass_->text().toFloat());
    ptr_->setVolume(samplevolume_->text().toFloat());
    ptr_->setConcentration(sampleconcentration_->text().toFloat());

    temp_ = (*ptr_);
  }

  void SampleVisualizer::undo_()
  {
    update_();
  }

}
