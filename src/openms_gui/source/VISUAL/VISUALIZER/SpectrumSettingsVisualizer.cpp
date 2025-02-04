// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

//QT
#include <QtWidgets/QComboBox>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  SpectrumSettingsVisualizer::SpectrumSettingsVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<SpectrumSettings>()
  {
    addLabel_("Modify the settings of the spectrum.");
    addSeparator_();
    addComboBox_(type_, "Type of spectrum");
    addLineEdit_(native_id_, "Native ID");
    addTextEdit_(comment_, "Comment");

    finishAdding_();
  }

  void SpectrumSettingsVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(type_, &temp_.NamesOfSpectrumType[temp_.getType()], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfSpectrumType, SpectrumSettings::SIZE_OF_SPECTRUMTYPE);
      type_->setCurrentIndex(temp_.getType());
    }

    native_id_->setText(temp_.getNativeID().c_str());
    comment_->setText(temp_.getComment().c_str());
  }

  void SpectrumSettingsVisualizer::store()
  {
    ptr_->setType((SpectrumSettings::SpectrumType)type_->currentIndex());
    ptr_->setNativeID(native_id_->text());
    ptr_->setComment(comment_->toPlainText());

    temp_ = (*ptr_);
  }

  void SpectrumSettingsVisualizer::undo_()
  {
    update_();
  }

}
