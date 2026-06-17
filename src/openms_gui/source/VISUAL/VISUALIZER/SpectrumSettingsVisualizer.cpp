// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

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
      fillComboBox_(type_, &temp_.NamesOfSpectrumType[static_cast<size_t>(temp_.getType())], 1);
    }
    else
    {
      fillComboBox_(type_, temp_.NamesOfSpectrumType, static_cast<int>(SpectrumSettings::SpectrumType::SIZE_OF_SPECTRUMTYPE));
      type_->setCurrentIndex(static_cast<int>(temp_.getType()));
    }

    native_id_->setText(temp_.getNativeID().c_str());
    comment_->setText(temp_.getComment().c_str());
  }

  void SpectrumSettingsVisualizer::store()
  {
    ptr_->setType((SpectrumSettings::SpectrumType)type_->currentIndex());
    ptr_->setNativeID(fromQString(native_id_->text()));
    ptr_->setComment(fromQString(comment_->toPlainText()));

    temp_ = (*ptr_);
  }

  void SpectrumSettingsVisualizer::undo_()
  {
    update_();
  }

}
