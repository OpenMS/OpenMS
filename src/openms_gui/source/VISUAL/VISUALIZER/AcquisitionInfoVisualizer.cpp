// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionInfoVisualizer.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

// QT
#include <QValidator>
#include <QtWidgets/QLineEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  AcquisitionInfoVisualizer::AcquisitionInfoVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<AcquisitionInfo>()
  {
    addLabel_("Show AcquisitionInfo information");
    addSeparator_();
    addIntLineEdit_(acquisitioninfo_method_, "Method of combination");

    finishAdding_();
  }

  void AcquisitionInfoVisualizer::update_()
  {
    acquisitioninfo_method_->setText(temp_.getMethodOfCombination().c_str());
  }

  void AcquisitionInfoVisualizer::store()
  {
    ptr_->setMethodOfCombination(fromQString(acquisitioninfo_method_->text()));

    temp_ = (*ptr_);
  }

  void AcquisitionInfoVisualizer::undo_()
  {
    update_();
  }

}
