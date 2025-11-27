// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QValidator>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  AcquisitionVisualizer::AcquisitionVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Acquisition>()
  {

    addLabel_("Show Acquisition information");
    addSeparator_();
    addIntLineEdit_(acquisitionnumber_, "Identifier of the scan");
    acquisitionnumber_->setReadOnly(true);

    finishAdding_();
  }

  void AcquisitionVisualizer::update_()
  {
    acquisitionnumber_->setText(temp_.getIdentifier().toQString());
  }

  void AcquisitionVisualizer::store()
  {
    ptr_->setIdentifier(acquisitionnumber_->text());

    temp_ = (*ptr_);
  }

  void AcquisitionVisualizer::undo_()
  {
    update_();
  }

}
