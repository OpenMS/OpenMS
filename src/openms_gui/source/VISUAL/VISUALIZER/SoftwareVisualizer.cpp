// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

  SoftwareVisualizer::SoftwareVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Software>()
  {
    addLabel_("Modify software information.");
    addSeparator_();
    addLineEdit_(software_name_, "Name");
    addLineEdit_(software_version_, "Version");

    finishAdding_();
  }

  void SoftwareVisualizer::update_()
  {
    software_name_->setText(temp_.getName().c_str());
    software_version_->setText(temp_.getVersion().c_str());
  }

  void SoftwareVisualizer::store()
  {
    ptr_->setName(software_name_->text());
    ptr_->setVersion(software_version_->text());

    temp_ = (*ptr_);
  }

  void SoftwareVisualizer::undo_()
  {
    update_();
  }

}
