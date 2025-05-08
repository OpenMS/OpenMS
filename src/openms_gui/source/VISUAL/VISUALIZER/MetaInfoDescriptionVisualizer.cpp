// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

using namespace std;

namespace OpenMS
{

  MetaInfoDescriptionVisualizer::MetaInfoDescriptionVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<MetaInfoDescription>()
  {
    addLabel_("Modify MetaInfoDescription information");
    addSeparator_();
    addLineEdit_(metainfodescription_name_, "Name of peak annotations");

    finishAdding_();
  }

  void MetaInfoDescriptionVisualizer::update_()
  {
    metainfodescription_name_->setText(temp_.getName().c_str());
  }

  void MetaInfoDescriptionVisualizer::store()
  {
    ptr_->setName(metainfodescription_name_->text().toStdString());

    temp_ = (*ptr_);
  }

  void MetaInfoDescriptionVisualizer::undo_()
  {
    update_();
  }

}
