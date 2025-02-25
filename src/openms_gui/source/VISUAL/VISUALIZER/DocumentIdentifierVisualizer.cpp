// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/DocumentIdentifierVisualizer.h>
#include <OpenMS/FORMAT/FileHandler.h>

//QT
#include <QtWidgets/QLineEdit>

using namespace std;

namespace OpenMS
{

  DocumentIdentifierVisualizer::DocumentIdentifierVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<DocumentIdentifier>()
  {
    addLabel_("Modify DocumentIdentifier information");
    addSeparator_();
    addLineEdit_(identifier_, "Identifier");
    addSeparator_();
    addLineEdit_(file_path_, "Loaded from file");
    addLineEdit_(file_type_, "File type");
    finishAdding_();
  }

  void DocumentIdentifierVisualizer::update_()
  {
    identifier_->setText(temp_.getIdentifier().c_str());
    file_path_->setText(temp_.getLoadedFilePath().c_str());
    file_type_->setText(FileTypes::typeToName(temp_.getLoadedFileType()).c_str());
    file_path_->setReadOnly(true);
    file_type_->setReadOnly(true);
  }

  void DocumentIdentifierVisualizer::store()
  {
    ptr_->setIdentifier(identifier_->text());

    temp_ = (*ptr_);
  }

  void DocumentIdentifierVisualizer::undo_()
  {
    update_();
  }

}
