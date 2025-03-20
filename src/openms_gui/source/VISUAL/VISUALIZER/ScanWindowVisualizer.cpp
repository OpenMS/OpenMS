// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/ScanWindowVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

  ScanWindowVisualizer::ScanWindowVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ScanWindow>()
  {
    addLabel_("Modify scan window information.");
    addSeparator_();
    addIntLineEdit_(begin_, "Begin");
    addIntLineEdit_(end_, "End");

    finishAdding_();
  }

  void ScanWindowVisualizer::update_()
  {
    begin_->setText(QString::number(temp_.begin));
    end_->setText(QString::number(temp_.end));
  }

  void ScanWindowVisualizer::store()
  {
    ptr_->begin = begin_->text().toDouble();
    ptr_->end = end_->text().toDouble();

    temp_ = (*ptr_);
  }

  void ScanWindowVisualizer::undo_()
  {
    update_();
  }

}
