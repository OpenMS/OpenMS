// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

  PeptideHitVisualizer::PeptideHitVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<PeptideHit>()
  {
    addLineEdit_(peptidehit_score_, "Score");
    addLineEdit_(peptidehit_charge_, "Charge");
    addLineEdit_(peptidehit_rank_, "Rank");
    addTextEdit_(peptidehit_sequence_, "Sequence");

    finishAdding_();
  }

  void PeptideHitVisualizer::update_()
  {
    peptidehit_score_->setText(String(temp_.getScore()).c_str());
    peptidehit_score_->setReadOnly(true);
    peptidehit_charge_->setText(String(temp_.getCharge()).c_str());
    peptidehit_charge_->setReadOnly(true);
    peptidehit_rank_->setText(String(temp_.getRank()).c_str());
    peptidehit_rank_->setReadOnly(true);
    peptidehit_sequence_->setText(temp_.getSequence().toString().c_str());
    peptidehit_sequence_->setReadOnly(true);
  }

  void PeptideHitVisualizer::store()
  {
    //TODO?
    (*ptr_) = temp_;
  }

  void PeptideHitVisualizer::undo_()
  {
    update_();
  }

}
