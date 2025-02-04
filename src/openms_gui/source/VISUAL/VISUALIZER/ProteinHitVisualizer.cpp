// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer:Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

  ProteinHitVisualizer::ProteinHitVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ProteinHit>()
  {
    addLineEdit_(proteinhit_score_, "Score");
    addLineEdit_(proteinhit_rank_, "Rank");
    addLineEdit_(proteinhit_accession_, "Accession");
    addTextEdit_(proteinhit_sequence_, "Sequence");

    finishAdding_();
  }

  void ProteinHitVisualizer::update_()
  {
    proteinhit_score_->setText(String(temp_.getScore()).c_str());
    proteinhit_score_->setReadOnly(true);
    proteinhit_rank_->setText(String(temp_.getRank()).c_str());
    proteinhit_rank_->setReadOnly(true);
    proteinhit_accession_->setText(temp_.getAccession().c_str());
    proteinhit_accession_->setReadOnly(true);
    proteinhit_sequence_->setText(temp_.getSequence().c_str());
    proteinhit_sequence_->setReadOnly(true);
  }

  void ProteinHitVisualizer::store()
  {
    //TODO?
    (*ptr_) = temp_;
  }

  void ProteinHitVisualizer::undo_()
  {
    update_();
  }

}
