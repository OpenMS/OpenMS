// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/DatabaseSequenceVisualizer.h>

#include <OpenMS/VISUAL/MISC/Qt5Port.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

namespace OpenMS
{

  DatabaseSequenceVisualizer::DatabaseSequenceVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IdentificationData::DatabaseSequence>()
  {
    addLineEdit_(accession_, "Accession");
    addLineEdit_(database_, "Database");
    addLineEdit_(target_decoy_, "Target/decoy");
    addTextEdit_(description_, "Description");
    addTextEdit_(sequence_, "Sequence");
    for (auto* field : {accession_, database_, target_decoy_}) field->setReadOnly(true);
    for (auto* field : {description_, sequence_}) field->setReadOnly(true);

    finishAdding_();
  }

  void DatabaseSequenceVisualizer::load(IdentificationData::DatabaseSequence & sequence, const IdentificationData::Run & run)
  {
    database_->setText(toQString(run.getDatabase(sequence.database).path));
    BaseVisualizer<IdentificationData::DatabaseSequence>::load(sequence);
  }

  void DatabaseSequenceVisualizer::update_()
  {
    accession_->setText(toQString(temp_.accession));
    switch (temp_.target_decoy)
    {
      case IdentificationData::TargetDecoy::TARGET: target_decoy_->setText("target"); break;
      case IdentificationData::TargetDecoy::DECOY: target_decoy_->setText("decoy"); break;
      case IdentificationData::TargetDecoy::BOTH: target_decoy_->setText("target+decoy"); break;
      default: target_decoy_->setText("unknown"); break;
    }
    description_->setText(toQString(temp_.description));
    sequence_->setText(toQString(temp_.sequence));
  }

  void DatabaseSequenceVisualizer::store()
  {
    // nothing to store: the sequence is shown, its metadata is edited in its own view
  }

  void DatabaseSequenceVisualizer::undo_()
  {
    update_();
  }

}
