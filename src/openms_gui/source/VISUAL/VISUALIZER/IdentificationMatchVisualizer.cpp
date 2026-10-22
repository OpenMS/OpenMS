// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/IdentificationMatchVisualizer.h>

#include <OpenMS/VISUAL/MISC/Qt5Port.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QTextEdit>

namespace OpenMS
{

  IdentificationMatchVisualizer::IdentificationMatchVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IdentificationData::MatchData>()
  {
    addTextEdit_(molecule_, "Molecule");
    addLineEdit_(charge_, "Charge");
    addLineEdit_(target_decoy_, "Target/decoy");
    addTextEdit_(scores_, "Scores");
    addTextEdit_(evidence_, "Sequence evidence");
    for (auto* field : {charge_, target_decoy_}) field->setReadOnly(true);
    for (auto* field : {molecule_, scores_, evidence_}) field->setReadOnly(true);

    finishAdding_();
  }

  void IdentificationMatchVisualizer::load(IdentificationData::MatchData & data, const IdentificationData::Run & run, const IdentificationData::Match & match)
  {
    QStringList scores;
    const auto values = run.getScores(match);
    const auto primary = run.getPrimaryScore();
    for (Size i = 0; i < values.size(); ++i)
    {
      QString text = toQString(run.getScoreDefinitions()[i].name) + ": " + (values[i] ? QString::number(*values[i]) : QString("-"));
      if (primary && primary->value == i) text += " (primary)";
      scores << text;
    }
    scores_->setText(scores.join("\n"));
    QStringList evidence;
    for (const auto& item : data.sequence_evidence)
    {
      QString text = toQString(item.accession.empty() ? std::string("(unknown entry)") : item.accession);
      if (item.start || item.end)
      {
        text += " " + (item.start ? QString::number(*item.start) : QString("?")) + "-" + (item.end ? QString::number(*item.end) : QString("?"));
      }
      evidence << text;
    }
    evidence_->setText(evidence.join("\n"));
    BaseVisualizer<IdentificationData::MatchData>::load(data);
  }

  void IdentificationMatchVisualizer::update_()
  {
    molecule_->setText(toQString(temp_.representation));
    charge_->setText(QString::number(temp_.charge));
    switch (temp_.target_decoy)
    {
      case IdentificationData::TargetDecoy::TARGET: target_decoy_->setText("target"); break;
      case IdentificationData::TargetDecoy::DECOY: target_decoy_->setText("decoy"); break;
      case IdentificationData::TargetDecoy::BOTH: target_decoy_->setText("target+decoy"); break;
      default: target_decoy_->setText("unknown"); break;
    }
  }

  void IdentificationMatchVisualizer::store()
  {
    // nothing to store: the match is shown, its metadata is edited in its own view
  }

  void IdentificationMatchVisualizer::undo_()
  {
    update_();
  }

}
