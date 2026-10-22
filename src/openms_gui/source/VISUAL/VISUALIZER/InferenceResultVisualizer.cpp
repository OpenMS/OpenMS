// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/InferenceResultVisualizer.h>

#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

//QT
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QTextEdit>

namespace OpenMS
{

  InferenceResultVisualizer::InferenceResultVisualizer(bool editable, QWidget * parent, MetaDataBrowser * caller) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<ProteinIdentification>(),
    caller_(caller)
  {
    addLineEdit_(identifier_, "Identifier");
    addTextEdit_(inputs_, "Runs");
    identifier_->setReadOnly(true);
    inputs_->setReadOnly(true);
    addSeparator_();
    addLineEdit_(score_type_, "Score type");
    addBooleanComboBox_(higher_better_, "Higher score is better");
    addDoubleLineEdit_(significance_threshold_, "Protein significance threshold");
    score_type_->setObjectName("score_type");

    addSeparator_();
    addLabel_("Show protein hits with score equal or better than a threshold.");
    QPushButton * button;
    addLineEditButton_("Score threshold", filter_threshold_, button, "Filter");
    // filtering the view is no edit
    filter_threshold_->setReadOnly(false);
    button->setEnabled(true);
    connect(button, &QPushButton::clicked, this, &InferenceResultVisualizer::updateTree_);

    finishAdding_();
  }

  void InferenceResultVisualizer::load(IdentificationData::InferenceResult & result, int tree_item_id)
  {
    tree_id_ = tree_item_id;
    identifier_->setText(toQString(result.identifier));
    QStringList inputs;
    for (const auto& input : result.inputs) inputs << toQString(input.run_identifier);
    inputs_->setText(inputs.join("\n"));
    BaseVisualizer<ProteinIdentification>::load(result.proteins);
  }

  void InferenceResultVisualizer::update_()
  {
    score_type_->setText(toQString(temp_.getScoreType()));
    higher_better_->setCurrentIndex(temp_.isHigherScoreBetter());
    significance_threshold_->setText(QString::number(temp_.getSignificanceThreshold()));
  }

  void InferenceResultVisualizer::updateTree_()
  {
    if (caller_ == nullptr) return;
    if (filter_threshold_->text() != "")
    {
      caller_->filterHits_(filter_threshold_->text().toDouble(), temp_.isHigherScoreBetter(), tree_id_);
    }
    else
    {
      caller_->showAllHits_(tree_id_);
    }
  }

  void InferenceResultVisualizer::store()
  {
    ptr_->setScoreType(fromQString(score_type_->text()));
    ptr_->setHigherScoreBetter(higher_better_->currentIndex());
    // (written only if its text changed, so a value that the text shows rounded stays exact)
    if (significance_threshold_->text() != QString::number(temp_.getSignificanceThreshold()))
    {
      ptr_->setSignificanceThreshold(significance_threshold_->text().toDouble());
    }

    temp_ = (*ptr_);
  }

  void InferenceResultVisualizer::undo_()
  {
    update_();
  }

}
