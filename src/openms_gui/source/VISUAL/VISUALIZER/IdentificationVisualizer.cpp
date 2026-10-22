// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/IdentificationVisualizer.h>

#include <OpenMS/VISUAL/MetaDataBrowser.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QtWidgets/QPushButton>

namespace OpenMS
{

  IdentificationVisualizer::IdentificationVisualizer(bool editable, QWidget * parent, MetaDataBrowser * caller) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<IdentificationData::Observation>(),
    caller_(caller)
  {
    addLineEdit_(run_, "Run");
    addLineEdit_(data_id_, "Data reference<br>(e.g. spectrum)");
    addLineEdit_(rt_, "Retention time");
    addLineEdit_(mz_, "m/z");
    addLineEdit_(score_, "Primary score");
    for (auto* field : {run_, data_id_, rt_, mz_, score_}) field->setReadOnly(true);

    addSeparator_();
    addLabel_("Show matches with score equal or better than a threshold.");
    QPushButton * button;
    addLineEditButton_("Score threshold", filter_threshold_, button, "Filter");
    // filtering the view is no edit
    filter_threshold_->setReadOnly(false);
    button->setEnabled(true);
    connect(button, &QPushButton::clicked, this, &IdentificationVisualizer::updateTree_);

    finishAdding_();
  }

  void IdentificationVisualizer::load(IdentificationData::Observation & observation, const IdentificationData::Run & run, int tree_item_id)
  {
    tree_id_ = tree_item_id;
    run_->setText(toQString(run.getIdentifier()));
    if (const auto primary = run.getPrimaryScore())
    {
      const auto& definition = run.getScoreDefinition(*primary);
      higher_better_ = definition.higher_better;
      score_->setText(toQString(definition.name) + (higher_better_ ? " (higher is better)" : " (lower is better)"));
    }
    BaseVisualizer<IdentificationData::Observation>::load(observation);
  }

  void IdentificationVisualizer::update_()
  {
    data_id_->setText(toQString(temp_.data_id));
    rt_->setText(temp_.rt ? QString::number(*temp_.rt) : QString());
    mz_->setText(temp_.mz ? QString::number(*temp_.mz) : QString());
  }

  void IdentificationVisualizer::updateTree_()
  {
    if (caller_ == nullptr) return;
    if (filter_threshold_->text() != "")
    {
      caller_->filterHits_(filter_threshold_->text().toDouble(), higher_better_, tree_id_);
    }
    else
    {
      caller_->showAllHits_(tree_id_);
    }
  }

  void IdentificationVisualizer::store()
  {
    // nothing to store: the observation is shown, its metadata is edited in its own view
  }

  void IdentificationVisualizer::undo_()
  {
    update_();
  }

}
