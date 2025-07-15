// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/PeptideIdentificationVisualizer.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

//QT
#include <QtWidgets/QLineEdit>
#include <QValidator>
#include <QtWidgets/QPushButton>
#include <QtWidgets/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

  PeptideIdentificationVisualizer::PeptideIdentificationVisualizer(bool editable, QWidget * parent, MetaDataBrowser * caller) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<PeptideIdentification>()
  {
    pidv_caller_ = caller;

    addLineEdit_(identifier_, "Identifier<br>(of corresponding ProteinIdentification)");
    addSeparator_();

    addLineEdit_(score_type_, "Score type");
    addBooleanComboBox_(higher_better_, "Higher score is better");
    addDoubleLineEdit_(identification_threshold_, "Peptide significance threshold");

    addSeparator_();
    addLabel_("Show peptide hits with score equal or better than a threshold.");
    QPushButton * button;
    addLineEditButton_("Score threshold", filter_threshold_, button, "Filter");
    connect(button, SIGNAL(clicked()), this, SLOT(updateTree_()));

    finishAdding_();
  }

  void PeptideIdentificationVisualizer::load(PeptideIdentification & s, int tree_item_id)
  {
    ptr_ = &s;
    temp_ = s;

    // id of the item in the tree
    tree_id_ = tree_item_id;

    identifier_->setText(temp_.getIdentifier().toQString());
    identification_threshold_->setText(QString::number(temp_.getSignificanceThreshold()));
    
    // Score metadata is now centralized in PeptideIdentificationList
    // For individual objects, we show placeholder values and disable editing
    score_type_->setText("N/A (see list level)");
    score_type_->setEnabled(false);
    higher_better_->setCurrentIndex(0);
    higher_better_->setEnabled(false);
  }

  void PeptideIdentificationVisualizer::updateTree_()
  {
    if (filter_threshold_->text() != "")
    {
      // For filtering, we need to determine score orientation from context
      // Default to higher-is-better = true as a reasonable fallback
      pidv_caller_->filterHits_(filter_threshold_->text().toDouble(), true, tree_id_);
    }
    else
    {
      pidv_caller_->showAllHits_(tree_id_);
    }
  }

  void PeptideIdentificationVisualizer::store()
  {
    ptr_->setIdentifier(identifier_->text());
    ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
    
    // Score metadata (setScoreType, setHigherScoreBetter) is now managed at
    // PeptideIdentificationList level and cannot be set on individual objects
    // These operations are disabled in this context

    temp_ = (*ptr_);
  }

  void PeptideIdentificationVisualizer::undo_()
  {
    load(*ptr_, tree_id_);
  }

}
