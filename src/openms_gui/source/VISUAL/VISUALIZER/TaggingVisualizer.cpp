// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>

#include <QtWidgets/QComboBox>
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QLineEdit>
#include <QValidator>

#include <iostream>

using namespace std;

namespace OpenMS
{

  TaggingVisualizer::TaggingVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Tagging>()
  {
    addLabel_("Modify Tagging information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(modificationname_, "Reagent name");
    addDoubleLineEdit_(modificationmass_, "Mass");

    addComboBox_(modificationspecificity_, "Specificity Type");
    addLineEdit_(modificationAA_, "Affected Amino Acids");

    addDoubleLineEdit_(taggingmass_shift_, "Mass_Shift");
    addComboBox_(taggingvariant_, "Variant");

    finishAdding_();
  }

  void TaggingVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(modificationspecificity_, &temp_.NamesOfSpecificityType[temp_.getSpecificityType()], 1);
      fillComboBox_(taggingvariant_, &temp_.NamesOfIsotopeVariant[temp_.getVariant()], 1);
    }
    else
    {
      fillComboBox_(modificationspecificity_, temp_.NamesOfSpecificityType, Tagging::SIZE_OF_SPECIFICITYTYPE);
      fillComboBox_(taggingvariant_, temp_.NamesOfIsotopeVariant, Tagging::SIZE_OF_ISOTOPEVARIANT);
      modificationspecificity_->setCurrentIndex(temp_.getSpecificityType());
      taggingvariant_->setCurrentIndex(temp_.getVariant());
    }
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    modificationname_->setText(temp_.getReagentName().c_str());
    modificationmass_->setText(String(temp_.getMass()).c_str());

    modificationAA_->setText(temp_.getAffectedAminoAcids().c_str());
    taggingmass_shift_->setText(String(temp_.getMassShift()).c_str());
  }

  void TaggingVisualizer::store()
  {
    ptr_->setComment(treatmentcomment_->toPlainText());
    ptr_->setReagentName(modificationname_->text());
    ptr_->setMass(modificationmass_->text().toDouble());
    ptr_->setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());
    ptr_->setAffectedAminoAcids(modificationAA_->text());
    ptr_->setMassShift(taggingmass_shift_->text().toFloat());
    ptr_->setVariant((Tagging::IsotopeVariant)taggingvariant_->currentIndex());

    temp_ = (*ptr_);
  }

  void TaggingVisualizer::undo_()
  {
    update_();
  }

}
