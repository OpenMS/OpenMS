// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>

//QT
#include <QtWidgets/QTextEdit>
#include <QtWidgets/QComboBox>
#include <QtWidgets/QLineEdit>

#include <iostream>
#include <vector>

using namespace std;

namespace OpenMS
{

  ModificationVisualizer::ModificationVisualizer(bool editable, QWidget * parent) :
    BaseVisualizerGUI(editable, parent),
    BaseVisualizer<Modification>()
  {
    addLabel_("Modify Modification information");
    addSeparator_();
    addLineEdit_(treatmenttype_, "Treatment type");
    addTextEdit_(treatmentcomment_, "Comment");
    addLineEdit_(modificationname_, "Reagent name");
    addDoubleLineEdit_(modificationmass_, "Mass change");

    addComboBox_(modificationspecificity_, "Specificity Type");
    addLineEdit_(modificationAA_, "Affected Amino Acids");

    finishAdding_();
  }

  void ModificationVisualizer::update_()
  {
    if (!isEditable())
    {
      fillComboBox_(modificationspecificity_, &temp_.NamesOfSpecificityType[temp_.getSpecificityType()], 1);
    }
    else
    {
      fillComboBox_(modificationspecificity_, temp_.NamesOfSpecificityType, Modification::SIZE_OF_SPECIFICITYTYPE);
      modificationspecificity_->setCurrentIndex(temp_.getSpecificityType());
    }
    treatmenttype_->setText(temp_.getType().c_str());
    treatmenttype_->setReadOnly(true);
    treatmentcomment_->setText(temp_.getComment().c_str());
    modificationname_->setText(temp_.getReagentName().c_str());
    modificationmass_->setText(String(temp_.getMass()).c_str());
    modificationAA_->setText(temp_.getAffectedAminoAcids().c_str());
  }

  void ModificationVisualizer::store()
  {
    try
    {
      ptr_->setComment(treatmentcomment_->toPlainText());
      ptr_->setReagentName(modificationname_->text());
      ptr_->setMass(modificationmass_->text().toFloat());
      ptr_->setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());
      ptr_->setAffectedAminoAcids(modificationAA_->text());
      temp_ = (*ptr_);
    }
    catch (exception & e)
    {
      std::cout << "Error while trying to store the new modification data. " << e.what() << endl;
    }
  }

  void ModificationVisualizer::undo_()
  {
    update_();
  }

}
