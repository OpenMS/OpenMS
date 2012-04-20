// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  this library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>

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
