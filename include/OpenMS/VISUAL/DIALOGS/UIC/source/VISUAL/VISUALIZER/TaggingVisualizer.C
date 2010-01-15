// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>

#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>

#include <iostream>

using namespace std;

namespace OpenMS
{

	TaggingVisualizer::TaggingVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Tagging>()
	{
		addLabel_("Modify Tagging information");		
		addSeparator_();
		addLineEdit_(treatmenttype_, "Treatment type" );
		addTextEdit_(treatmentcomment_, "Comment" );
		addLineEdit_(modificationname_, "Reagent name" );
		addDoubleLineEdit_(modificationmass_, "Mass" );
		 
		addComboBox_(modificationspecificity_, "Specificity Type");
		addLineEdit_(modificationAA_, "Affected Amino Acids" );
		
		addDoubleLineEdit_(taggingmass_shift_, "Mass_Shift" );
		addComboBox_(taggingvariant_, "Variant");
		
		finishAdding_();
	}
	
	void TaggingVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(modificationspecificity_,& temp_.NamesOfSpecificityType[temp_.getSpecificityType()], 1);
	  	fillComboBox_(taggingvariant_,& temp_.NamesOfIsotopeVariant[temp_.getVariant()], 1);
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
		modificationmass_->setText(String(temp_.getMass()).c_str() );
		
		modificationAA_->setText(temp_.getAffectedAminoAcids().c_str() ); 
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
