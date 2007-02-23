// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free software; you can redistribute it and/or
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
// $Maintainer: stefan_heess  $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>

//QT

#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>

#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
TaggingVisualizer::TaggingVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
	type_="Tagging";
	
	addLabel("Modify Tagging information");		
	addSeperator();
	addLineEdit(treatmenttype_, "Treatment type" );
	addTextEdit(treatmentcomment_, "Comment" );
	addLineEdit(modificationname_, "Reagent name" );
	addLineEdit(modificationmass_, "Mass" );
	 
	addComboBox(modificationspecificity_, "Specificity Type");
	addLineEdit(modificationAA_, "Affected Amino Acids" );
	
	addLineEdit(taggingmass_shift_, "Mass_Shift" );
	addComboBox(taggingvariant_, "Variant");
	
	finishAdding_();
	
	massvali_ = new QDoubleValidator(modificationmass_);
	modificationmass_->setValidator(massvali_);
	
	shiftvali_ = new QDoubleValidator(taggingmass_shift_);
	taggingmass_shift_->setValidator(shiftvali_);
}


void TaggingVisualizer::load(Tagging &t)
{
  ptr_ = &t;
	
	//Copy of current object for restoring the original values
	temptag_=t;
  
	fillComboBox(modificationspecificity_, t.NamesOfSpecificityType, Tagging::SIZE_OF_SPECIFICITYTYPE);
  fillComboBox(taggingvariant_, t.NamesOfIsotopeVariant, Tagging::SIZE_OF_ISOTOPEVARIANT);
	
	
	updateTag_();
}

void TaggingVisualizer::updateTag_()
{
	treatmenttype_->setText(temptag_.getType().c_str());
	treatmenttype_->setReadOnly(true);
	treatmentcomment_->setText(temptag_.getComment().c_str());
	modificationname_->setText(temptag_.getReagentName().c_str());
	modificationmass_->setText(String(temptag_.getMass()).c_str() );
	modificationspecificity_->setCurrentIndex(temptag_.getSpecificityType());
	modificationAA_->setText(temptag_.getAffectedAminoAcids().c_str() ); 
	taggingmass_shift_->setText(String(temptag_.getMassShift()).c_str());
	taggingvariant_->setCurrentIndex(temptag_.getVariant());			

}

void TaggingVisualizer::store()
{
	try
	{
		(*ptr_).setComment(treatmentcomment_->toPlainText().toStdString());
		(*ptr_).setReagentName(modificationname_->text().toStdString());
				
		String m(modificationmass_->text().toStdString());
		(*ptr_).setMass(m.toFloat());
				
		(*ptr_).setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());
		(*ptr_).setAffectedAminoAcids(modificationAA_->text().toStdString());
		(*ptr_).setMassShift(taggingmass_shift_->text().toFloat());
		(*ptr_).setVariant((Tagging::IsotopeVariant)taggingvariant_->currentIndex());
		
		temptag_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new tagging data. "<<e.what()<<endl;
	}  
}

void TaggingVisualizer::reject()
{
	try
	{
		updateTag_();
		
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original tagging data. "<<e.what()<<endl;
	} 
}

}
