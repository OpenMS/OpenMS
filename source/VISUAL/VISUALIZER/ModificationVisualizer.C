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

//Constructor
ModificationVisualizer::ModificationVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="Modification";
	
	addLabel("Modify Modification information");		
	addSeperator();
	addLineEdit(treatmenttype_, "Treatment type" );
	addTextEdit(treatmentcomment_, "Comment" );
	addLineEdit(modificationname_, "Reagent name" );
	addLineEdit(modificationmass_, "Mass change" );
	 
	addComboBox(modificationspecificity_, "Specificity Type");
	addLineEdit(modificationAA_, "Affected Amino Acids" );
	
	finishAdding_();
	
	// A Validator to check the input for the mass
	QDoubleValidator *massvali_ = new QDoubleValidator(modificationmass_);
	modificationmass_->setValidator(massvali_);
}


void ModificationVisualizer::load(Modification &m)
{
  ptr_ = &m;
	
	//Copy of current object for restoring the original values
	tempmod_=m;
	
	
	fillComboBox(modificationspecificity_, m.NamesOfSpecificityType, Modification::SIZE_OF_SPECIFICITYTYPE);
  
	updateMod_();
			
}

void ModificationVisualizer::updateMod_()
{
	treatmenttype_->setText(tempmod_.getType().c_str());
	treatmenttype_->setReadOnly(true);
	treatmentcomment_->setText(tempmod_.getComment().c_str());
  modificationname_->setText(tempmod_.getReagentName().c_str());
	modificationmass_->setText(String(tempmod_.getMass()).c_str() );
	modificationspecificity_->setCurrentIndex(tempmod_.getSpecificityType());
	modificationAA_->setText(tempmod_.getAffectedAminoAcids().c_str() ); 

}


void ModificationVisualizer::store()
{
try
	{
		(*ptr_).setComment(treatmentcomment_->toPlainText().toStdString());
		
		(*ptr_).setReagentName(modificationname_->text().toStdString());
				
		String m(modificationmass_->text().toStdString()) ;
		
				
		(*ptr_).setMass(m.toFloat() );
		
		(*ptr_).setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentIndex());		
		
				
		(*ptr_).setAffectedAminoAcids(modificationAA_->text().toStdString());
		
		tempmod_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new modification data. "<<e.what()<<endl;
	}
	  
}

void ModificationVisualizer::reject()
{
	try
	{
		updateMod_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original modification data. "<<e.what()<<endl;
	}
}

}
