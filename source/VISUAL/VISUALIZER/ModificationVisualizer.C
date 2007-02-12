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
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/VISUALIZER/ModificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/METADATA/Modification.h>

//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <qvalidator.h>
#include <iostream>
#include <vector>


//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
ModificationVisualizer::ModificationVisualizer(bool editable, QWidget *parent, const char *name) : BaseVisualizer(editable, parent, name)
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
	treatmenttype_->setText(tempmod_.getType());
	treatmenttype_->setReadOnly(true);
	treatmentcomment_->setText(tempmod_.getComment());
  modificationname_->setText(tempmod_.getReagentName());
	modificationmass_->setText(String(tempmod_.getMass()) );
	modificationspecificity_->setCurrentItem(tempmod_.getSpecificityType());
	modificationAA_->setText(tempmod_.getAffectedAminoAcids() ); 

}


void ModificationVisualizer::store()
{
try
	{
		(*ptr_).setComment(string((const char*) treatmentcomment_->text()));
		
		(*ptr_).setReagentName(string((const char*) modificationname_->text()));
				
		String m((const char*) modificationmass_->text()) ;
		
				
		(*ptr_).setMass(m.toFloat() );
		
		(*ptr_).setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentItem());		
		
				
		(*ptr_).setAffectedAminoAcids(string((const char*) modificationAA_->text()) );
		
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

