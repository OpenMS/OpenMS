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
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/TaggingVisualizer.h>
#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/METADATA/Tagging.h>

//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>


using namespace OpenMS;
using namespace std;

//Constructor
TaggingVisualizer::TaggingVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="Tagging";
	
	addLabel("Modify Tagging information");		
	addSeperator();
	addLineEdit(treatmenttype_, "Treatment type" );
	addLineEdit(modificationname_, "Reagent name" );
	addLineEdit(modificationmass_, "Mass" );
	 
	addComboBox(modificationspecificity_, "Specificity Type");
	addLineEdit(modificationAA_, "Affected Amino Acids" );
	
	addLineEdit(taggingmass_shift_, "Mass_Shift" );
	addComboBox(taggingvariant_, "Variant");
	addVSpacer();	
	addSeperator();
	addLabel("Save changes or restore original data");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
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
	
	
	updateTag();
}

void TaggingVisualizer::updateTag()
{
	modificationname_->setText(temptag_.getReagentName());
	modificationmass_->setText(String(temptag_.getMass()) );
	modificationspecificity_->setCurrentItem(temptag_.getSpecificityType());
	modificationAA_->setText(temptag_.getAffectedAminoAcids() ); 
	taggingmass_shift_->setText(String(temptag_.getMassShift()));
	taggingvariant_->setCurrentItem(temptag_.getVariant());			

}

void TaggingVisualizer::store()
{
	try
	{
		
		(*ptr_).setReagentName(string((const char*) modificationname_->text()));
				
		String m((const char*) modificationmass_->text()) ;
		(*ptr_).setMass(m.toFloat() );
				
		(*ptr_).setSpecificityType((Modification::SpecificityType)modificationspecificity_->currentItem());		
		
		(*ptr_).setAffectedAminoAcids(string((const char*) modificationAA_->text()) );
		
		(*ptr_).setMassShift(String((const char*)taggingmass_shift_->text()).toFloat() );
		
		(*ptr_).setVariant((Tagging::IsotopeVariant)taggingvariant_->currentItem());		
		
		
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
		updateTag();
		//load(temptag_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original tagging data. "<<e.what()<<endl;
	} 
}
