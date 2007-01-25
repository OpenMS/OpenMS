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
// $Maintainer:  stefan_heess $
// --------------------------------------------------------------------------s


#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Sample.h>



//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qcombobox.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qlistview.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <qstring.h>
#include <qvalidator.h>

#include <iostream>
#include <vector>
#include <string>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
SampleVisualizer::SampleVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="Sample";
  
	addLabel("Modify Sample information");		
	addSeperator();
  addLineEdit(samplename_, "Name" );
	addLineEdit(samplenumber_, "Number" );
	addLineEdit(sampleorganism_, "Organism" );
  addTextEdit(samplecomment_, "Comment");
	
	addComboBox(samplestate_, "State");
	
	addEmptyLine();
	
	addLineEdit(samplemass_,"Mass (in mg)");
	addLineEdit(samplevolume_, "Volume (in ml)");
	addLineEdit(sampleconcentration_, "Concentration (in mg/ml)");
	addVSpacer();
	addSeperator();
	addLabel("Save changes or restore original data");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
	// A validator to check the input for the mass
	QDoubleValidator *sample_massvali_ = new QDoubleValidator(samplemass_);
	samplemass_->setValidator(sample_massvali_);
	// A validator to check the input for the volume
	QDoubleValidator *volumevali_ = new QDoubleValidator(samplevolume_);
	samplevolume_->setValidator(volumevali_);
	// A Validator to check the input for the concentration
	QDoubleValidator *concentrationvali_ = new QDoubleValidator(sampleconcentration_);
	sampleconcentration_->setValidator(concentrationvali_);
	
	
}


void SampleVisualizer::load(Sample &s)
{
  //Pointer to current object	 to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempsample_=s;
		
	fillComboBox(samplestate_, s.NamesOfSampleState , Sample::SIZE_OF_SAMPLESTATE);
	
  update();		
}

void SampleVisualizer::update()
{
	samplename_->setText(tempsample_.getName());
	samplenumber_->setText(tempsample_.getNumber());
	sampleorganism_->setText(tempsample_.getOrganism());
  samplecomment_->setText(tempsample_.getComment());
	samplestate_->setCurrentItem(tempsample_.getState()); 
	samplemass_->setText(String( tempsample_.getMass())   );
	samplevolume_->setText(String( tempsample_.getVolume() ));
	sampleconcentration_->setText(String( tempsample_.getConcentration() ) );
	
}


void SampleVisualizer::store()
{
	try
	{
		(*ptr_).setName(string((const char*) samplename_->text()) );
		(*ptr_).setNumber(string((const char*) samplenumber_->text()) );
		(*ptr_).setOrganism(string((const char*) sampleorganism_->text()) );
		(*ptr_).setComment(string((const char*) samplecomment_->text()) );
				
		(*ptr_).setState((Sample::SampleState)samplestate_->currentItem() );		
	
		String m((const char*) samplemass_->text()) ;
		(*ptr_).setMass(m.toFloat() );
		
		String v((const char*) samplevolume_->text()) ;
		(*ptr_).setVolume(v.toFloat() );
		
		String c((const char*) sampleconcentration_->text()) ;
		(*ptr_).setConcentration(c.toFloat() );
		
		tempsample_=(*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new sample data. "<<e.what()<<endl;
	}
	
}

void SampleVisualizer::reject()
{
	
	try
	{
		//load(tempsample_);
		update();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original sample data. "<<e.what()<<endl;
	}
	
}

