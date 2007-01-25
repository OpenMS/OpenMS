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
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/DataTable.h>
#include <OpenMS/METADATA/Digestion.h>

//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qaction.h>
#include <qcombobox.h>
#include <qfiledialog.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qmenubar.h>
#include <qmessagebox.h>
#include <qpopupmenu.h>
#include <qsettings.h>
#include <qstatusbar.h>
#include <qapplication.h>
#include <qlistview.h>
#include <qtextedit.h>
#include <qhbox.h>
#include <qgroupbox.h>
#include <qpushbutton.h>
#include <qvalidator.h>

//STL
#include <iostream>
#include <vector>


//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
DigestionVisualizer::DigestionVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
  type_="Digestion";
	
	addLabel("Modify Digestion information");		
	addSeperator();
	addLineEdit(treatmenttype_, "Treatment type" );
	addLineEdit(digestionenzyme_, "Enzyme" );
	addLineEdit(digestiontime_, "Digestion time (in minutes)" );
	addLineEdit(digestiontemperature_, "Temperature (in °C)" );
	addLineEdit(digestionPH_, "PH" );
	addVSpacer();	
	addSeperator();
	addLabel("Save changes or restore original data.");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
	// A validator to check the input for the time
	QDoubleValidator *timevali_ = new QDoubleValidator(digestiontime_);
	digestiontime_->setValidator(timevali_);
	// A validator to check the input for the temperature
	QDoubleValidator *tempvali_ = new QDoubleValidator(digestiontemperature_);
	digestiontemperature_->setValidator(tempvali_);
	// A validator to check the input for the ph value
	QDoubleValidator *phvali_ = new QDoubleValidator(digestionPH_);
	digestionPH_->setValidator(phvali_);
	
}


void DigestionVisualizer::load(Digestion &d)
{
  ptr_ = &d;
	
	//Copy of current object for restoring the original values
	tempdig_=d;
  digestionenzyme_->setText(d.getEnzyme());
	digestiontime_->setText(String(d.getDigestionTime()) );
  digestiontemperature_->setText(String(d.getTemperature()));
	digestionPH_->setText(String(d.getPh())); 
	
			
}

void DigestionVisualizer::store()
{
	try
	{		
		(*ptr_).setEnzyme(string((const char*)digestionenzyme_->text()));
				
		(*ptr_).setDigestionTime(String((const char*)digestiontime_->text()).toFloat() );
		
		(*ptr_).setTemperature(String((const char*)digestiontime_->text()).toFloat() );
		
		(*ptr_).setPh(String((const char*)digestiontime_->text()).toFloat() );
		
		tempdig_ = (*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new digestion data. "<<e.what()<<endl;
	}
}

void DigestionVisualizer::reject()
{
	try
	{
		load(tempdig_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original digestion data. "<<e.what()<<endl;
	} 
}
