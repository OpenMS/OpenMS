// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free instrument; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Instrument Foundation; either
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


#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Instrument.h>



//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qstring.h>

//STL
#include <iostream>
#include <vector>
#include <string>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
InstrumentVisualizer::InstrumentVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="Instrument";
  
	addLabel("Modify instrument information.");	
	addSeperator();        
	addLineEdit(instrument_name_, "Name" );
	addLineEdit(instrument_vendor_, "Vendor" );	
	addLineEdit(instrument_model_, "Model");
	addTextEdit(instrument_customizations_, "Customizations" );
	
	addEmptyLine();
	addSeperator();
	addLabel("Save changes or restore original data.");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
			
}


void InstrumentVisualizer::load(Instrument &s)
{
        //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempinstrument_=s;
			
  instrument_name_->setText(s.getName());
	instrument_vendor_->setText(s.getVendor());
	instrument_model_->setText(s.getModel());
  instrument_customizations_->setText(s.getCustomizations()); 
				
}

void InstrumentVisualizer::store()
{
	try
	{
		
		(*ptr_).setName(string((const char*) instrument_name_->text()) );
		(*ptr_).setVendor(string((const char*) instrument_vendor_->text()) );
		(*ptr_).setModel(string((const char*) instrument_model_->text()));
		(*ptr_).setCustomizations(string((const char*) instrument_customizations_->text()) );
		
		tempinstrument_ = (*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new instrument data. "<<e.what()<<endl;
	}
	
}

void InstrumentVisualizer::reject()
{
	
	try
	{
		load(tempinstrument_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original instrument data. "<<e.what()<<endl;
	}
	
}

