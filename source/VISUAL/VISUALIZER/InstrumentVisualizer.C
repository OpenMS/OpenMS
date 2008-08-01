// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
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
// $Maintainer:  Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>
#include <string>

using namespace std;

namespace OpenMS
{

//Constructor
InstrumentVisualizer::InstrumentVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="Instrument";
  
	addLabel("Modify instrument information.");	
	addSeperator();        
	addLineEdit(instrument_name_, "Name" );
	addLineEdit(instrument_vendor_, "Vendor" );	
	addLineEdit(instrument_model_, "Model");
	addTextEdit(instrument_customizations_, "Customizations" );
	
	finishAdding_();
}



void InstrumentVisualizer::load(Instrument &s)
{
        //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempinstrument_=s;
			
  instrument_name_->setText(s.getName().c_str());
	instrument_vendor_->setText(s.getVendor().c_str());
	instrument_model_->setText(s.getModel().c_str());
  instrument_customizations_->setText(s.getCustomizations().c_str()); 
				
}

void InstrumentVisualizer::store_()
{
	try
	{
		
		(*ptr_).setName(instrument_name_->text().toStdString());
		(*ptr_).setVendor(instrument_vendor_->text().toStdString());
		(*ptr_).setModel(instrument_model_->text().toStdString());
		(*ptr_).setCustomizations(instrument_customizations_->toPlainText().toStdString());
		
		tempinstrument_ = (*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new instrument data. "<<e.what()<<endl;
	}
	
}

void InstrumentVisualizer::reject_()
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

}
