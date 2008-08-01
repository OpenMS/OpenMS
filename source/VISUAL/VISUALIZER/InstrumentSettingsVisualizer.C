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
//  License as published by the Free InstrumentSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
InstrumentSettingsVisualizer::InstrumentSettingsVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
	type_="InstrumentSettings";
  
	addLabel("Modify the settings of the instrument.");	
	addSeperator();  
	addComboBox(instrumentsettings_scan_mode_, "Scan mode");
	addComboBox(instrumentsettings_polarity_, "Polarity");
	addDoubleLineEdit(instrumentsettings_mz_range_start_, "Scan begin (in m/z dimension)");
	addDoubleLineEdit(instrumentsettings_mz_range_stop_, "Scan stop (in m/z dimension)");
		
	finishAdding_();
		
}



void InstrumentSettingsVisualizer::load(InstrumentSettings &is)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &is;
	
	//Copy of current object for restoring the original values
	tempinstrumentsettings_=is;
		
	update_();
}

void InstrumentSettingsVisualizer::update_()
{
		if(! isEditable())
		{
			fillComboBox(instrumentsettings_scan_mode_, &tempinstrumentsettings_.NamesOfScanMode[tempinstrumentsettings_.getScanMode()] , 1);
			fillComboBox(instrumentsettings_polarity_, &IonSource::NamesOfPolarity[tempinstrumentsettings_.getPolarity()] , 1);
		}
		else
		{
			fillComboBox(instrumentsettings_scan_mode_, InstrumentSettings::NamesOfScanMode , InstrumentSettings::SIZE_OF_SCANMODE);
			fillComboBox(instrumentsettings_polarity_, IonSource::NamesOfPolarity , IonSource::SIZE_OF_POLARITY);
			
			instrumentsettings_scan_mode_->setCurrentIndex(tempinstrumentsettings_.getScanMode()); 
			instrumentsettings_polarity_->setCurrentIndex(tempinstrumentsettings_.getPolarity()); 
		}
		
		
		instrumentsettings_mz_range_start_->setText(String(tempinstrumentsettings_.getMzRangeStart() ).c_str() );
		instrumentsettings_mz_range_stop_->setText(String(tempinstrumentsettings_.getMzRangeStop() ).c_str() );
}

void InstrumentSettingsVisualizer::store_()
{
	try
	{
			
		(*ptr_).setScanMode((InstrumentSettings::ScanMode)instrumentsettings_scan_mode_->currentIndex());		
		(*ptr_).setPolarity((IonSource::Polarity)instrumentsettings_polarity_->currentIndex());		
		(*ptr_).setMzRangeStart(instrumentsettings_mz_range_start_->text().toFloat());
		(*ptr_).setMzRangeStop(instrumentsettings_mz_range_stop_->text().toFloat());
		
		tempinstrumentsettings_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new InstrumentSettings data. "<<e.what()<<endl;
	}
	
}

void InstrumentSettingsVisualizer::reject_()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original InstrumentSettings data. "<<e.what()<<endl;
	}
	
}

}
