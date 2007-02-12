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
// $Maintainer: stefan_heess   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
//#include <OpenMS/VISUAL/DataTable.h>

//QT
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>
#include <qvalidator.h>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
HPLCVisualizer::HPLCVisualizer(bool editable, QWidget *parent, const char *name) 
	: BaseVisualizer(editable, parent, name)
{
  
	addLabel("Modify HPLC information");		
	addSeperator();
	addLineEdit(hplcinstrument_, "Instrument" );
	addLineEdit(hplccolumn_, "Column" );
	addLineEdit(hplctemperature_, "Temperature (in °C)" );
	addLineEdit(hplcpressure_, "Pressure (in bar)" );
	addLineEdit(hplcflux_, "Flux (in µl/sec)" );
	addTextEdit(hplccomment_, "Comment");
		
	finishAdding_();
	
	// A validator to check the input for the temperature.
	QIntValidator *hplc_temperature_vali_ = new QIntValidator(hplctemperature_);
	hplctemperature_->setValidator(hplc_temperature_vali_);
	// A validator to check the input for the pressure.
	QIntValidator *hplc_pressure_vali_ = new QIntValidator(hplcpressure_);
	hplcpressure_->setValidator(hplc_pressure_vali_); 
	// A validator to check the input for the flux.
	QIntValidator *hplc_flux_vali_ = new QIntValidator(hplcflux_);
	hplcflux_->setValidator(hplc_flux_vali_);
	
}


void HPLCVisualizer::load(HPLC &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempHPLC_=h;
  hplcinstrument_->setText(h.getInstrument());
	hplccolumn_->setText(h.getColumn() );
  hplctemperature_->setText(String(h.getTemperature()));
	hplcpressure_->setText(String(h.getPressure()));
	hplcflux_->setText(String(h.getFlux()));
	hplccomment_->setText(h.getComment()); 
	
			
}

void HPLCVisualizer::store()
{
	try
	{
				
		
		(*ptr_).setInstrument(string((const char*)hplcinstrument_->text()));
				
		(*ptr_).setColumn(string((const char*)hplccolumn_->text()) );
		
		(*ptr_).setTemperature(String((const char*)hplctemperature_->text()).toInt() );
		
		(*ptr_).setPressure(String((const char*)hplcpressure_->text()).toInt() );
		
		(*ptr_).setFlux(String((const char*)hplcflux_->text()).toInt() );
		
		(*ptr_).setComment(string((const char*)hplccomment_->text()) );
		
		tempHPLC_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new HPLC data. "<<e.what()<<endl;
	}
}

void HPLCVisualizer::reject()
{
	try
	{
		load(tempHPLC_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original HPLC data. "<<e.what()<<endl;
	} 
}
