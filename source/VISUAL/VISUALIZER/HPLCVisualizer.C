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
// $Maintainer: Marc Sturm   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QValidator>
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
HPLCVisualizer::HPLCVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
  
	addLabel("Modify HPLC information");		
	addSeperator();
	addLineEdit(hplcinstrument_, "Instrument" );
	addLineEdit(hplccolumn_, "Column" );
	addIntLineEdit(hplctemperature_, "Temperature (in °C)" );
	addIntLineEdit(hplcpressure_, "Pressure (in bar)" );
	addIntLineEdit(hplcflux_, "Flux (in µl/sec)" );
	addTextEdit(hplccomment_, "Comment");
		
	finishAdding_();
	
	
	
}


void HPLCVisualizer::load(HPLC &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempHPLC_=h;
  hplcinstrument_->setText(h.getInstrument().c_str());
	hplccolumn_->setText(h.getColumn().c_str() );
  hplctemperature_->setText(String(h.getTemperature()).c_str());
	hplcpressure_->setText(String(h.getPressure()).c_str());
	hplcflux_->setText(String(h.getFlux()).c_str());
	hplccomment_->setText(h.getComment().c_str()); 
	
			
}

void HPLCVisualizer::store_()
{
	try
	{
				
		
		(*ptr_).setInstrument(hplcinstrument_->text().toStdString());
				
		(*ptr_).setColumn(hplccolumn_->text().toStdString());
		
		(*ptr_).setTemperature(hplctemperature_->text().toInt() );
		
		(*ptr_).setPressure(hplcpressure_->text().toInt() );
		
		(*ptr_).setFlux(hplcflux_->text().toInt());
		
		(*ptr_).setComment(hplccomment_->toPlainText().toStdString());
		
		tempHPLC_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new HPLC data. "<<e.what()<<endl;
	}
}

void HPLCVisualizer::reject_()
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

}
