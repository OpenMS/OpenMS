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
//  License as published by the Free IonDetector Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/IonDetectorVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>


//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
IonDetectorVisualizer::IonDetectorVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="IonDetector";
  
	addLabel("Modify iondetector information.");	
	addSeperator();  
	
	addComboBox(iondetector_type_, "Type");
	addComboBox(iondetector_ac_mode_, "Acquisition mode");
	addDoubleLineEdit(iondetector_res_, "Resolution (in ns)" );
	addDoubleLineEdit(iondetector_freq_, "ADC sampling frequency (in MHz)" );
	
	finishAdding_();			
}



void IonDetectorVisualizer::load(IonDetector &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempiondetector_=s;
			
	update_();
}

void IonDetectorVisualizer::update_()
{
		if(! isEditable())
		{
			fillComboBox(iondetector_type_, &tempiondetector_.NamesOfType[tempiondetector_.getType()] , 1);
			fillComboBox(iondetector_ac_mode_, &tempiondetector_.NamesOfAcquisitionMode[tempiondetector_.getAcquisitionMode()] , 1);
		}
		else
		{
			fillComboBox(iondetector_type_, tempiondetector_.NamesOfType , IonDetector::SIZE_OF_TYPE);
			fillComboBox(iondetector_ac_mode_, tempiondetector_.NamesOfAcquisitionMode , IonDetector::SIZE_OF_ACQUISITIONMODE);
			iondetector_type_->setCurrentIndex(tempiondetector_.getType()); 
			iondetector_ac_mode_->setCurrentIndex(tempiondetector_.getAcquisitionMode()); 
		}
		
		iondetector_res_->setText(String( tempiondetector_.getResolution() ).c_str());
		iondetector_freq_->setText(String( tempiondetector_.getADCSamplingFrequency() ).c_str());
		
		
}

void IonDetectorVisualizer::store_()
{
	try
	{

		String m(iondetector_res_->text().toStdString());
		(*ptr_).setResolution(m.toFloat());
		String n(iondetector_freq_->text().toStdString());
		(*ptr_).setADCSamplingFrequency(n.toFloat());
		(*ptr_).setType((IonDetector::Type)iondetector_type_->currentIndex());		
		(*ptr_).setAcquisitionMode((IonDetector::AcquisitionMode)iondetector_ac_mode_->currentIndex());		
		
		
		tempiondetector_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new ion detector data. "<<e.what()<<endl;
	}
	
}

void IonDetectorVisualizer::reject_()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original ion detector data. "<<e.what()<<endl;
	}
	
}

}
