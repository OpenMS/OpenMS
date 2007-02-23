// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free processingmethod; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free ProcessingMethod Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/ProcessingMethodVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
ProcessingMethodVisualizer::ProcessingMethodVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="ProcessingMethod";
  
	addLabel("Modify processing method information.");	
	addSeperator();  
	
	addComboBox(processingmethod_deisotoping_, "Deisotoping");
	addComboBox(processingmethod_charge_deconvolution_, "Charge deconvolution");
	addComboBox(processingmethod_method_, "Method");
	addLineEdit(processingmethod_intensity_cutoff_, "Intensity cutoff");	
	
	finishAdding_();
	
	// A validator to check the input for intensity cutoff
	QDoubleValidator *processingmethod_intensity_cutoff_vali_ = new QDoubleValidator(processingmethod_intensity_cutoff_);
	processingmethod_intensity_cutoff_->setValidator(processingmethod_intensity_cutoff_vali_);
}


void ProcessingMethodVisualizer::load(ProcessingMethod &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempprocessingmethod_=s;
	
	//An array for bool values
	String bool_values_[3]= {"FALSE","TRUE"};		
	
	fillComboBox(processingmethod_method_, SpectrumSettings::NamesOfSpectrumType , SpectrumSettings::SIZE_OF_SPECTRUMTYPE);
	fillComboBox(processingmethod_deisotoping_, bool_values_ , 2);
	fillComboBox(processingmethod_charge_deconvolution_, bool_values_ , 2);
	
	update_();
}

void ProcessingMethodVisualizer::update_()
{		
		
		//update deisotoping
		if(tempprocessingmethod_.getDeisotoping())
		{
			processingmethod_deisotoping_->setCurrentIndex(1);
		}
		else 
		{ 
			processingmethod_deisotoping_->setCurrentIndex(0);
		}
		
		
		//update charge_deconvolution
		if(tempprocessingmethod_.getChargeDeconvolution())
		{
			processingmethod_charge_deconvolution_->setCurrentIndex(1);
		}
		else 
		{ 
			processingmethod_charge_deconvolution_->setCurrentIndex(0);
		}
		
		processingmethod_intensity_cutoff_->setText(String( tempprocessingmethod_.getIntensityCutoff() ).c_str() );		
		processingmethod_method_->setCurrentIndex(tempprocessingmethod_.getSpectrumType()); 
	
}

void ProcessingMethodVisualizer::store()
{
	try
	{
		(*ptr_).setSpectrumType((SpectrumSettings::SpectrumType)processingmethod_method_->currentIndex());		
		(*ptr_).setDeisotoping(processingmethod_deisotoping_->currentIndex());		
		(*ptr_).setChargeDeconvolution(processingmethod_charge_deconvolution_->currentIndex());
		(*ptr_).setIntensityCutoff(processingmethod_intensity_cutoff_->text().toFloat());
			
		tempprocessingmethod_=(*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new processing method data. "<<e.what()<<endl;
	}
	
}

void ProcessingMethodVisualizer::reject()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original processing method data. "<<e.what()<<endl;
	}
	
}

}
