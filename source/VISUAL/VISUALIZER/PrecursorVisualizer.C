// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free precursor; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Precursor Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
PrecursorVisualizer::PrecursorVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
	type_="Precursor";
  
	addLabel("Modify processing method information.");	
	addSeperator();  
	
	addComboBox(precursor_activation_method_, "Activation method");
	addLineEdit(precursor_activation_energy_, "Activation energy");
	addComboBox(precursor_energy_units_, "Energy unit");
	addLineEdit(precursor_window_size_, "Window size");	
	
	finishAdding_();
	
	// A validator to check the input for the energy.
	QDoubleValidator *precursor_energy_vali_ = new QDoubleValidator(precursor_activation_energy_);
	precursor_activation_energy_->setValidator(precursor_energy_vali_);
	
	// A validator to check the input for window size
	QDoubleValidator *precursor_window_size_vali_ = new QDoubleValidator(precursor_window_size_);
	precursor_window_size_->setValidator(precursor_window_size_vali_);
}


void PrecursorVisualizer::load(Precursor &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempprecursor_=s;
	
		
	fillComboBox(precursor_activation_method_, Precursor::NamesOfActivationMethod , Precursor::SIZE_OF_ACTIVATIONMETHOD);
	fillComboBox(precursor_energy_units_, Precursor::NamesOfEnergyUnits , Precursor::SIZE_OF_ENERGYUNITS);
	
	
	update_();
}

void PrecursorVisualizer::update_()
{		
		
	precursor_activation_method_->setCurrentIndex(tempprecursor_.getActivationMethod()); 
	precursor_activation_energy_->setText(String( tempprecursor_.getActivationEnergy() ).c_str() );
	precursor_energy_units_->setCurrentIndex(tempprecursor_.getActivationEnergyUnit()); 
 	precursor_window_size_->setText(String( tempprecursor_.getWindowSize() ).c_str() );		
}

void PrecursorVisualizer::store()
{
	try
	{
		(*ptr_).setActivationMethod((Precursor::ActivationMethod)precursor_activation_method_->currentIndex());		
		(*ptr_).setActivationEnergy(precursor_activation_energy_->text().toFloat());		
		(*ptr_).setActivationEnergyUnit((Precursor::EnergyUnits)precursor_energy_units_->currentIndex());	
		(*ptr_).setWindowSize(precursor_window_size_->text().toFloat());		
		
		tempprecursor_=(*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new processing method data. "<<e.what()<<endl;
	}
	
}

void PrecursorVisualizer::reject()
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
