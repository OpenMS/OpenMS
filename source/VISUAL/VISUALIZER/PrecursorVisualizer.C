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
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Precursor.h>



//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qstring.h>
#include <qvalidator.h>

//STL
#include <iostream>
#include <vector>
#include <string>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
PrecursorVisualizer::PrecursorVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="Precursor";
  
	addLabel("Modify processing method information.");	
	addSeperator();  
	
	addComboBox(precursor_activation_method_, "Activation method");
	addLineEdit(precursor_activation_energy_, "Activation energy");
	addComboBox(precursor_energy_units_, "Energy unit");
		
	addVSpacer();
	addSeperator();
	addLabel("Save changes or restore original data.");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
	// A validator to check the input for the energy.
	QDoubleValidator *precursor_energy_vali_ = new QDoubleValidator(precursor_activation_energy_);
	precursor_activation_energy_->setValidator(precursor_energy_vali_);
}


void PrecursorVisualizer::load(Precursor &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempprecursor_=s;
	
		
	fillComboBox(precursor_activation_method_, Precursor::NamesOfActivationMethod , Precursor::SIZE_OF_ACTIVATIONMETHOD);
	fillComboBox(precursor_energy_units_, Precursor::NamesOfEnergyUnits , Precursor::SIZE_OF_ENERGYUNITS);
	
	
	update();
}

void PrecursorVisualizer::update()
{		
		
	precursor_activation_method_->setCurrentItem(tempprecursor_.getActivationMethod()); 
	precursor_activation_energy_->setText(String( tempprecursor_.getActivationEnergy() ) );
	precursor_energy_units_->setCurrentItem(tempprecursor_.getActivationEnergyUnit()); 
		
}

void PrecursorVisualizer::store()
{
	try
	{
		(*ptr_).setActivationMethod((Precursor::ActivationMethod)precursor_activation_method_->currentItem());		
		(*ptr_).setActivationEnergy(String((const char*)precursor_activation_energy_->text()  ).toFloat() );		
		(*ptr_).setActivationEnergyUnit((Precursor::EnergyUnits)precursor_energy_units_->currentItem());	
		
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

		update();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original processing method data. "<<e.what()<<endl;
	}
	
}

