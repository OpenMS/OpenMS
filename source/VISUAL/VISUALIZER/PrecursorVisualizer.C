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
// $Maintainer:  Marc Sturm $
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
	
	PrecursorVisualizer::PrecursorVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Precursor>()
	{
		addLabel_("Modify processing method information.");	
		
		addSeparator_();  
		
		addComboBox_(precursor_activation_method_, "Activation method");
		addDoubleLineEdit_(precursor_activation_energy_, "Activation energy");
		addComboBox_(precursor_energy_units_, "Energy unit");
		addDoubleLineEdit_(precursor_window_size_, "Window size");	
		
		finishAdding_();
	}
	
	void PrecursorVisualizer::update_()
	{		
		if(! isEditable())
		{
			fillComboBox_(precursor_activation_method_,& temp_.NamesOfActivationMethod[temp_.getActivationMethod()] , 1);
			fillComboBox_(precursor_energy_units_,& temp_.NamesOfEnergyUnits[temp_.getActivationEnergyUnit()] , 1);
		}
		else
		{
			fillComboBox_(precursor_activation_method_, Precursor::NamesOfActivationMethod , Precursor::SIZE_OF_ACTIVATIONMETHOD);
			fillComboBox_(precursor_energy_units_, Precursor::NamesOfEnergyUnits , Precursor::SIZE_OF_ENERGYUNITS);
			precursor_activation_method_->setCurrentIndex(temp_.getActivationMethod()); 
			precursor_energy_units_->setCurrentIndex(temp_.getActivationEnergyUnit()); 
		}
		
		precursor_activation_energy_->setText(String( temp_.getActivationEnergy() ).c_str() );
		precursor_window_size_->setText(String( temp_.getWindowSize() ).c_str() );		
	}
	
	void PrecursorVisualizer::store()
	{
		ptr_->setActivationMethod((Precursor::ActivationMethod)precursor_activation_method_->currentIndex());		
		ptr_->setActivationEnergy(precursor_activation_energy_->text().toFloat());		
		ptr_->setActivationEnergyUnit((Precursor::EnergyUnits)precursor_energy_units_->currentIndex());	
		ptr_->setWindowSize(precursor_window_size_->text().toFloat());		
		
		temp_=(*ptr_);
	}
	
	void PrecursorVisualizer::undo_()
	{
		update_();
	}

}
