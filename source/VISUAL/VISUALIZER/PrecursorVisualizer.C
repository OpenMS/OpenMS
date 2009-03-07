// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: $
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
		addDoubleLineEdit_(precursor_window_low_, "Isolation window lower bound");
		addDoubleLineEdit_(precursor_window_up_, "Isolation window upper bound");
		
		finishAdding_();
	}
	
	void PrecursorVisualizer::update_()
	{		
		if(! isEditable())
		{
			fillComboBox_(precursor_activation_method_,& temp_.NamesOfActivationMethod[temp_.getActivationMethod()] , 1);
		}
		else
		{
			fillComboBox_(precursor_activation_method_, Precursor::NamesOfActivationMethod , Precursor::SIZE_OF_ACTIVATIONMETHOD);
			precursor_activation_method_->setCurrentIndex(temp_.getActivationMethod()); 
		}
		
		precursor_activation_energy_->setText(String( temp_.getActivationEnergy() ).c_str() );
		precursor_window_low_->setText(String( temp_.getIsolationWindowLowerBound() ).c_str() );
		precursor_window_up_->setText(String( temp_.getIsolationWindowUpperBound() ).c_str() );
	}
	
	void PrecursorVisualizer::store()
	{
		ptr_->setActivationMethod((Precursor::ActivationMethod)precursor_activation_method_->currentIndex());		
		ptr_->setActivationEnergy(precursor_activation_energy_->text().toFloat());
		ptr_->setIsolationWindowLowerBound(precursor_window_low_->text().toFloat());		
		ptr_->setIsolationWindowUpperBound(precursor_window_up_->text().toFloat());		
		
		temp_=(*ptr_);
	}
	
	void PrecursorVisualizer::undo_()
	{
		update_();
	}

}
