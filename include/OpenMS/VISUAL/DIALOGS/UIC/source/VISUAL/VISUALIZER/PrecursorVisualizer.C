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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/PrecursorVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>
#include <QtGui/QListWidget>

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

		addDoubleLineEdit_(mz_, "m/z");
		addDoubleLineEdit_(int_, "intensity");
		addIntLineEdit_(charge_, "charge");

		addDoubleLineEdit_(window_low_, "Lower offset from target m/z");
		addDoubleLineEdit_(window_up_, "Upper offset from target m/z");
			
		addListView_(activation_methods_, "Activation methods");
		addDoubleLineEdit_(activation_energy_, "Activation energy");

		finishAdding_();
	}
	
	void PrecursorVisualizer::update_()
	{
		mz_->setText(String( temp_.getMZ() ).c_str() );
		int_->setText(String( temp_.getIntensity() ).c_str() );
		charge_->setText(String( temp_.getCharge() ).c_str() );
		
		window_low_->setText(String( temp_.getIsolationWindowLowerOffset() ).c_str() );
		window_up_->setText(String( temp_.getIsolationWindowUpperOffset() ).c_str() );
		
		//actions
		activation_methods_->clear();
		for (Size i=0; i<Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
		{
			QListWidgetItem* item = new QListWidgetItem(activation_methods_);
			item->setText(QString::fromStdString(Precursor::NamesOfActivationMethod[i]));
			if (temp_.getActivationMethods().count(Precursor::ActivationMethod(i))==1)
			{
				item->setCheckState(Qt::Checked);
			}
			else
			{
				item->setCheckState(Qt::Unchecked);
			}
			if (isEditable())
			{
				item->setFlags(Qt::ItemIsEnabled | Qt::ItemIsUserCheckable);
			}
			else
			{
				item->setFlags(Qt::ItemIsEnabled);
			}
			activation_methods_->addItem(item);
		}
		
		activation_energy_->setText(String( temp_.getActivationEnergy() ).c_str() );
	}
	
	void PrecursorVisualizer::store()
	{
		ptr_->setMZ(mz_->text().toFloat());
		ptr_->setIntensity(int_->text().toFloat());
		ptr_->setCharge(charge_->text().toInt());
		
		ptr_->setIsolationWindowLowerOffset(window_low_->text().toFloat());		
		ptr_->setIsolationWindowUpperOffset(window_up_->text().toFloat());		

		ptr_->getActivationMethods().clear();
		for (UInt i=0; i<Precursor::SIZE_OF_ACTIVATIONMETHOD; ++i)
		{
			if (activation_methods_->item(i)->checkState()==Qt::Checked)
			{
				ptr_->getActivationMethods().insert(Precursor::ActivationMethod(i));
			}
		}
		ptr_->setActivationEnergy(activation_energy_->text().toFloat());
		
		temp_=(*ptr_);
	}
	
	void PrecursorVisualizer::undo_()
	{
		update_();
	}

}
