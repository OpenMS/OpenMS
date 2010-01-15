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
// $Maintainer: $
// $Authors: Marc Sturm $
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
	
	IonDetectorVisualizer::IonDetectorVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<IonDetector>()
	{
		addLabel_("Modify iondetector information.");	
		addSeparator_();  
		
		addIntLineEdit_(order_, "Order" );
		addComboBox_(type_, "Type");
		addComboBox_(ac_mode_, "Acquisition mode");
		addDoubleLineEdit_(res_, "Resolution (in ns)" );
		addDoubleLineEdit_(freq_, "ADC sampling frequency (in Hz)" );
		
		finishAdding_();			
	}
	
	void IonDetectorVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(type_,& temp_.NamesOfType[temp_.getType()] , 1);
			fillComboBox_(ac_mode_,& temp_.NamesOfAcquisitionMode[temp_.getAcquisitionMode()] , 1);
		}
		else
		{
			fillComboBox_(type_, temp_.NamesOfType , IonDetector::SIZE_OF_TYPE);
			fillComboBox_(ac_mode_, temp_.NamesOfAcquisitionMode , IonDetector::SIZE_OF_ACQUISITIONMODE);
			type_->setCurrentIndex(temp_.getType()); 
			ac_mode_->setCurrentIndex(temp_.getAcquisitionMode()); 
		}

		order_->setText(String(temp_.getOrder()).c_str());
		res_->setText(String( temp_.getResolution() ).c_str());
		freq_->setText(String( temp_.getADCSamplingFrequency() ).c_str());
	}
	
	void IonDetectorVisualizer::store()
	{
		ptr_->setOrder(order_->text().toInt());
		ptr_->setResolution(res_->text().toDouble());
		ptr_->setADCSamplingFrequency(freq_->text().toDouble());
		ptr_->setType((IonDetector::Type)type_->currentIndex());		
		ptr_->setAcquisitionMode((IonDetector::AcquisitionMode)ac_mode_->currentIndex());		
		
		temp_=(*ptr_);
	}
	
	void IonDetectorVisualizer::undo_()
	{
		update_();
	}

}
