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
	
	IonDetectorVisualizer::IonDetectorVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<IonDetector>()
	{
		addLabel_("Modify iondetector information.");	
		addSeparator_();  
		
		addComboBox_(iondetector_type_, "Type");
		addComboBox_(iondetector_ac_mode_, "Acquisition mode");
		addDoubleLineEdit_(iondetector_res_, "Resolution (in ns)" );
		addDoubleLineEdit_(iondetector_freq_, "ADC sampling frequency (in MHz)" );
		
		finishAdding_();			
	}
	
	void IonDetectorVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(iondetector_type_,& temp_.NamesOfType[temp_.getType()] , 1);
			fillComboBox_(iondetector_ac_mode_,& temp_.NamesOfAcquisitionMode[temp_.getAcquisitionMode()] , 1);
		}
		else
		{
			fillComboBox_(iondetector_type_, temp_.NamesOfType , IonDetector::SIZE_OF_TYPE);
			fillComboBox_(iondetector_ac_mode_, temp_.NamesOfAcquisitionMode , IonDetector::SIZE_OF_ACQUISITIONMODE);
			iondetector_type_->setCurrentIndex(temp_.getType()); 
			iondetector_ac_mode_->setCurrentIndex(temp_.getAcquisitionMode()); 
		}
		
		iondetector_res_->setText(String( temp_.getResolution() ).c_str());
		iondetector_freq_->setText(String( temp_.getADCSamplingFrequency() ).c_str());
	}
	
	void IonDetectorVisualizer::store()
	{
		String m(iondetector_res_->text().toStdString());
		ptr_->setResolution(m.toFloat());
		String n(iondetector_freq_->text().toStdString());
		ptr_->setADCSamplingFrequency(n.toFloat());
		ptr_->setType((IonDetector::Type)iondetector_type_->currentIndex());		
		ptr_->setAcquisitionMode((IonDetector::AcquisitionMode)iondetector_ac_mode_->currentIndex());		
		
		temp_=(*ptr_);
	}
	
	void IonDetectorVisualizer::undo_()
	{
		update_();
	}

}
