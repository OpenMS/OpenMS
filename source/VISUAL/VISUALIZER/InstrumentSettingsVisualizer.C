// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free InstrumentSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/InstrumentSettingsVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	InstrumentSettingsVisualizer::InstrumentSettingsVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<InstrumentSettings>()
	{
		addLabel_("Modify the settings of the instrument.");	
		addSeparator_();  
		addComboBox_(instrumentsettings_scan_mode_, "Scan mode");
		addBooleanComboBox_(zoom_scan_,"Zoom scan"); 
		addComboBox_(instrumentsettings_polarity_, "Polarity");
			
		finishAdding_();
	}

	void InstrumentSettingsVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(instrumentsettings_scan_mode_,& temp_.NamesOfScanMode[temp_.getScanMode()] , 1);
			fillComboBox_(instrumentsettings_polarity_,& IonSource::NamesOfPolarity[temp_.getPolarity()] , 1);

		}
		else
		{
			fillComboBox_(instrumentsettings_scan_mode_, InstrumentSettings::NamesOfScanMode , InstrumentSettings::SIZE_OF_SCANMODE);
			fillComboBox_(instrumentsettings_polarity_, IonSource::NamesOfPolarity , IonSource::SIZE_OF_POLARITY);
			

			instrumentsettings_scan_mode_->setCurrentIndex(temp_.getScanMode()); 
			zoom_scan_->setCurrentIndex(temp_.getZoomScan());
			instrumentsettings_polarity_->setCurrentIndex(temp_.getPolarity()); 
		}
	}
	
	void InstrumentSettingsVisualizer::store()
	{
		ptr_->setScanMode((InstrumentSettings::ScanMode)instrumentsettings_scan_mode_->currentIndex());
		ptr_->setZoomScan(zoom_scan_->currentIndex());		
		ptr_->setPolarity((IonSource::Polarity)instrumentsettings_polarity_->currentIndex());		
		
		temp_=(*ptr_);
	}
	
	void InstrumentSettingsVisualizer::undo_()
	{
		update_();
	}

}
