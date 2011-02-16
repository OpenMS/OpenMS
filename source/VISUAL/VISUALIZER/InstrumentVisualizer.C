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
//  License as published by the Free Instrument Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/InstrumentVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>
#include <string>

using namespace std;

namespace OpenMS
{
	
	InstrumentVisualizer::InstrumentVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Instrument>()
	{
		addLabel_("Modify instrument information.");	
		addSeparator_();        
		addLineEdit_(name_, "Name" );
		addLineEdit_(vendor_, "Vendor" );	
		addLineEdit_(model_, "Model");
		addTextEdit_(customizations_, "Customizations" );
		addComboBox_(ion_optics_,"Ion optics");
		
		finishAdding_();
	}
	
	void InstrumentVisualizer::update_()
	{
	  name_->setText(temp_.getName().c_str());
		vendor_->setText(temp_.getVendor().c_str());
		model_->setText(temp_.getModel().c_str());
	  customizations_->setText(temp_.getCustomizations().c_str()); 

		if(! isEditable())
		{
			fillComboBox_(ion_optics_, &temp_.NamesOfIonOpticsType[temp_.getIonOptics()]  , 1);
		}
		else
		{
			fillComboBox_(ion_optics_, temp_.NamesOfIonOpticsType  , Instrument::SIZE_OF_IONOPTICSTYPE);
		}
	}
	
	void InstrumentVisualizer::store()
	{
		ptr_->setName(name_->text());
		ptr_->setVendor(vendor_->text());
		ptr_->setModel(model_->text());
		ptr_->setCustomizations(customizations_->toPlainText());
		ptr_->setIonOptics((Instrument::IonOpticsType)ion_optics_->currentIndex());		
		
		temp_ = (*ptr_);		
	}
	
	void InstrumentVisualizer::undo_()
	{
		update_();
	}

}
