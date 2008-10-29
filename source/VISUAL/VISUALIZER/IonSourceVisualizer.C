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
//  License as published by the Free IonSource Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>

//QT
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	IonSourceVisualizer::IonSourceVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<IonSource>()
	{
		addLabel_("Modify ionsource information.");	
		addSeparator_();  
		addComboBox_(ionsource_inlet_type_, "Inlet type");
		addComboBox_(ionsource_ionization_method_, "Ionization method");
		addComboBox_(ionsource_polarity_, "Polarity");      
		
		finishAdding_();
	}
	
	void IonSourceVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(ionsource_inlet_type_,& temp_.NamesOfInletType[temp_.getInletType()]  , 1);
			fillComboBox_(ionsource_ionization_method_,& temp_.NamesOfIonizationMethod[temp_.getIonizationMethod()] , 1);
			fillComboBox_(ionsource_polarity_,& temp_.NamesOfPolarity[temp_.getPolarity()] , 1);	
		}
		else
		{
			fillComboBox_(ionsource_inlet_type_, temp_.NamesOfInletType  , IonSource::SIZE_OF_INLETTYPE);
			fillComboBox_(ionsource_ionization_method_, temp_.NamesOfIonizationMethod , IonSource::SIZE_OF_IONIZATIONMETHOD);
			fillComboBox_(ionsource_polarity_, temp_.NamesOfPolarity , IonSource::SIZE_OF_POLARITY);
			
			ionsource_inlet_type_->setCurrentIndex(temp_.getInletType()); 
			ionsource_ionization_method_->setCurrentIndex(temp_.getIonizationMethod()); 
			ionsource_polarity_->setCurrentIndex(temp_.getPolarity()); 
		}
	}
	
	void IonSourceVisualizer::store()
	{
		ptr_->setInletType((IonSource::InletType)ionsource_inlet_type_->currentIndex());		
		ptr_->setIonizationMethod((IonSource::IonizationMethod)ionsource_ionization_method_->currentIndex());		
		ptr_->setPolarity((IonSource::Polarity)ionsource_polarity_->currentIndex());		
		
		temp_=(*ptr_);
	}
	
	void IonSourceVisualizer::undo_()
	{
		update_();
	}

}
