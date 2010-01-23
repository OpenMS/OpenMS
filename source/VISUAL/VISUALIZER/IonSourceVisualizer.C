// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>

//QT
#include <QtGui/QComboBox>
#include <QtGui/QLineEdit>

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
		
		addIntLineEdit_(order_, "Order" );
		addComboBox_(inlet_type_, "Inlet type");
		addComboBox_(ionization_method_, "Ionization method");
		addComboBox_(polarity_, "Polarity");      
		
		finishAdding_();
	}
	
	void IonSourceVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(inlet_type_,& temp_.NamesOfInletType[temp_.getInletType()]  , 1);
			fillComboBox_(ionization_method_,& temp_.NamesOfIonizationMethod[temp_.getIonizationMethod()] , 1);
			fillComboBox_(polarity_,& temp_.NamesOfPolarity[temp_.getPolarity()] , 1);	
		}
		else
		{
			fillComboBox_(inlet_type_, temp_.NamesOfInletType  , IonSource::SIZE_OF_INLETTYPE);
			fillComboBox_(ionization_method_, temp_.NamesOfIonizationMethod , IonSource::SIZE_OF_IONIZATIONMETHOD);
			fillComboBox_(polarity_, temp_.NamesOfPolarity , IonSource::SIZE_OF_POLARITY);
			
			inlet_type_->setCurrentIndex(temp_.getInletType()); 
			ionization_method_->setCurrentIndex(temp_.getIonizationMethod()); 
			polarity_->setCurrentIndex(temp_.getPolarity()); 
		}

		order_->setText(String(temp_.getOrder()).c_str());
	}
	
	void IonSourceVisualizer::store()
	{
		ptr_->setOrder(order_->text().toInt());
		ptr_->setInletType((IonSource::InletType)inlet_type_->currentIndex());		
		ptr_->setIonizationMethod((IonSource::IonizationMethod)ionization_method_->currentIndex());		
		ptr_->setPolarity((IonSource::Polarity)polarity_->currentIndex());		
		
		temp_=(*ptr_);
	}
	
	void IonSourceVisualizer::undo_()
	{
		update_();
	}

}
