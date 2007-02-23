// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free ionsource; you can redistribute it and/or
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
// $Maintainer:  stefan_heess $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/IonSourceVisualizer.h>

//QT
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
IonSourceVisualizer::IonSourceVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="IonSource";
  
	addLabel("Modify ionsource information.");	
	addSeperator();  
	addComboBox(ionsource_inlet_type_, "Inlet type");
	addComboBox(ionsource_ionization_method_, "Ionization method");
	addComboBox(ionsource_polarity_, "Polarity");      
	
	finishAdding_();
}


void IonSourceVisualizer::load(IonSource &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempionsource_=s;
			
  fillComboBox(ionsource_inlet_type_, s.NamesOfInletType , IonSource::SIZE_OF_INLETTYPE);
	fillComboBox(ionsource_ionization_method_, s.NamesOfIonizationMethod , IonSource::SIZE_OF_IONIZATIONMETHOD);
	fillComboBox(ionsource_polarity_, s.NamesOfPolarity , IonSource::SIZE_OF_POLARITY);
	
	update_();
}

void IonSourceVisualizer::update_()
{
		ionsource_inlet_type_->setCurrentIndex(tempionsource_.getInletType()); 
		ionsource_ionization_method_->setCurrentIndex(tempionsource_.getIonizationMethod()); 
		ionsource_polarity_->setCurrentIndex(tempionsource_.getPolarity()); 
}

void IonSourceVisualizer::store()
{
	try
	{
	
		(*ptr_).setInletType((IonSource::InletType)ionsource_inlet_type_->currentIndex());		
		(*ptr_).setIonizationMethod((IonSource::IonizationMethod)ionsource_ionization_method_->currentIndex());		
		(*ptr_).setPolarity((IonSource::Polarity)ionsource_polarity_->currentIndex());		
		
		tempionsource_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new ion source data. "<<e.what()<<endl;
	}
	
}

void IonSourceVisualizer::reject()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original ion source data. "<<e.what()<<endl;
	}
	
}

}
