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
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/IonSource.h>



//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qpushbutton.h>
#include <qstring.h>

//STL
#include <iostream>
#include <vector>
#include <string>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
IonSourceVisualizer::IonSourceVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="IonSource";
  
	addLabel("Modify ionsource information.");	
	addSeperator();  
	addComboBox(ionsource_inlet_type_, "Inlet type");
	addComboBox(ionsource_ionization_method_, "Ionization method");
	addComboBox(ionsource_polarity_, "Polarity");      
	
	addVSpacer();
	addSeperator();
	addLabel("Save changes or restore original data.");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
			
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
	
	update();
}

void IonSourceVisualizer::update()
{
		ionsource_inlet_type_->setCurrentItem(tempionsource_.getInletType()); 
		ionsource_ionization_method_->setCurrentItem(tempionsource_.getIonizationMethod()); 
		ionsource_polarity_->setCurrentItem(tempionsource_.getPolarity()); 
}

void IonSourceVisualizer::store()
{
	try
	{
	
		(*ptr_).setInletType((IonSource::InletType)ionsource_inlet_type_->currentItem());		
		(*ptr_).setIonizationMethod((IonSource::IonizationMethod)ionsource_ionization_method_->currentItem());		
		(*ptr_).setPolarity((IonSource::Polarity)ionsource_polarity_->currentItem());		
		
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

		update();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original ion source data. "<<e.what()<<endl;
	}
	
}

