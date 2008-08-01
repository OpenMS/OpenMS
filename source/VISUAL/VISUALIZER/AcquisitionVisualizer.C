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
//  License as published by the Free Software Foundation; either
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
// $Maintainer: Marc Sturm   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
AcquisitionVisualizer::AcquisitionVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
  
	addLabel("Show Acquisition information");		
	addSeperator();
	addIntLineEdit(acquisitionnumber_, "Index/Number of the scan" );
	acquisitionnumber_->setReadOnly(true);
		
	finishAdding_();
	
}


void AcquisitionVisualizer::load(Acquisition &a)
{
  ptr_ = &a;
	
	//Copy of current object for restoring the original values
	tempAcquisition_=a;
  acquisitionnumber_->setText(String(tempAcquisition_.getNumber()).c_str() );
				
}

void AcquisitionVisualizer::store_()
{
	try
	{
				
		//(*ptr_).setNumber(String((const char*)acquisitionnumber_->text()).toInt() );
					
		tempAcquisition_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new Acquisition data. "<<e.what()<<endl;
	}
}

void AcquisitionVisualizer::reject_()
{
	try
	{
		load(tempAcquisition_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original Acquisition data. "<<e.what()<<endl;
	} 
}

}
