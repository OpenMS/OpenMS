// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free software; you can redistribute it and/or
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
// $Maintainer: stefan_heess   $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionInfoVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>


//QT
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>
#include <qvalidator.h>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
AcquisitionInfoVisualizer::AcquisitionInfoVisualizer(bool editable, QWidget *parent, const char *name) 
	: BaseVisualizer(editable, parent, name)
{
  
	addLabel("Show AcquisitionInfo information");		
	addSeperator();
	addLineEdit(acquisitioninfo_method_, "Method of combination" );
	
	finishAdding_();
	
	// A validator to check the input for the temperature.
	QIntValidator *acquisitioninfo_method_vali_ = new QIntValidator(acquisitioninfo_method_);
	acquisitioninfo_method_->setValidator(acquisitioninfo_method_vali_);
	
}


void AcquisitionInfoVisualizer::load(AcquisitionInfo &a)
{
  ptr_ = &a;
	
	//Copy of current object for restoring the original values
	tempAcquisitionInfo_=a;
  acquisitioninfo_method_->setText(tempAcquisitionInfo_.getMethodOfCombination() );
				
}

void AcquisitionInfoVisualizer::store()
{
	try
	{
				
		(*ptr_).setMethodOfCombination(String((const char*)acquisitioninfo_method_->text()) );
					
		tempAcquisitionInfo_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new AcquisitionInfo data. "<<e.what()<<endl;
	}
}

void AcquisitionInfoVisualizer::reject()
{
	try
	{
		load(tempAcquisitionInfo_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original AcquisitionInfo data. "<<e.what()<<endl;
	} 
}
