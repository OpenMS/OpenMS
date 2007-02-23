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
// $Maintainer: stefan_heess  $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>

//QT
#include <QtGui/QValidator>
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
DigestionVisualizer::DigestionVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
  type_="Digestion";
	
	addLabel("Modify Digestion information");		
	addSeperator();
	addLineEdit(treatmenttype_, "Treatment type" );
	addTextEdit(treatmentcomment_, "Comment" );
	addLineEdit(digestionenzyme_, "Enzyme" );
	addLineEdit(digestiontime_, "Digestion time (in minutes)" );
	addLineEdit(digestiontemperature_, "Temperature (in °C)" );
	addLineEdit(digestionPH_, "PH" );
	
	finishAdding_();
	
	// A validator to check the input for the time
	QDoubleValidator *timevali_ = new QDoubleValidator(digestiontime_);
	digestiontime_->setValidator(timevali_);
	// A validator to check the input for the temperature
	QDoubleValidator *tempvali_ = new QDoubleValidator(digestiontemperature_);
	digestiontemperature_->setValidator(tempvali_);
	// A validator to check the input for the ph value
	QDoubleValidator *phvali_ = new QDoubleValidator(digestionPH_);
	digestionPH_->setValidator(phvali_);
	
}


void DigestionVisualizer::load(Digestion &d)
{
  ptr_ = &d;
	
	//Copy of current object for restoring the original values
	tempdig_=d;
	treatmenttype_->setText(tempdig_.getType().c_str());
	treatmenttype_->setReadOnly(true);
	treatmentcomment_->setText(tempdig_.getComment().c_str());
  digestionenzyme_->setText(tempdig_.getEnzyme().c_str());
	digestiontime_->setText(String(tempdig_.getDigestionTime()).c_str() );
  digestiontemperature_->setText(String(tempdig_.getTemperature()).c_str());
	digestionPH_->setText(String(tempdig_.getPh()).c_str()); 
	
			
}

void DigestionVisualizer::store()
{
	try
	{		
		(*ptr_).setComment(treatmentcomment_->toPlainText().toStdString());
		(*ptr_).setEnzyme(digestionenzyme_->text().toStdString());
		(*ptr_).setDigestionTime(digestiontime_->text().toFloat());
		(*ptr_).setTemperature(digestiontime_->text().toFloat());
		(*ptr_).setPh(digestiontime_->text().toFloat());
		
		tempdig_ = (*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new digestion data. "<<e.what()<<endl;
	}
}

void DigestionVisualizer::reject()
{
	try
	{
		load(tempdig_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original digestion data. "<<e.what()<<endl;
	} 
}

}
