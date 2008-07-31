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
//  License as published by the Free ExperimentalSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/ExperimentalSettingsVisualizer.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
ExperimentalSettingsVisualizer::ExperimentalSettingsVisualizer(bool editable, QWidget *parent) 
: BaseVisualizer(editable, parent)
{
	type_="ExperimentalSettings";
  
	addLabel("Modify the settings of the experiment.");	
	addSeperator();  
	addComboBox(experimentalsettings_type_, "Type of the experiment");
	addLineEdit(experimentalsettings_date_, "Date of experiment");
	addTextEdit(experimentalsettings_comment_, "Comment");
	
	finishAdding_();
	
	
}



void ExperimentalSettingsVisualizer::load(ExperimentalSettings &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempexperimentalsettings_=s;
			
  
		
	update_();
}

void ExperimentalSettingsVisualizer::update_()
{
		if(! isEditable())
		{
			
			fillComboBox(experimentalsettings_type_, &tempexperimentalsettings_.NamesOfExperimentType[tempexperimentalsettings_.getType()] ,1); 
			
		}
		else
		{
			fillComboBox(experimentalsettings_type_, tempexperimentalsettings_.NamesOfExperimentType , ExperimentalSettings::SIZE_OF_EXPERIMENTTYPE);
			experimentalsettings_type_->setCurrentIndex(tempexperimentalsettings_.getType()); 
		}
		String str;
    tempexperimentalsettings_.getDate().get(str);
	  experimentalsettings_date_->setText(str.c_str()); 
		experimentalsettings_comment_->setText(tempexperimentalsettings_.getComment().c_str());
}

void ExperimentalSettingsVisualizer::store_()
{
	try
	{
		
		(*ptr_).setType((ExperimentalSettings::ExperimentType)experimentalsettings_type_->currentIndex());		
		Date date;
		String n(experimentalsettings_date_->text().toStdString());
		try
		{
			date.set(n);
			(*ptr_).setDate(date);
		}
		catch(exception& e)
		{
			if(date.isNull())
			{
				std::string status= "Format of date in EXPERIMENTALSETTINGS is not correct.";
				emit sendStatus(status);
			}
		}
		
		(*ptr_).setComment(experimentalsettings_comment_->toPlainText().toStdString());
		
		tempexperimentalsettings_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new ExperimentalSettings data. "<<e.what()<<endl;
	}
	
}

void ExperimentalSettingsVisualizer::reject_()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original ExperimentalSettings data. "<<e.what()<<endl;
	}
	
}

}
