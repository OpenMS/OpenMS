// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free spectrumsettings; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free SpectrumSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

//QT
#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
SpectrumSettingsVisualizer::SpectrumSettingsVisualizer(bool editable, QWidget *parent) 
	: BaseVisualizer(editable, parent)
{
	type_="SpectrumSettings";
  
	addLabel("Modify the settings of the spectrum.");	
	addSeperator();  
	addComboBox(spectrumsettings_type_, "Type of spectrum");
	addTextEdit(spectrumsettings_comment_, "Comment");
		
	finishAdding_();
}


void SpectrumSettingsVisualizer::load(SpectrumSettings &s)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempspectrumsettings_=s;
			
  fillComboBox(spectrumsettings_type_, s.NamesOfSpectrumType , SpectrumSettings::SIZE_OF_SPECTRUMTYPE);
		
	update_();
}

void SpectrumSettingsVisualizer::update_()
{
		spectrumsettings_type_->setCurrentIndex(tempspectrumsettings_.getType()); 
		spectrumsettings_comment_->setText(tempspectrumsettings_.getComment().c_str());
	  
}

void SpectrumSettingsVisualizer::store()
{
	try
	{
			
		(*ptr_).setType((SpectrumSettings::SpectrumType)spectrumsettings_type_->currentIndex());			
		(*ptr_).setComment(spectrumsettings_comment_->toPlainText().toStdString());
		
		tempspectrumsettings_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new SpectrumSettings data. "<<e.what()<<endl;
	}
	
}

void SpectrumSettingsVisualizer::reject()
{
	
	try
	{

		update_();
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original SpectrumSettings data. "<<e.what()<<endl;
	}
	
}

}
