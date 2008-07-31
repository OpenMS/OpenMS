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
// $Maintainer:  Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
SoftwareVisualizer::SoftwareVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	type_="Software";
  
	addLabel("Modify software information.");	
	addSeperator();        
	addLineEdit(software_name_, "Name" );
	addLineEdit(software_version_, "Version" );	
	addTextEdit(software_comment_, "Comment");
	addLineEdit(software_completion_time_, "Completion time" );

	finishAdding_();
	
	
}



void SoftwareVisualizer::load(Software &s)
{
        //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempsoftware_=s;
			
  software_name_->setText(s.getName().c_str());
	software_version_->setText(s.getVersion().c_str());
	software_comment_->setText(s.getComment().c_str());
  String str;
  s.getCompletionTime().get(str);
	software_completion_time_->setText(str.c_str()); 
				
}

void SoftwareVisualizer::store_()
{
	try
	{
		(*ptr_).setName(software_name_->text().toStdString());
		(*ptr_).setVersion(software_version_->text().toStdString());
		String m(software_completion_time_->text().toStdString());
		DateTime date;
		
		try
		{
			date.set(m);
			(*ptr_).setCompletionTime(date);
		}
		catch(exception& e)
		{
			if(date.isNull())
			{
				std::string status= "Format of date in SOFTWARE is not correct.";
				emit sendStatus(status);
			}
		}
		
		
		(*ptr_).setComment(software_comment_->toPlainText().toStdString());
		
		tempsoftware_=(*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new software data. "<<e.what()<<endl;
	}
	
}

void SoftwareVisualizer::reject_()
{
	
	try
	{
		load(tempsoftware_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original software data. "<<e.what()<<endl;
	}
	
}

}
