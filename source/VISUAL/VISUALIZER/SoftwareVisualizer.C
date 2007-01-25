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
// $Maintainer:  stefan_heess $
// --------------------------------------------------------------------------s


#include <OpenMS/VISUAL/VISUALIZER/SoftwareVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Software.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>


//QT
#include <qlayout.h>
#include <qwidget.h>
#include <qcombobox.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qlistview.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <qstring.h>
#include <qvalidator.h>

#include <iostream>
#include <vector>
#include <string>

//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
SoftwareVisualizer::SoftwareVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
	type_="Software";
  
	addLabel("Modify software information.");	
	addSeperator();        
	addLineEdit(software_name_, "Name" );
	addLineEdit(software_version_, "Version" );	
	addTextEdit(software_comment_, "Comment");
	addLineEdit(software_completion_time_, "Completion time" );

	addEmptyLine();
	addSeperator();
	addLabel("Save changes or restore original data.");
	addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
			
}


void SoftwareVisualizer::load(Software &s)
{
        //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempsoftware_=s;
			
  software_name_->setText(s.getName());
	software_version_->setText(s.getVersion());
	software_comment_->setText(s.getComment());
  String str;
  s.getCompletionTime().get(str);
	software_completion_time_->setText(str); 
				
}

void SoftwareVisualizer::store()
{
	try
	{
		(*ptr_).setName(string((const char*) software_name_->text()) );
		(*ptr_).setVersion(string((const char*) software_version_->text()) );
		String m((const char*) software_completion_time_->text()) ;
		DateTime date;
		date.set(m);
		(*ptr_).setCompletionTime(date);
		(*ptr_).setComment(string((const char*) software_comment_->text()) );
		
		tempsoftware_=(*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new software data. "<<e.what()<<endl;
	}
	
}

void SoftwareVisualizer::reject()
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

