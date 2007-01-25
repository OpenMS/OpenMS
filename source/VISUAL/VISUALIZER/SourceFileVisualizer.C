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
#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>


//QT
#include <qwidget.h>
#include <qlabel.h> 
#include <qlineedit.h>
#include <qtextedit.h>
#include <qpushbutton.h>
#include <iostream>
#include <vector>


//using namespace std;
using namespace OpenMS;
using namespace std;

//Constructor
SourceFileVisualizer::SourceFileVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
        addLabel("Modify source file information");
        addSeperator();	
        addLineEdit(name_of_file_, "Name of file" );
	addLineEdit(path_to_file_, "Path to file" );
	addLineEdit(file_type_, "File type" );
	addVSpacer();		
	addSeperator();
        addLabel("Save changes or restore original data");	
        addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
        connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
		
}


void SourceFileVisualizer::load(SourceFile &s)
{
  ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempSourceFile_=s;
        name_of_file_->setText(s.getNameOfFile());
	path_to_file_->setText(s.getPathToFile() );
        file_type_->setText(String(s.getFileType()));
		
			
}

void SourceFileVisualizer::store()
{
	try
	{
				
		(*ptr_).setNameOfFile(string((const char*)name_of_file_->text()));
				
		(*ptr_).setPathToFile(string((const char*)path_to_file_->text()) );
		
		(*ptr_).setFileType(string((const char*)file_type_->text() ));
				
		tempSourceFile_=(*ptr_);
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new source file data. "<<e.what()<<endl;
	}
}

        void SourceFileVisualizer::reject()
{
	try
	{
		load(tempSourceFile_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original source file data. "<<e.what()<<endl;
	} 
}
