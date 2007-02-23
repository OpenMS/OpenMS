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

// OpenMS
#include <OpenMS/VISUAL/VISUALIZER/SourceFileVisualizer.h>

// QT
#include <QtGui/QLineEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
SourceFileVisualizer::SourceFileVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
	addLabel("Modify source file information");
	addSeperator();	
	addLineEdit(name_of_file_, "Name of file" );
	addLineEdit(path_to_file_, "Path to file" );
	addLineEdit(file_size_, "File size (in MB)" );
	addLineEdit(file_type_, "File type" );
	addLineEdit(sha1_, "SHA1 hash value" );
	
	finishAdding_();
		
}


void SourceFileVisualizer::load(SourceFile &s)
{
  ptr_ = &s;
	
	//Copy of current object for restoring the original values
	tempSourceFile_=s;
  name_of_file_->setText(s.getNameOfFile().c_str());
	path_to_file_->setText(s.getPathToFile().c_str() );
	file_size_->setText(String(s.getFileSize()).c_str());
  file_type_->setText(String(s.getFileType()).c_str());
	sha1_->setText(String(s.getSha1()).c_str());
		
			
}

void SourceFileVisualizer::store()
{
	try
	{
				
		(*ptr_).setNameOfFile(name_of_file_->text().toStdString());
		(*ptr_).setPathToFile(path_to_file_->text().toStdString());
		(*ptr_).setFileSize(file_size_->text().toFloat());
		(*ptr_).setFileType(file_type_->text().toStdString());
		(*ptr_).setSha1(sha1_->text().toStdString());
				
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

}
