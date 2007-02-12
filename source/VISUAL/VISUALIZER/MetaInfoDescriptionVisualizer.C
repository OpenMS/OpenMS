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
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>
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
MetaInfoDescriptionVisualizer::MetaInfoDescriptionVisualizer(bool editable, QWidget *parent, const char *name) : BaseVisualizer(editable, parent, name)
{
  
	addLabel("Modify MetaInfoDescription information");		
	addSeperator();
	addLineEdit(metainfodescription_name_, "Name of peak annotations" );
	addTextEdit(metainfodescription_comment_, "Comment" );
		
	finishAdding_();
}


void MetaInfoDescriptionVisualizer::load(MetaInfoDescription &a)
{
  ptr_ = &a;
	
	//Copy of current object for restoring the original values
	tempMetaInfoDescription_=a;
	
  metainfodescription_name_->setText(tempMetaInfoDescription_.getName() );
	metainfodescription_comment_->setText(tempMetaInfoDescription_.getComment() );
				
}

void MetaInfoDescriptionVisualizer::store()
{
	try
	{
				
		(*ptr_).setName(String((const char*)metainfodescription_name_->text()) );
		(*ptr_).setComment(String((const char*)metainfodescription_comment_->text()) );
					
		tempMetaInfoDescription_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new MetaInfoDescription data. "<<e.what()<<endl;
	}
}

void MetaInfoDescriptionVisualizer::reject()
{
	try
	{
		load(tempMetaInfoDescription_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original MetaInfoDescription data. "<<e.what()<<endl;
	} 
}
