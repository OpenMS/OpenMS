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
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

using namespace std;

namespace OpenMS
{

//Constructor
MetaInfoDescriptionVisualizer::MetaInfoDescriptionVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
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
	
  metainfodescription_name_->setText(tempMetaInfoDescription_.getName().c_str() );
	metainfodescription_comment_->setText(tempMetaInfoDescription_.getComment().c_str() );
				
}

void MetaInfoDescriptionVisualizer::store_()
{
	try
	{
				
		(*ptr_).setName(metainfodescription_name_->text().toStdString());
		(*ptr_).setComment(metainfodescription_comment_->toPlainText().toStdString());
					
		tempMetaInfoDescription_ = (*ptr_);
		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new MetaInfoDescription data. "<<e.what()<<endl;
	}
}

void MetaInfoDescriptionVisualizer::reject_()
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

}
