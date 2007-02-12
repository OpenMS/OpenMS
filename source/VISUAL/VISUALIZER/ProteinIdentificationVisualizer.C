// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free proteinidentification; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free ProteinIdentification Foundation; either
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


#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
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
ProteinIdentificationVisualizer::ProteinIdentificationVisualizer(bool editable, QWidget *parent, MSMetaDataExplorer *caller, const char *name) : BaseVisualizer(editable, parent, name)
{
	type_="ProteinIdentification";
	pidv_caller_= caller;
	
	addLabel("Modify protein identification information.");	
	addSeperator();        
	addLineEdit(proteinidentification_date_, "Date of search" );
	addLineEdit(proteinidentification_threshold_, "Protein significance threshold" );	
	addSeperator();       
	addLabel("Show protein hits with score equal or higher than current threshold.");
	addButton(updatebutton_, "Show protein hits");
	
	finishAdding_();
	
	connect(updatebutton_, SIGNAL(clicked()), this, SLOT(updateTree()) );
		
	// A validator to check the input for the protein significance threshold
	QDoubleValidator *proteinidentification_threshold_vali_ = new QDoubleValidator(proteinidentification_threshold_);
	proteinidentification_threshold_->setValidator(proteinidentification_threshold_vali_);	
}


void ProteinIdentificationVisualizer::load(ProteinIdentification &s, int tree_item_id)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	// id of the item in the tree
	tree_id_ = tree_item_id;
	
	//Copy of current object for restoring the original values
	tempproteinidentification_=s;
  
  String str;
  tempproteinidentification_.getDateTime().get(str);
	proteinidentification_date_->setText(str); 
	proteinidentification_threshold_->setText(String ( tempproteinidentification_.getProteinSignificanceThreshold() ) );			
}


void ProteinIdentificationVisualizer::updateTree()
{
	String m((const char*) proteinidentification_threshold_->text()) ;
	tempproteinidentification_.setProteinSignificanceThreshold(m.toFloat() );
	
	pidv_caller_->updateProteinHits_(tempproteinidentification_, tree_id_);
	
}


void ProteinIdentificationVisualizer::store()
{
	try
	{
		String m((const char*) proteinidentification_threshold_->text()) ;
		(*ptr_).setProteinSignificanceThreshold(m.toFloat() );
		
		DateTime date;
		String n((const char*) proteinidentification_date_->text());
		date.set(n);
		(*ptr_).setDateTime(date);
				
		tempproteinidentification_=(*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new protein identification data. "<<e.what()<<endl;
	}
	
}

void ProteinIdentificationVisualizer::reject()
{
	
	try
	{
		//load(tempproteinidentification_);
		load(*ptr_, tree_id_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original protein identification data. "<<e.what()<<endl;
	}
	
}

