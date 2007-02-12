// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  this library is free identification; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Identification Foundation; either
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


#include <OpenMS/VISUAL/VISUALIZER/IdentificationVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/METADATA/Identification.h>
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
IdentificationVisualizer::IdentificationVisualizer(bool editable, QWidget *parent, MSMetaDataExplorer *caller, const char *name) : BaseVisualizer(editable,parent, name)
{
	type_="Identification";
	pidv_caller_= caller;
	
	addLabel("Modify identification information.");	
	addSeperator();        
	addLineEdit(identification_date_, "Date and Time of DB search" );
	addLineEdit(identification_threshold_, "Peptide significance threshold" );	
	addSeperator();       
	addLabel("Show peptide hits with score equal or higher than current threshold.");
	addLabel("(To show all peptide hits set threshold to 0).");
	addButton(updatebutton_, "Show peptide hits");
	addVSpacer();
	addSeperator();
	addLabel("Show peptide hits referencing a certain protein.");
	addLineEdit(identification_ref_date_, "Date and Time of DB search (YYYY-MM-DD hh:mm:ss)" );
	addLineEdit(identification_acc_, "Accession number of the protein." );
	addButton(updatebutton2_, "Show peptide hits");
	addVSpacer();
	addSeperator();
	addLabel("Show peptide hits NOT referencing any protein.");
	addButton(updatebutton3_, "Show peptide hits");
	
	finishAdding_();
	
	connect(updatebutton_, SIGNAL(clicked()), this, SLOT(updateTree()) );
	connect(updatebutton2_, SIGNAL(clicked()), this, SLOT(searchRefPeptides()) );
	connect(updatebutton3_, SIGNAL(clicked()), this, SLOT(searchNonRefPeptides()) );
		
	// A validator to check the input for the protein significance threshold
	QDoubleValidator *identification_threshold_vali_ = new QDoubleValidator(identification_threshold_);
	identification_threshold_->setValidator(identification_threshold_vali_);	
}


void IdentificationVisualizer::load(Identification &s, int tree_item_id)
{
  //Pointer to current object to keep track of the actual object
	ptr_ = &s;
	
	// id of the item in the tree
	tree_id_ = tree_item_id;
	
	//Copy of current object for restoring the original values
	tempidentification_=s;
  
  String str;
  tempidentification_.getDateTime().get(str);
	identification_date_->setText(str); 
	identification_threshold_->setText(String ( tempidentification_.getPeptideSignificanceThreshold() ) );					
}


void IdentificationVisualizer::updateTree()
{
	String m((const char*) identification_threshold_->text()) ;
	tempidentification_.setPeptideSignificanceThreshold(m.toFloat() );
			
	pidv_caller_->updatePeptideHits_(tempidentification_, tree_id_ );
	
}

void IdentificationVisualizer::searchRefPeptides()
{
	
	String ref_date((const char*) identification_ref_date_->text()) ;
	String ref_acc((const char*) identification_acc_->text()) ;
	
	pidv_caller_->updateRefPeptideHits_(tempidentification_, tree_id_, ref_date, ref_acc);
	
}

void IdentificationVisualizer::searchNonRefPeptides()
{
	pidv_caller_->updateNonRefPeptideHits_(tempidentification_ , tree_id_ );
}

void IdentificationVisualizer::store()
{
	try
	{
		String m((const char*) identification_threshold_->text()) ;
		(*ptr_).setProteinSignificanceThreshold(m.toFloat() );
		
		DateTime date;
		String o((const char*) identification_date_->text());
		date.set(o);
		(*ptr_).setDateTime(date);
				
		tempidentification_=(*ptr_);		
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new protein identification data. "<<e.what()<<endl;
	}
	
}

void IdentificationVisualizer::reject()
{
	
	try
	{
		//load(tempidentification_);
		load(*ptr_, tree_id_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original protein identification data. "<<e.what()<<endl;
	}
	
}

