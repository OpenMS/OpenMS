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
//  License as published by the Free PeptideIdentification Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/PeptideIdentificationVisualizer.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>

#include <iostream>

using namespace std;

namespace OpenMS
{

	PeptideIdentificationVisualizer::PeptideIdentificationVisualizer(bool editable, QWidget *parent, MSMetaDataExplorer *caller) : BaseVisualizer(editable,parent)
	{
		type_="PeptideIdentification";
		pidv_caller_= caller;
		
		addLabel("Modify ProteinIdentification information.");	
		addSeperator();        
//		addLineEdit(identification_date_, "Date and Time of DB search" );
		addDoubleLineEdit(identification_threshold_, "Peptide significance threshold" );	
		addSeperator();       
		addLabel("Show peptide hits with score equal or higher than current threshold.");
		addLabel("(To show all peptide hits set threshold to 0).");
		addButton(updatebutton_, "Show peptide hits");
		addSeperator();
		addLabel("Show peptide hits referencing a certain protein.");
//		addLineEdit(identification_ref_date_, "Date and Time of DB search (YYYY-MM-DD hh:mm:ss)" );
//		addLineEdit(identification_acc_, "Accession number of the protein." );
		addButton(updatebutton2_, "Show peptide hits");
		addSeperator();
		addLabel("Show peptide hits NOT referencing any protein.");
		addButton(updatebutton3_, "Show peptide hits");
		
		finishAdding_();
		
		connect(updatebutton_, SIGNAL(clicked()), this, SLOT(updateTree_()) );
		connect(updatebutton2_, SIGNAL(clicked()), this, SLOT(searchRefPeptides_()) );
		connect(updatebutton3_, SIGNAL(clicked()), this, SLOT(searchNonRefPeptides_()) );
	}

	void PeptideIdentificationVisualizer::load(PeptideIdentification &s, int tree_item_id)
	{
	  //Pointer to current object to keep track of the actual object
		ptr_ = &s;
		
		// id of the item in the tree
		tree_id_ = tree_item_id;
		
		//Copy of current object for restoring the original values
		tempidentification_=s;
	  
		identification_threshold_->setText(String ( tempidentification_.getSignificanceThreshold() ).c_str() );					
	}
	
	void PeptideIdentificationVisualizer::updateTree_()
	{
		String m(identification_threshold_->text().toStdString()) ;
		tempidentification_.setSignificanceThreshold(m.toFloat());
				
//	pidv_caller_->updatePeptideHits_(tempidentification_, tree_id_ );
		
	}
	
	void PeptideIdentificationVisualizer::searchRefPeptides_()
	{
		
//		String ref_date(identification_ref_date_->text().toStdString());
//		String ref_acc(identification_acc_->text().toStdString());
		
//		pidv_caller_->updateRefPeptideHits_(tempidentification_, tree_id_, ref_date, ref_acc);
		
	}
	
	void PeptideIdentificationVisualizer::searchNonRefPeptides_()
	{
		//pidv_caller_->updateNonRefPeptideHits_(tempidentification_ , tree_id_ );
	}
	
	void PeptideIdentificationVisualizer::store_()
	{
		
		
		try
		{
			String m(identification_threshold_->text().toStdString());
			(*ptr_).setSignificanceThreshold(m.toFloat());
			DateTime date;
//			String o(identification_date_->text().toStdString());
			String o("");
			try
			{
				date.set(o);
				//(*ptr_).setDateTime(date);
			}
			catch(exception& e)
			{
				if(date.isNull())
				{
					std::string status= "Format of date in ProteinIdentification is not correct.";
					emit sendStatus(status);
				}
			}
			
			tempidentification_=(*ptr_);		
		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to store the ProteinIdentification data. "<<e.what()<<endl;
		}
		
	}
	
	void PeptideIdentificationVisualizer::reject_()
	{
		
		try
		{
			load(*ptr_, tree_id_);
		}
		catch(exception e)
		{
			cout<<"Error while trying to restore original protein ProteinIdentification data. "<<e.what()<<endl;
		}
		
	}

}
