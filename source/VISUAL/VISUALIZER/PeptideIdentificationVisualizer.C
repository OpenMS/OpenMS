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
#include <QtGui/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

	PeptideIdentificationVisualizer::PeptideIdentificationVisualizer(bool editable, QWidget *parent, MSMetaDataExplorer *caller) 
		: BaseVisualizer(editable,parent)
	{
		type_="PeptideIdentification";
		pidv_caller_= caller;
		
		addLineEdit(identifier_, "Identifier<br>(of corresponding ProteinIdentification)" );
		addSeperator();   
		
		addLineEdit(score_type_, "Score type" );
		addBooleanComboBox(higher_better_,"Higher score is better"); 
		addDoubleLineEdit(identification_threshold_, "Peptide significance threshold" );	
		
		addSeperator();       
		addLabel("Show peptide hits with score equal or better than a threshold.");
		QPushButton* button;
		addLineEditButton("Score threshold", filter_threshold_, button, "Filter");
		connect(button, SIGNAL(clicked()), this, SLOT(updateTree_()) );
		
		finishAdding_();
	}

	void PeptideIdentificationVisualizer::load(PeptideIdentification &s, int tree_item_id)
	{
	  //Pointer to current object to keep track of the actual object
		ptr_ = &s;
		
		// id of the item in the tree
		tree_id_ = tree_item_id;
		
		//Copy of current object for restoring the original values
		tempidentification_=s;
	  
	  identifier_->setText(tempidentification_.getIdentifier().toQString());
		identification_threshold_->setText(QString::number(tempidentification_.getSignificanceThreshold()));					
		score_type_->setText(tempidentification_.getScoreType().toQString());
		higher_better_->setCurrentIndex(tempidentification_.isHigherScoreBetter());
	}
	
	void PeptideIdentificationVisualizer::updateTree_()
	{
		if (filter_threshold_->text()!="")
		{
			pidv_caller_->filterHits_(filter_threshold_->text().toDouble(),tempidentification_.isHigherScoreBetter(),tree_id_ );
		}
		else
		{
			pidv_caller_->showAllHits_(tree_id_);
		}
	}
	
	void PeptideIdentificationVisualizer::store_()
	{
		try
		{
			ptr_->setIdentifier(identifier_->text());
			ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
			ptr_->setScoreType(score_type_->text());
			ptr_->setHigherScoreBetter(higher_better_->currentIndex());
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
