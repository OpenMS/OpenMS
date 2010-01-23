// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------s

#include <OpenMS/VISUAL/VISUALIZER/PeptideIdentificationVisualizer.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{

	PeptideIdentificationVisualizer::PeptideIdentificationVisualizer(bool editable, QWidget* parent, MetaDataBrowser* caller) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<PeptideIdentification>()
	{
		pidv_caller_= caller;
		
		addLineEdit_(identifier_, "Identifier<br>(of corresponding ProteinIdentification)" );
		addSeparator_();   
		
		addLineEdit_(score_type_, "Score type" );
		addBooleanComboBox_(higher_better_,"Higher score is better"); 
		addDoubleLineEdit_(identification_threshold_, "Peptide significance threshold" );	
		
		addSeparator_();       
		addLabel_("Show peptide hits with score equal or better than a threshold.");
		QPushButton* button;
		addLineEditButton_("Score threshold", filter_threshold_, button, "Filter");
		connect(button, SIGNAL(clicked()), this, SLOT(updateTree_()) );
		
		finishAdding_();
	}

	void PeptideIdentificationVisualizer::load(PeptideIdentification& s, int tree_item_id)
	{
		ptr_ = &s;
		temp_ = s;
		
		// id of the item in the tree
		tree_id_ = tree_item_id;
		
	  identifier_->setText(temp_.getIdentifier().toQString());
		identification_threshold_->setText(QString::number(temp_.getSignificanceThreshold()));					
		score_type_->setText(temp_.getScoreType().toQString());
		higher_better_->setCurrentIndex(temp_.isHigherScoreBetter());
	}
	
	void PeptideIdentificationVisualizer::updateTree_()
	{
		if (filter_threshold_->text()!="")
		{
			pidv_caller_->filterHits_(filter_threshold_->text().toDouble(),temp_.isHigherScoreBetter(),tree_id_ );
		}
		else
		{
			pidv_caller_->showAllHits_(tree_id_);
		}
	}
	
	void PeptideIdentificationVisualizer::store()
	{
		ptr_->setIdentifier(identifier_->text());
		ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
		ptr_->setScoreType(score_type_->text());
		ptr_->setHigherScoreBetter(higher_better_->currentIndex());
		
		temp_=(*ptr_);		
	}
	
	void PeptideIdentificationVisualizer::undo_()
	{
		load(*ptr_, tree_id_);
	}

}
