// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software you can redistribute it and/or
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
// $Maintainer:  Marc Sturm $
// --------------------------------------------------------------------------s


#include <OpenMS/VISUAL/VISUALIZER/ProteinIdentificationVisualizer.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>
#include <OpenMS/VISUAL/MSMetaDataExplorer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>
#include <QtGui/QPushButton>
#include <QtGui/QComboBox>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	//Constructor
	ProteinIdentificationVisualizer::ProteinIdentificationVisualizer(bool editable, QWidget *parent, MSMetaDataExplorer *caller) : BaseVisualizer(editable, parent)
	{
		type_="Identification";
		pidv_caller_= caller;
	     
		addLineEdit(identifier_, "Identifier<br>(of corresponding PeptideIdentifications)" );
		
		addSeperator();
		addLineEdit(engine_, "Search engine" );
		addLineEdit(engine_version_, "Search engine version" );
		addLineEdit(identification_date_, "Date of search" );
		addLineEdit(score_type_, "Score type" );
		addBooleanComboBox(higher_better_,"Higher score is better"); 
		addDoubleLineEdit(identification_threshold_, "Protein significance threshold" );	
		
		addSeperator();
		addLabel("Search Parameters:");
		addLineEdit(db_, "Database name" );
		addLineEdit(db_version_, "Database version" );
		addLineEdit(taxonomy_, "Taxonomy restriction" );
		addLineEdit(charges_, "Allowed charges" );
		addIntLineEdit(missed_cleavages_, "Missed Cleavages" );
		addDoubleLineEdit(peak_tolerance_, "Fragment ion mass tolerance" );
		addDoubleLineEdit(precursor_tolerance_, "Precursor ion mass tolerance" );
		addComboBox(mass_type_, "Mass type" );
		addComboBox(enzyme_, "Digestion enzyme" );
		
		addSeperator();       
		addLabel("Show protein hits with score equal or better than a threshold.");
		QPushButton* button;
		addLineEditButton("Score threshold", filter_threshold_, button, "Filter");
		connect(button, SIGNAL(clicked()), this, SLOT(updateTree_()) );
	
		finishAdding_();
	}



	void ProteinIdentificationVisualizer::load(ProteinIdentification &s, int tree_item_id)
	{
	  //Pointer to current object to keep track of the actual object
		ptr_ = &s;
		
		// id of the item in the tree
		tree_id_ = tree_item_id;
		
		//Copy of current object for restoring the original values
		tempidentification_=s;
	  
	  String str;
	  tempidentification_.getDateTime().get(str);
		identification_date_->setText(str.toQString()); 
		identification_threshold_->setText(QString::number(tempidentification_.getSignificanceThreshold()));			
		identifier_->setText(tempidentification_.getIdentifier().toQString());
		engine_->setText(tempidentification_.getSearchEngine().toQString());
		engine_version_->setText(tempidentification_.getSearchEngineVersion().toQString());
		score_type_->setText(tempidentification_.getScoreType().toQString());
		higher_better_->setCurrentIndex(tempidentification_.isHigherScoreBetter());
		
		db_->setText(tempidentification_.getSearchParameters().db.toQString());
		db_version_->setText(tempidentification_.getSearchParameters().db_version.toQString());
		taxonomy_->setText(tempidentification_.getSearchParameters().taxonomy.toQString());
		charges_->setText(tempidentification_.getSearchParameters().charges.toQString());
		missed_cleavages_->setText(QString::number(tempidentification_.getSearchParameters().missed_cleavages));
		peak_tolerance_->setText(QString::number(tempidentification_.getSearchParameters().peak_mass_tolerance));
		precursor_tolerance_->setText(QString::number(tempidentification_.getSearchParameters().precursor_tolerance));

		if(! isEditable())
		{
			fillComboBox(mass_type_, &ProteinIdentification::NamesOfPeakMassType[tempidentification_.getSearchParameters().mass_type], 1);
			fillComboBox(enzyme_, &ProteinIdentification::NamesOfDigestionEnzyme[tempidentification_.getSearchParameters().enzyme], 1);
		}
		else
		{
			fillComboBox(mass_type_, ProteinIdentification::NamesOfPeakMassType, ProteinIdentification::SIZE_OF_PEAKMASSTYPE);
			fillComboBox(enzyme_, ProteinIdentification::NamesOfDigestionEnzyme, ProteinIdentification::SIZE_OF_DIGESTIONENZYME);
			
			enzyme_->setCurrentIndex(tempidentification_.getSearchParameters().enzyme); 
			mass_type_->setCurrentIndex(tempidentification_.getSearchParameters().mass_type); 
		}
	}
	
	void ProteinIdentificationVisualizer::updateTree_()
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
	
	void ProteinIdentificationVisualizer::store_()
	{
		try
		{
			ptr_->setSearchEngine(engine_->text());
			ptr_->setSearchEngineVersion(engine_version_->text());
			ptr_->setIdentifier(identifier_->text());
			ptr_->setSignificanceThreshold(identification_threshold_->text().toFloat());
			ptr_->setScoreType(score_type_->text());
			ptr_->setHigherScoreBetter(higher_better_->currentIndex());
			//date
			DateTime date;
			try
			{
				date.set(identification_date_->text());
				ptr_->setDateTime(date);
			}
			catch(exception& e)
			{
				if(date.isNull())
				{
					std::string status= "Format of date in PROTEINIDENTIFICATION is not correct.";
					emit sendStatus(status);
				}
			}
			
			//search parameters
			ProteinIdentification::SearchParameters tmp = ptr_->getSearchParameters();
			tmp.db = db_->text();
			tmp.db_version = db_version_->text();
			tmp.taxonomy = taxonomy_->text();
			tmp.charges = charges_->text();
			tmp.missed_cleavages = missed_cleavages_->text().toInt();
			tmp.peak_mass_tolerance = peak_tolerance_->text().toFloat();
			tmp.precursor_tolerance = precursor_tolerance_->text().toFloat();
			tmp.enzyme = (ProteinIdentification::DigestionEnzyme)(enzyme_->currentIndex());
			tmp.mass_type = (ProteinIdentification::PeakMassType)(mass_type_->currentIndex());
			ptr_->setSearchParameters(tmp);
			
			tempidentification_=(*ptr_);		
		}
		catch(exception& e)
		{
			std::cout<<"Error while trying to store the new protein ProteinIdentification data. "<<e.what()<<endl;
		}
		
	}
	
	void ProteinIdentificationVisualizer::reject_()
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
