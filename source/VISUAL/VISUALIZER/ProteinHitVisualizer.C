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
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
//#include <OpenMS/VISUAL/DataTable.h>

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
ProteinHitVisualizer::ProteinHitVisualizer(bool editable, QWidget *parent, const char *name) : BaseVisualizer(editable, parent, name)
{
  
	addLabel("Show ProteinHit information.");		
	addSeperator();
	addLineEdit(proteinhit_score_, "Score" );
	addLineEdit(proteinhit_score_type_, "Score type" );
	addLineEdit(proteinhit_rank_, "Rank" );
	addLineEdit(proteinhit_accession_, "Accession" );
	addLineEdit(proteinhit_accession_type_, "Accession type" );
	addTextEdit(proteinhit_sequence_, "Sequence" );
	
	
	addSeperator();
	
}


void ProteinHitVisualizer::load(ProteinHit &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempProteinHit_=h;
  proteinhit_score_->setText(String(tempProteinHit_.getScore()) );
	proteinhit_score_->setReadOnly(true);
	proteinhit_score_type_->setText(tempProteinHit_.getScoreType() );
	proteinhit_score_type_->setReadOnly(true);
  proteinhit_rank_->setText(String(tempProteinHit_.getRank()));
	proteinhit_rank_->setReadOnly(true);
	proteinhit_accession_->setText(tempProteinHit_.getAccession());
	proteinhit_accession_->setReadOnly(true);
	proteinhit_accession_type_->setText(tempProteinHit_.getAccessionType());
	proteinhit_accession_type_->setReadOnly(true);
	proteinhit_sequence_->setText(tempProteinHit_.getSequence()); 
	proteinhit_sequence_->setReadOnly(true);
	
			
}

void ProteinHitVisualizer::store()
{
	try
	{
	
		//Information of ProteinHit is not to be changed.
		//(*ptr_) =tempProteinHit_ ;
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new ProteinHit data. "<<e.what()<<endl;
	}
}

void ProteinHitVisualizer::reject()
{
	try
	{
		load(tempProteinHit_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original ProteinHit data. "<<e.what()<<endl;
	} 
}
