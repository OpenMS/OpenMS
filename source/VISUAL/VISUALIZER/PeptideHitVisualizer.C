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
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>
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
PeptideHitVisualizer::PeptideHitVisualizer(QWidget *parent, const char *name) : BaseVisualizer(parent, name)
{
  
	addLabel("Show PeptideHit information.");		
	addSeperator();
	addLineEdit(peptidehit_score_, "Score" );
	addLineEdit(peptidehit_score_type_, "Score type" );
	addLineEdit(peptidehit_rank_, "Rank" );
	addTextEdit(peptidehit_sequence_, "Sequence" );
	
	
	addSeperator();
	//addLabel("Save changes or restore original data.");
	//addHorizontalButtons(savebutton_, "Save",  cancelbutton_, "Cancel");
	
  //connect(savebutton_, SIGNAL(clicked()), this, SLOT(store()) );
	//connect(cancelbutton_, SIGNAL(clicked()), this, SLOT(reject()) );
	
}


void PeptideHitVisualizer::load(PeptideHit &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempPeptideHit_=h;
  peptidehit_score_->setText(String(tempPeptideHit_.getScore()) );
	peptidehit_score_->setReadOnly(true);
	peptidehit_score_type_->setText(tempPeptideHit_.getScoreType() );
	peptidehit_score_type_->setReadOnly(true);
  peptidehit_rank_->setText(String(tempPeptideHit_.getRank()));
	peptidehit_rank_->setReadOnly(true);
	peptidehit_sequence_->setText(tempPeptideHit_.getSequence()); 
	peptidehit_sequence_->setReadOnly(true);
	
			
}

void PeptideHitVisualizer::store()
{
	try
	{
	
		//Information of PeptideHit is not to be changed.
		//(*ptr_) =tempPeptideHit_ ;
	}
	catch(exception& e)
	{
		std::cout<<"Error while trying to store the new PeptideHit data. "<<e.what()<<endl;
	}
}

void PeptideHitVisualizer::reject()
{
	try
	{
		load(tempPeptideHit_);
	}
	catch(exception e)
	{
		cout<<"Error while trying to restore original PeptideHit data. "<<e.what()<<endl;
	} 
}
