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
#include <OpenMS/VISUAL/VISUALIZER/PeptideHitVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{

//Constructor
PeptideHitVisualizer::PeptideHitVisualizer(bool editable, QWidget *parent) : BaseVisualizer(editable, parent)
{
  
	addLabel("Show PeptideHit information.");		
	addSeperator();
	addLineEdit(peptidehit_score_, "Score" );
	addLineEdit(peptidehit_score_type_, "Score type" );
	addLineEdit(peptidehit_charge_, "Charge" );
	addLineEdit(peptidehit_rank_, "Rank" );
	addTextEdit(peptidehit_sequence_, "Sequence" );
	
	
	addSeperator();
	
	
}


void PeptideHitVisualizer::load(PeptideHit &h)
{
  ptr_ = &h;
	
	//Copy of current object for restoring the original values
	tempPeptideHit_=h;
  peptidehit_score_->setText(String(tempPeptideHit_.getScore()).c_str() );
	peptidehit_score_->setReadOnly(true);
	peptidehit_score_type_->setText(tempPeptideHit_.getScoreType().c_str() );
	peptidehit_score_type_->setReadOnly(true);
	peptidehit_charge_->setText(String(tempPeptideHit_.getCharge()).c_str() );
	peptidehit_charge_->setReadOnly(true);
  peptidehit_rank_->setText(String(tempPeptideHit_.getRank()).c_str());
	peptidehit_rank_->setReadOnly(true);
	peptidehit_sequence_->setText(tempPeptideHit_.getSequence().c_str()); 
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

}
