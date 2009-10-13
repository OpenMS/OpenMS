// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
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

	PeptideHitVisualizer::PeptideHitVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<PeptideHit>()
	{
		addLineEdit_(peptidehit_score_, "Score" );
		addLineEdit_(peptidehit_charge_, "Charge" );
		addLineEdit_(peptidehit_rank_, "Rank" );
		addTextEdit_(peptidehit_sequence_, "Sequence" );
		
		finishAdding_();
	}
	
	void PeptideHitVisualizer::update_()
	{
	  peptidehit_score_->setText(String(temp_.getScore()).c_str() );
		peptidehit_score_->setReadOnly(true);
		peptidehit_charge_->setText(String(temp_.getCharge()).c_str() );
		peptidehit_charge_->setReadOnly(true);
	  peptidehit_rank_->setText(String(temp_.getRank()).c_str());
		peptidehit_rank_->setReadOnly(true);
		peptidehit_sequence_->setText(temp_.getSequence().toString().c_str()); 
		peptidehit_sequence_->setReadOnly(true);			
	}
	
	void PeptideHitVisualizer::store()
	{
		//TODO?
		(*ptr_) =temp_ ;
	}
	
	void PeptideHitVisualizer::undo_()
	{
		update_();
	}

}
