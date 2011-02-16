// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/VISUALIZER/ProteinHitVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QLineEdit>

#include <iostream>

using namespace std;

namespace OpenMS
{

	ProteinHitVisualizer::ProteinHitVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<ProteinHit>()
	{
		addLineEdit_(proteinhit_score_, "Score" );
		addLineEdit_(proteinhit_rank_, "Rank" );
		addLineEdit_(proteinhit_accession_, "Accession" );
		addTextEdit_(proteinhit_sequence_, "Sequence" );
			
		finishAdding_();
	}
	
	void ProteinHitVisualizer::update_()
	{
	  proteinhit_score_->setText(String(temp_.getScore()).c_str() );
		proteinhit_score_->setReadOnly(true);
	  proteinhit_rank_->setText(String(temp_.getRank()).c_str());
		proteinhit_rank_->setReadOnly(true);
		proteinhit_accession_->setText(temp_.getAccession().c_str());
		proteinhit_accession_->setReadOnly(true);
		proteinhit_sequence_->setText(temp_.getSequence().c_str()); 
		proteinhit_sequence_->setReadOnly(true);
	}
	
	void ProteinHitVisualizer::store()
	{
		//TODO?
		(*ptr_) =temp_ ;
	}
	
	void ProteinHitVisualizer::undo_()
	{
		update_();
	}

}
