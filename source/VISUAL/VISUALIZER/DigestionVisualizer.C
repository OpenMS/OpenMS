// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/VISUALIZER/DigestionVisualizer.h>

//QT
#include <QtGui/QValidator>
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	DigestionVisualizer::DigestionVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Digestion>()
	{
		addLabel_("Modify Digestion information");		
		addSeparator_();
		addLineEdit_(treatmenttype_, "Treatment type" );
		addTextEdit_(treatmentcomment_, "Comment" );
		addLineEdit_(digestionenzyme_, "Enzyme" );
		addDoubleLineEdit_(digestiontime_, "Digestion time (in minutes)" );
		addDoubleLineEdit_(digestiontemperature_, "Temperature (in °C)" );
		addDoubleLineEdit_(digestionPH_, "PH" );
		
		finishAdding_();
	}
	
	void DigestionVisualizer::update_()
	{
		treatmenttype_->setText(temp_.getType().c_str());
		treatmenttype_->setReadOnly(true);
		treatmentcomment_->setText(temp_.getComment().c_str());
	  digestionenzyme_->setText(temp_.getEnzyme().c_str());
		digestiontime_->setText(String(temp_.getDigestionTime()).c_str() );
	  digestiontemperature_->setText(String(temp_.getTemperature()).c_str());
		digestionPH_->setText(String(temp_.getPh()).c_str()); 
	}
	
	void DigestionVisualizer::store()
	{
		ptr_->setComment(treatmentcomment_->toPlainText());
		ptr_->setEnzyme(digestionenzyme_->text());
		ptr_->setDigestionTime(digestiontime_->text().toFloat());
		ptr_->setTemperature(digestiontime_->text().toFloat());
		ptr_->setPh(digestiontime_->text().toFloat());
		
		temp_ = (*ptr_);
	}
	
	void DigestionVisualizer::undo_()
	{
		update_();
	}

}
