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


#include <OpenMS/VISUAL/VISUALIZER/SampleVisualizer.h>

//QT
#include <QtGui/QTextEdit>
#include <QtGui/QValidator>
#include <QtGui/QLineEdit>
#include <QtGui/QComboBox>

#include <iostream>

using namespace std;

namespace OpenMS
{
	
	SampleVisualizer::SampleVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Sample>()
	{
		addLabel("Modify Sample information");		
		addSeparator();
	  addLineEdit(samplename_, "Name" );
		addLineEdit(samplenumber_, "Number" );
		addLineEdit(sampleorganism_, "Organism" );
	  addTextEdit(samplecomment_, "Comment");
		addComboBox(samplestate_, "State");
		addDoubleLineEdit(samplemass_,"Mass (in mg)");
		addDoubleLineEdit(samplevolume_, "Volume (in ml)");
		addDoubleLineEdit(sampleconcentration_, "Concentration (in mg/ml)");
		
		finishAdding_();
	}
	
	void SampleVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox(samplestate_,& temp_.NamesOfSampleState[temp_.getState()] ,1);
		}
		else
		{
			fillComboBox(samplestate_, temp_.NamesOfSampleState , Sample::SIZE_OF_SAMPLESTATE);
			samplestate_->setCurrentIndex(temp_.getState());
		}
		
		samplename_->setText(temp_.getName().c_str());
		samplenumber_->setText(temp_.getNumber().c_str());
		sampleorganism_->setText(temp_.getOrganism().c_str());
	  samplecomment_->setText(temp_.getComment().c_str());
		
		samplemass_->setText(String(temp_.getMass()).c_str()   );
		samplevolume_->setText(String(temp_.getVolume()).c_str());
		sampleconcentration_->setText(String(temp_.getConcentration()).c_str() );
	}
	
	void SampleVisualizer::store()
	{
		ptr_->setName(samplename_->text().toStdString());
		ptr_->setNumber(samplenumber_->text().toStdString());
		ptr_->setOrganism(sampleorganism_->text().toStdString());
		ptr_->setComment(samplecomment_-> toPlainText().toStdString());
		ptr_->setState((Sample::SampleState)samplestate_->currentIndex());		
		String m(samplemass_->text().toStdString());
		ptr_->setMass(m.toFloat());
		String v(samplevolume_->text().toStdString()) ;
		ptr_->setVolume(v.toFloat());
		String c(sampleconcentration_->text().toStdString()) ;
		ptr_->setConcentration(c.toFloat());
		
		temp_=(*ptr_);
	}
	
	void SampleVisualizer::undo_()
	{
		update_();
	}

}
