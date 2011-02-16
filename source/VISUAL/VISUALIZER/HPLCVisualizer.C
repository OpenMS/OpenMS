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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/HPLCVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QValidator>
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	HPLCVisualizer::HPLCVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<HPLC>()
	{
		addLabel_("Modify HPLC information");		
		addSeparator_();
		addLineEdit_(hplcinstrument_, "Instrument" );
		addLineEdit_(hplccolumn_, "Column" );
		addIntLineEdit_(hplctemperature_, "Temperature (in °C)" );
		addIntLineEdit_(hplcpressure_, "Pressure (in bar)" );
		addIntLineEdit_(hplcflux_, "Flux (in µl/sec)" );
		addTextEdit_(hplccomment_, "Comment");
			
		finishAdding_();
	}
	
	
	void HPLCVisualizer::update_()
	{
	  hplcinstrument_->setText(temp_.getInstrument().c_str());
		hplccolumn_->setText(temp_.getColumn().c_str() );
	  hplctemperature_->setText(String(temp_.getTemperature()).c_str());
		hplcpressure_->setText(String(temp_.getPressure()).c_str());
		hplcflux_->setText(String(temp_.getFlux()).c_str());
		hplccomment_->setText(temp_.getComment().c_str()); 
	}
	
	void HPLCVisualizer::store()
	{
		ptr_->setInstrument(hplcinstrument_->text());
		ptr_->setColumn(hplccolumn_->text());
		ptr_->setTemperature(hplctemperature_->text().toInt() );
		ptr_->setPressure(hplcpressure_->text().toInt() );
		ptr_->setFlux(hplcflux_->text().toInt());
		ptr_->setComment(hplccomment_->toPlainText());
		temp_ = (*ptr_);
	}
	
	void HPLCVisualizer::undo_()
	{
		update_();
	}

}
