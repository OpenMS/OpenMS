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
//  License as published by the Free ExperimentalSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/ExperimentalSettingsVisualizer.h>
#include <OpenMS/DATASTRUCTURES/Date.h>
#include <OpenMS/VISUAL/MetaDataBrowser.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>
#include <QtGui/QComboBox>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	ExperimentalSettingsVisualizer::ExperimentalSettingsVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<ExperimentalSettings>()
	{
		addLabel_("Modify the settings of the experiment.");	
		addSeparator_();
		addLineEdit_(datetime_, "Date and time of experiment");
		addTextEdit_(comment_, "Comment");
		addLineEdit_(fraction_identifier_, "Fraction identifier");
		
		finishAdding_();
	}
	
	void ExperimentalSettingsVisualizer::update_()
	{		
	  datetime_->setText(temp_.getDateTime().get().c_str()); 
		comment_->setText(temp_.getComment().c_str());
		fraction_identifier_->setText(temp_.getFractionIdentifier().c_str());
	}
	
	void ExperimentalSettingsVisualizer::store()
	{
		DateTime date;
		try
		{
			date.set(datetime_->text());
			ptr_->setDateTime(date);
		}
		catch(exception& /*e*/)
		{
			if(date.isNull())
			{
				std::string status= "Format of date in EXPERIMENTALSETTINGS is not correct.";
				emit sendStatus(status);
			}
		}
		
		ptr_->setComment(comment_->toPlainText());
		ptr_->setFractionIdentifier(fraction_identifier_->text());
		
		temp_=(*ptr_);
	}
	
	void ExperimentalSettingsVisualizer::undo_()
	{
		update_();
	}

}
