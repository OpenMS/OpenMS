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
//  License as published by the Free SpectrumSettings Foundation; either
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

#include <OpenMS/VISUAL/VISUALIZER/SpectrumSettingsVisualizer.h>

//QT
#include <QtGui/QComboBox>
#include <QtGui/QTextEdit>

//STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	SpectrumSettingsVisualizer::SpectrumSettingsVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<SpectrumSettings>()
	{
		addLabel_("Modify the settings of the spectrum.");	
		addSeparator_();  
		addComboBox_(spectrumsettings_type_, "Type of spectrum");
		addTextEdit_(spectrumsettings_comment_, "Comment");
			
		finishAdding_();
	}
	
	void SpectrumSettingsVisualizer::update_()
	{
		if(! isEditable())
		{
			fillComboBox_(spectrumsettings_type_,& temp_.NamesOfSpectrumType[temp_.getType()] , 1);
		}
		else
		{
			fillComboBox_(spectrumsettings_type_, temp_.NamesOfSpectrumType , SpectrumSettings::SIZE_OF_SPECTRUMTYPE);
			spectrumsettings_type_->setCurrentIndex(temp_.getType()); 
		}
		
		spectrumsettings_comment_->setText(temp_.getComment().c_str());
	}
	
	void SpectrumSettingsVisualizer::store()
	{
		ptr_->setType((SpectrumSettings::SpectrumType)spectrumsettings_type_->currentIndex());			
		ptr_->setComment(spectrumsettings_comment_->toPlainText().toStdString());
		
		temp_=(*ptr_);
	}
	
	void SpectrumSettingsVisualizer::undo_()
	{
		update_();
	}

}
