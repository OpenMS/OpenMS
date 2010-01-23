// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copymain (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/VISUAL/VISUALIZER/MetaInfoDescriptionVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QTextEdit>

using namespace std;

namespace OpenMS
{
	
	MetaInfoDescriptionVisualizer::MetaInfoDescriptionVisualizer(bool editable, QWidget* parent)
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<MetaInfoDescription>()
	{
		addLabel_("Modify MetaInfoDescription information");		
		addSeparator_();
		addLineEdit_(metainfodescription_name_, "Name of peak annotations" );
			
		finishAdding_();
	}
	
	void MetaInfoDescriptionVisualizer::update_()
	{
	  metainfodescription_name_->setText(temp_.getName().c_str() );
	}
	
	void MetaInfoDescriptionVisualizer::store()
	{
		ptr_->setName(metainfodescription_name_->text().toStdString());
					
		temp_ = (*ptr_);
	}
	
	void MetaInfoDescriptionVisualizer::undo_()
	{
		update_();
	}

}
