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
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionInfoVisualizer.h>

// QT
#include <QtGui/QValidator>
#include <QtGui/QLineEdit>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	AcquisitionInfoVisualizer::AcquisitionInfoVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<AcquisitionInfo>()
	{
		addLabel_("Show AcquisitionInfo information");		
		addSeparator_();
		addIntLineEdit_(acquisitioninfo_method_, "Method of combination" );
		
		finishAdding_();
	}
	
	void AcquisitionInfoVisualizer::update_()
	{
	  acquisitioninfo_method_->setText( temp_.getMethodOfCombination().c_str() );
	}
	
	void AcquisitionInfoVisualizer::store()
	{
		ptr_->setMethodOfCombination(acquisitioninfo_method_->text());
					
		temp_ = (*ptr_);
	}
	
	void AcquisitionInfoVisualizer::undo_()
	{
		update_();
	}
}
