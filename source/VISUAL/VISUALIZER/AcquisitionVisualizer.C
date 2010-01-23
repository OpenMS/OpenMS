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
#include <OpenMS/VISUAL/VISUALIZER/AcquisitionVisualizer.h>

//QT
#include <QtGui/QLineEdit>
#include <QtGui/QValidator>

// STL
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	AcquisitionVisualizer::AcquisitionVisualizer(bool editable, QWidget* parent) 
		: BaseVisualizerGUI(editable, parent),
			BaseVisualizer<Acquisition>()
	{
	  
		addLabel_("Show Acquisition information");		
		addSeparator_();
		addIntLineEdit_(acquisitionnumber_, "Identifier of the scan" );
		acquisitionnumber_->setReadOnly(true);
			
		finishAdding_();
	}
	
	
	void AcquisitionVisualizer::update_()
	{
	  acquisitionnumber_->setText(temp_.getIdentifier().toQString() );
	}
	
	void AcquisitionVisualizer::store()
	{
		ptr_->setIdentifier(acquisitionnumber_->text());
					
		temp_ = (*ptr_);
	}
	
	void AcquisitionVisualizer::undo_()
	{
		update_();
	}

}
