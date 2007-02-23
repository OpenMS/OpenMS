// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DWindow.h>
#include <OpenMS/VISUAL/Spectrum3DWidget.h>

namespace OpenMS 
{
	
	Spectrum3DWindow::Spectrum3DWindow(QWidget* parent)
		: SpectrumWindow(parent)
	{
		setWidget_(new Spectrum3DWidget(this));
		setCentralWidget(widget());
		connectWidgetSignals(widget());
	}	
	
	/// Destructor
	Spectrum3DWindow::~Spectrum3DWindow()
	{
	}
	
	void Spectrum3DWindow::showGoToDialog()
	{
	}
		
	PreferencesDialogPage* Spectrum3DWindow::createPreferences(QWidget* parent)
	{  
	  return widget()->createPreferences(parent);
	}
	
	Spectrum3DWidget* Spectrum3DWindow::widget()
	{
	 return static_cast<Spectrum3DWidget*>(widget_);
	}

	

} //namespace
