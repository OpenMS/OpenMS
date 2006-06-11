// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: Spectrum3DWindow.C,v 1.7 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DWindow.h>
#include<OpenMS/VISUAL/Spectrum3DWidget.h>
#include <iostream.h>
#include <qpopupmenu.h>
namespace OpenMS 
{
	
	Spectrum3DWindow::Spectrum3DWindow(QWidget* parent, const char* name, WFlags f)
		: SpectrumWindow(parent,name,f)
	{
		setWidget_(new Spectrum3DWidget(this));
		setCentralWidget(widget());
		connectWidgetSignals(widget());
		connect(widget(), SIGNAL(contextMenu(QPoint)), this, SLOT(showContextMenu_(QPoint)));
		
	}	
	
	/// Destructor
	Spectrum3DWindow::~Spectrum3DWindow()
	{
	}


	void Spectrum3DWindow::setMainPreferences(const Param& pref)
	{ 
	  widget()->setMainPreferences(pref);
	}
		void Spectrum3DWindow::showGoToDialog()
	{
	}	

	void Spectrum3DWindow::createContextMenu_()
	{
		SpectrumWindow::createContextMenu_();
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
