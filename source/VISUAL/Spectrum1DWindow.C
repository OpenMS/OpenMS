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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

//STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWindow.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

using namespace std;

namespace OpenMS
{

	Spectrum1DWindow::Spectrum1DWindow(QWidget* parent)
		: SpectrumWindow(parent)
	{
		setWidget_(new Spectrum1DWidget(this));
		widget()->PreferencesManager::setParent(this);
		widget()->show();
		setCentralWidget(widget());
		connectWidgetSignals(widget());
	}
	
	Spectrum1DWidget* Spectrum1DWindow::widget()
	{
		return static_cast<Spectrum1DWidget*>(widget_);
	}
	
	Spectrum1DWindow::~Spectrum1DWindow()
	{
	
	}
	
	PreferencesDialogPage* Spectrum1DWindow::createPreferences(QWidget* parent)
	{
		return widget()->createPreferences(parent);
	}
	
	void Spectrum1DWindow::showGoToDialog()
	{
	  Spectrum1DGoToDialog goto_dialog(this);
	  goto_dialog.setMin(widget()->canvas()->getDataRange().minX());
	  goto_dialog.setMax(widget()->canvas()->getDataRange().maxX());
	  if (goto_dialog.exec())
	  {
	  	widget()->canvas()->setVisibleArea(SpectrumCanvas::AreaType(goto_dialog.getMin(),0,goto_dialog.getMax(),0));
		}
	}

} //namespace


