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
// $Id: Spectrum1DWindow.C,v 1.19 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
//


// QT
#include <qpopupmenu.h>

//STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWindow.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>

using namespace std;

namespace OpenMS
{

	Spectrum1DWindow::Spectrum1DWindow(QWidget* parent, const char* name, WFlags f)
		: SpectrumWindow(parent,name,f)
	{
		setWidget_(new Spectrum1DWidget(this));
		widget()->setParent(this);
		widget()->show();
		setCentralWidget(widget());
		connectWidgetSignals(widget());
		//connect(widget(), SIGNAL(visibleAreaChanged(double, double)), this, SIGNAL(loXHiXChanged(double, double)));
	}
	
	Spectrum1DWidget* Spectrum1DWindow::widget()
	{
		return static_cast<Spectrum1DWidget*>(widget_);
	}
	
	Spectrum1DWindow::~Spectrum1DWindow()
	{
	
	}
	
	void Spectrum1DWindow::createContextMenu_()
	{
		SpectrumWindow::createContextMenu_();
		
		SignedInt item;
		//axis modes
		QPopupMenu* axis_menu = new QPopupMenu(context_menu_);
		item = axis_menu->insertItem("absolute",widget(),SLOT(intensityAxisAbsolute()));
		if (widget()->getLabelMode() == Spectrum1DCanvas::LM_XABSOLUTE_YABSOLUTE) axis_menu->setItemEnabled(item,false);
		item = axis_menu->insertItem("relative",widget(),SLOT(intensityAxisRelative()));
		if (widget()->getLabelMode() == Spectrum1DCanvas::LM_XABSOLUTE_YPERCENT) axis_menu->setItemEnabled(item,false);
		context_menu_->insertItem("y axis label",axis_menu);
		context_menu_->insertSeparator();
	}
	
	void Spectrum1DWindow::switchAxis(bool b)
	{
		widget()->switchAxis(b);
	}
	
	void Spectrum1DWindow::setMirroredXAxis(bool b)
	{
		widget()->setMirroredXAxis(b);
	}
	
	void Spectrum1DWindow::setMirroredYAxis(bool b)
	{
		widget()->setMirroredYAxis(b);
	}
	
	int Spectrum1DWindow::getDrawMode()
	{
		return widget()->canvas()->getDrawMode();
	}
	
	void Spectrum1DWindow::setDrawMode(QAction* a)
	{
		widget()->setDrawMode(a);
	}
	
	void Spectrum1DWindow::setMainPreferences(const Param& prefs)
	{
		widget()->setMainPreferences(prefs);
	}
	
	PreferencesDialogPage* Spectrum1DWindow::createPreferences(QWidget* parent)
	{
		return widget()->createPreferences(parent);
	}
	
	void Spectrum1DWindow::showGoToDialog()
	{
	  Spectrum1DGoToDialog goToDialog(this, "Spectrum1DGoToDialog");
	  const SpectrumCanvas::AreaType& visible_area = widget()->canvas()->getVisibleArea();
	  goToDialog.setMinPosition(visible_area.minX());
	  goToDialog.setMaxPosition(visible_area.maxX());
	  goToDialog.exec();
	  widget()->setVisibleArea(goToDialog.getMinPosition(),goToDialog.getMaxPosition());
	}
	
	bool Spectrum1DWindow::getSnapToMax()
	{
		return ((Spectrum1DWidget*)widget_)->getSnapToMax();
	}
	
	void Spectrum1DWindow::setSnapToMax(bool b)
	{
		((Spectrum1DWidget*)widget_)->setSnapToMax(b);
	}

} //namespace


