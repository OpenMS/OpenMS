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
	
	void Spectrum3DWindow::showGoToDialog()
	{
	}	

	void Spectrum3DWindow::createContextMenu_()
	{
		if (context_menu_!=0)
		{
			delete(context_menu_);
		}
		
		//create menu
		context_menu_ = new QPopupMenu(this);
	 	SignedInt item;
	
		//intensity distrubution
		context_menu_->insertItem("intensity distribution",widget_,SLOT(showIntensityDistribution()));
		context_menu_->insertSeparator();
	
		//legend menu
		QPopupMenu* legend_menu = new QPopupMenu(context_menu_);
		item = legend_menu->insertItem("shown",widget(),SLOT(showLegend(int)),0,1);
		if (widget()->isLegendShown()) legend_menu->setItemEnabled(item,false);
		item = legend_menu->insertItem("hidden",widget(),SLOT(showLegend(int)),0,0);
		if (!widget()->isLegendShown()) legend_menu->setItemEnabled(item,false);
		context_menu_->insertItem("legend",legend_menu);
		context_menu_->insertSeparator();	
	
		//Preferences
		context_menu_->insertItem("Preferences",this,SLOT(showPreferences_()));
		context_menu_->insertSeparator();	
		//SpectrumWindow::createContextMenu_();
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
