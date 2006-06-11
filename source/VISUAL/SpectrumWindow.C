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
// $Id: SpectrumWindow.C,v 1.13 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream.h>

#include <OpenMS/VISUAL/SpectrumWindow.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>

#include <qpopupmenu.h>

namespace OpenMS
{

	SpectrumWindow::SpectrumWindow(QWidget* parent, const char* name, WFlags f)  
		: QMainWindow(parent,name,f),
		PreferencesManager(),
		context_menu_(0)
	{
		setMinimumSize(200,200);	// prevents errors caused by too small width,height values
	}
	
	SpectrumWindow::~SpectrumWindow()
	{
		if (context_menu_!=0) delete(context_menu_);
	}
	
	void SpectrumWindow::setWidget_(SpectrumWidget* widget)
	{
		widget_ = widget;
		widget_->setSpectrumWindow(this);
	}
	
	void SpectrumWindow::resetZoom()
	{
		widget_->canvas()->resetZoom();
	}
	
	void SpectrumWindow::showStatusMessage(std::string msg,OpenMS::UnsignedInt time)
	{
		emit sendStatusMessage(msg,time);
	}
	
	void SpectrumWindow::showCursorStatus(double mz, double intens, double rt)
	{
		emit sendCursorStatus(mz,intens,rt);
	}
	
	SpectrumWidget* SpectrumWindow::widget()
	{
		return widget_;
	}
	
	void SpectrumWindow::modesChangedSlot(QWidget* /*w*/)
	{
		emit modesChanged(this);
	}
	
	void SpectrumWindow::connectWidgetSignals(SpectrumWidget* sw)
	{
		connect(sw,SIGNAL(sendStatusMessage(std::string,OpenMS::UnsignedInt)),this,SLOT(showStatusMessage(std::string,OpenMS::UnsignedInt)));
		connect(sw,SIGNAL(sendCursorStatus(double,double,double)),this,SLOT(showCursorStatus(double,double,double)));
		connect(sw,SIGNAL(modesChanged(QWidget*)),this,SLOT(modesChangedSlot(QWidget*)));
		connect(sw, SIGNAL(contextMenu(QPoint)), this, SLOT(showContextMenu_(QPoint)));
	}
	
	void SpectrumWindow::setActionMode(QAction* a)
	{
		emit sendCursorStatus();
		widget_->setActionMode(a);
	}
	
	int SpectrumWindow::getActionMode()
	{
		return widget_->getActionMode();
	}
	
	bool SpectrumWindow::getGridMode()
	{
		return widget_->canvas()->getGridMode();
	}
	
	void SpectrumWindow::showGridLines(bool b)
	{
		widget_->showGridLines(b);
	}
	
	void SpectrumWindow::createContextMenu_()
	{
		if (context_menu_!=0)
		{
			delete(context_menu_);
		}
		
		//create menu
		context_menu_ = new QPopupMenu(this);
	 	SignedInt item;
	
		//intensity modification
		QPopupMenu* intensity_menu = new QPopupMenu(context_menu_);
	
		item = intensity_menu->insertItem("linear",widget_,SLOT(setIntensityModificationNone()));
		if (widget_->canvas()->getIntensityModification() == SpectrumCanvas::IM_NONE) intensity_menu->setItemEnabled(item,false);
	
		item = intensity_menu->insertItem("logarithmic",widget_,SLOT(setIntensityModificationLog()));
		if (widget_->canvas()->getIntensityModification() == SpectrumCanvas::IM_LOG) intensity_menu->setItemEnabled(item,false);
	
		context_menu_->insertItem("intensity scale",intensity_menu);
		context_menu_->insertSeparator();	
	
		//intensity modification
		context_menu_->insertItem("intensity distribution",widget_,SLOT(showIntensityDistribution()));
		context_menu_->insertSeparator();
	
		//legend menu
		QPopupMenu* legend_menu = new QPopupMenu(context_menu_);
		item = legend_menu->insertItem("on",widget_,SLOT(showLegend()));
		if (widget_->getShowLegend()) legend_menu->setItemEnabled(item,false);
		item = legend_menu->insertItem("off",widget_,SLOT(showNoLegend()));
		if (!widget_->getShowLegend()) legend_menu->setItemEnabled(item,false);
		context_menu_->insertItem("legend",legend_menu);
		context_menu_->insertSeparator();	
	
		//Preferences
		context_menu_->insertItem("Preferences",this,SLOT(showPreferences_()));
		context_menu_->insertSeparator();	
	
	}
	
	void SpectrumWindow::showContextMenu_(QPoint p)
	{
		createContextMenu_();
		context_menu_->popup(p);
	}
	
	void SpectrumWindow::showPreferences_()
	{
		setActive(true);
		emit openPreferences();
	}

}//namespace OpenMS
