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

#include <iostream.h>

#include <OpenMS/VISUAL/SpectrumWindow.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>

namespace OpenMS
{

	SpectrumWindow::SpectrumWindow(QWidget* parent)  
		: QMainWindow(parent),
			PreferencesManager(),
			window_id(-1)
	{
		setAttribute(Qt::WA_DeleteOnClose);
		setMinimumSize(300,300);	// prevents errors caused by too small width,height values
	}
	
	SpectrumWindow::~SpectrumWindow()
	{
		emit aboutToBeDestroyed(window_id);
	}
	
	void SpectrumWindow::setWidget_(SpectrumWidget* widget)
	{
		widget_ = widget;
		widget_->setSpectrumWindow(this);
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
	}

	void SpectrumWindow::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
		widget()->setMainPreferences(prefs);
	}

}//namespace OpenMS
