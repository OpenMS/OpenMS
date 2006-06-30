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
// $Id: Spectrum2DWindow.C,v 1.22 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// QT
#include <qpopupmenu.h>
#include <qlayout.h>

// STL

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DWindow.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

namespace OpenMS
{

	Spectrum2DWindow::Spectrum2DWindow(QWidget* parent, const char* name, WFlags f)  
		: SpectrumWindow(parent,name,f)
	{
		QWidget* w = new QWidget(this,"centralWidget");
		setCentralWidget(w);
		grid_ = new QGridLayout(w, 2, 2, 0, 0, "Spectrum2DGridLayout");
		grid_->setRowStretch(1, 3);
		grid_->setColStretch(1, 3);
		
		tic_ = new 	Spectrum1DWidget(w);
		tic_->hideAxes();
		tic_->mzToXAxis(true);
		tic_->hide();
	
		projection_ = new Spectrum1DWidget(w);
		projection_->hideAxes();
		projection_->hide();
		
		setWidget_(new Spectrum2DWidget(w));
		
		grid_->addWidget(projection_, 0, 1);
		grid_->addWidget(tic_, 1, 0);
		grid_->addWidget(widget(), 1, 1);
		connectWidgetSignals(widget());
		
		show1DProjections(false);
		
		connect(widget(), SIGNAL(contextMenu(QPoint)), this, SLOT(showContextMenu_(QPoint)));
		
		connect(widget()->canvas(), SIGNAL(selectedHorz(const DSpectrum<1>&)), this, SLOT(horizontalSpectrum(const DSpectrum<1>&)));
		connect(widget()->canvas(), SIGNAL(selectedVert(const DSpectrum<1>&)), this, SLOT(horizontalSpectrum(const DSpectrum<1>&)));
	
	}
	
	Spectrum2DWidget* Spectrum2DWindow::widget()
	{
		return static_cast<Spectrum2DWidget*>(widget_);
	}
	
	Spectrum2DWindow::~Spectrum2DWindow()
	{
		
	}
	
	void Spectrum2DWindow::createContextMenu_()
	{
		SpectrumWindow::createContextMenu_();
		
		SignedInt item;
	 	
		QPopupMenu* proj_menu = new QPopupMenu(context_menu_);
		item = proj_menu->insertItem("on",this,SLOT(changeShow1DProjections()));
		if (projection_->isVisible() || tic_->isVisible()) proj_menu->setItemEnabled(item,false);
		item = proj_menu->insertItem("off",this,SLOT(changeShow1DProjections()));
		if (!(projection_->isVisible() || tic_->isVisible())) proj_menu->setItemEnabled(item,false);
		context_menu_->insertItem("1D projections",proj_menu);
	}
	
	void Spectrum2DWindow::show1DProjections(bool on)
	{
	// 	show_1D_projections_ = on;
		if (!on) {
			projection_->hide();
			tic_->hide();
	/*		grid_->remove(projection_);
			grid_->remove(tic_);
			grid_->setRowStretch(0,0);
			grid_->setColStretch(0,0);*/
		} else {
	/*		grid_->setRowStretch(0,1);
			grid_->setColStretch(0,1);
			grid_->addWidget(projection_,0,1);
			grid_->addWidget(tic_,1,0);*/
			projection_->show();
			tic_->show();
		}
	}
	
	void Spectrum2DWindow::changeShow1DProjections()
	{
		show1DProjections(!(projection_->isVisible() || tic_->isVisible()));
	}
	
	void Spectrum2DWindow::horizontalSpectrum(const DSpectrum<1>&)
	{
		projection_->show();
	}
	
	void Spectrum2DWindow::verticalSpectrum(const DSpectrum<1>&)
	{
		tic_->show();
	}
	
	PreferencesDialogPage* Spectrum2DWindow::createPreferences(QWidget* parent)
	{
		return widget()->createPreferences(parent);
	}
	
	SignedInt Spectrum2DWindow::getDotMode()
	{
		return widget()->canvas()->getDotMode();
	}
	
	void Spectrum2DWindow::showGoToDialog()
	{
	  Spectrum2DGoToDialog goToDialog(this, "Spectrum2DGoToDialog");
	  const DRange<3>& area = widget()->canvas()->getDataRange();
	  goToDialog.setMinX(area.minX());
	  goToDialog.setMaxX(area.maxX());
	  goToDialog.setMinY(area.minY());
	  goToDialog.setMaxY(area.maxY());  
	  goToDialog.exec();
	  widget()->canvas()->setVisibleArea(SpectrumCanvas::AreaType( goToDialog.getMinX(), goToDialog.getMinY(), goToDialog.getMaxX(), goToDialog.getMaxY()));
	}

}//namespace OpenMS


