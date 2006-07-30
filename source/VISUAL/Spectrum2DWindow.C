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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// QT
#include <qpopupmenu.h>
#include <qlayout.h>

// STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum2DWindow.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

using namespace std;

namespace OpenMS
{

	Spectrum2DWindow::Spectrum2DWindow(QWidget* parent, const char* name, WFlags f)  
		: SpectrumWindow(parent,name,f)
	{
		QWidget* w = new QWidget(this,"centralWidget");
		setCentralWidget(w);

		grid_ = new QGridLayout(w, 2, 2, 0, 0, "Spectrum2DGridLayout");
		grid_->setRowStretch(1, 3);
		grid_->setColStretch(0, 3);
		
		projection_vert_ = new 	Spectrum1DWidget(w);
		projection_vert_->hide();
	
		projection_horz_ = new Spectrum1DWidget(w);
		projection_horz_->hide();
		
		setWidget_(new Spectrum2DWidget(w));
		
		grid_->addWidget(projection_horz_, 0, 0);
		grid_->addWidget(projection_vert_, 1, 1);
		grid_->addWidget(widget(), 1, 0);
		
		connectWidgetSignals(widget());
		
		connect(widget(), SIGNAL(contextMenu(QPoint)), this, SLOT(showContextMenu_(QPoint)));
		connect(widget()->canvas(), SIGNAL(showProjectionHorizontal(const MSExperiment<>&)), this, SLOT(horizontalProjection(const MSExperiment<>&)));
		connect(widget()->canvas(), SIGNAL(showProjectionVertical(const MSExperiment<>&)), this, SLOT(verticalProjection(const MSExperiment<>&)));
	
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
		context_menu_->insertItem("Hide Projections",this,SLOT(hideProjections()));
	}
	
	void Spectrum2DWindow::hideProjections()
	{
		projection_horz_->hide();
		projection_vert_->hide();
	}
	
	void Spectrum2DWindow::horizontalProjection(const MSExperiment<>& exp)
	{
		if (exp[0].size()<3)
		{
			projection_horz_->hide();
			return;
		}
		projection_horz_->setMainPreferences(prefs_);
		projection_horz_->mzToXAxis(true);
		projection_horz_->showLegend(false);
		projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_horz_->canvas()->removeDataSet(0);
		projection_horz_->canvas()->addDataSet(exp);
		projection_horz_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_horz_->show();
	}
	
	void Spectrum2DWindow::verticalProjection(const MSExperiment<>& exp)
	{
		if (exp[0].size()<3)
		{
			projection_vert_->hide();
			return;
		}
		projection_vert_->setMainPreferences(prefs_);
		projection_vert_->mzToXAxis(false);
		projection_vert_->showLegend(false);
		projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_vert_->canvas()->removeDataSet(0);
		projection_vert_->canvas()->addDataSet(exp);
		projection_vert_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_vert_->show();
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

	const Spectrum1DWidget* Spectrum2DWindow::getHorizontalProjection() const
	{
		return projection_horz_;
	}

	const Spectrum1DWidget* Spectrum2DWindow::getVerticalProjection() const
	{
		return projection_vert_;
	}

}//namespace OpenMS


