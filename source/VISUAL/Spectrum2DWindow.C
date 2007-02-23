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

// QT
#include <QtGui/QGridLayout>

// STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DWindow.h>
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

using namespace std;

namespace OpenMS
{

	Spectrum2DWindow::Spectrum2DWindow(QWidget* parent)  
		: SpectrumWindow(parent)
	{
		QWidget* w = new QWidget(this);
		setCentralWidget(w);

		grid_ = new QGridLayout(w);
		grid_->setRowStretch(1, 3);
		grid_->setColumnStretch(0, 3);
		
		projection_vert_ = new 	Spectrum1DWidget(w);
		projection_vert_->hide();
		grid_->addWidget(projection_vert_, 1, 1);
		
		projection_horz_ = new Spectrum1DWidget(w);
		projection_horz_->hide();
		grid_->addWidget(projection_horz_, 0, 0);
		
		setWidget_(new Spectrum2DWidget(w));
		grid_->addWidget(widget(), 1, 0);
		connectWidgetSignals(widget());
		connect(widget()->canvas(), SIGNAL(showProjectionHorizontal(const MSExperiment<>&)), this, SLOT(horizontalProjection(const MSExperiment<>&)));
		connect(widget()->canvas(), SIGNAL(showProjectionVertical(const MSExperiment<>&)), this, SLOT(verticalProjection(const MSExperiment<>&)));
		
		hide_button_ = new QPushButton("Hide projections", w);
		hide_button_->hide();
		grid_->addWidget(hide_button_, 0, 1, Qt::AlignLeft | Qt::AlignBottom);
		connect(hide_button_, SIGNAL(clicked()), this, SLOT(hideProjections()));
	}
	
	Spectrum2DWidget* Spectrum2DWindow::widget()
	{
		return static_cast<Spectrum2DWidget*>(widget_);
	}
	
	Spectrum2DWindow::~Spectrum2DWindow()
	{
		
	}
	
	void Spectrum2DWindow::hideProjections()
	{
		hide_button_->hide();
		projection_horz_->hide();
		projection_vert_->hide();
	}
	
	void Spectrum2DWindow::horizontalProjection(const MSExperiment<>& exp)
	{
		projection_horz_->setMainPreferences(prefs_);
		projection_horz_->mzToXAxis(true);
		projection_horz_->showLegend(false);
		projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_horz_->canvas()->removeLayer(0);
		projection_horz_->canvas()->addLayer(exp);
		projection_horz_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_horz_->show();
		hide_button_->show();
	}
	
	void Spectrum2DWindow::verticalProjection(const MSExperiment<>& exp)
	{
		projection_vert_->setMainPreferences(prefs_);
		projection_vert_->mzToXAxis(false);
		projection_vert_->showLegend(false);
		projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_vert_->canvas()->removeLayer(0);
		projection_vert_->canvas()->addLayer(exp);
		projection_vert_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_vert_->show();
		hide_button_->show();
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
	  Spectrum2DGoToDialog goto_dialog(this);
	  const DRange<3>& area = widget()->canvas()->getDataRange();
	  goto_dialog.setMinRT(area.minY());
	  goto_dialog.setMaxRT(area.maxY());
	  goto_dialog.setMinMZ(area.minX());
	  goto_dialog.setMaxMZ(area.maxX());  
	  if(goto_dialog.exec())
	  {
	  	widget()->canvas()->setVisibleArea(SpectrumCanvas::AreaType( goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT()));
		}
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


