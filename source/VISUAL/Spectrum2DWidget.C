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

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

#include <QtGui/QPushButton>
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum2DWidget::Spectrum2DWidget(QWidget* parent)
		: SpectrumWidget(parent)
	{
		setCanvas_(new Spectrum2DCanvas(this),1,2);

		x_axis_->setLegend(String(RawDataPoint2D::shortDimensionName(RawDataPoint2D::MZ))+" ["+String(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::MZ))+"]");
		y_axis_->setLegend(String(RawDataPoint2D::shortDimensionName(RawDataPoint2D::RT))+" ["+String(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::RT))+"]");
		y_axis_->setMinimumWidth(50);
		
		addClient(canvas(), "Canvas", true);
		
		projection_vert_ = new 	Spectrum1DWidget(this);
		projection_vert_->hide();
		grid_->addWidget(projection_vert_,1,3,2,1);
		
		projection_horz_ = new Spectrum1DWidget(this);
		projection_horz_->hide();
		grid_->addWidget(projection_horz_,0,1,1,2);
		connect(canvas(), SIGNAL(showProjectionHorizontal(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)), this, SLOT(horizontalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)));
		connect(canvas(), SIGNAL(showProjectionVertical(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)), this, SLOT(verticalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)));
		
		hide_button_ = new QPushButton("Hide projections", this);
		hide_button_->hide();
		grid_->addWidget(hide_button_, 0, 3, Qt::AlignLeft | Qt::AlignBottom);
		connect(hide_button_, SIGNAL(clicked()), this, SLOT(hideProjections()));
	}
	
	Spectrum2DCanvas* Spectrum2DWidget::canvas()
	{
		return static_cast<Spectrum2DCanvas*>(canvas_);
	}
	
	Spectrum2DWidget::~Spectrum2DWidget()
	{
		
	}
	
	void Spectrum2DWidget::recalculateAxes_()
	{
		const SpectrumCanvas::AreaType area = canvas()->getVisibleArea();
		
		if (canvas()->isMzToXAxis())
		{
			x_axis_->setAxisBounds(area.minX(), area.maxX());
			y_axis_->setAxisBounds(area.minY(), area.maxY());
		}
		else
		{
			x_axis_->setAxisBounds(area.minY(), area.maxY());
			y_axis_->setAxisBounds(area.minX(), area.maxX());
		}
	}
	
	Histogram<UInt,float> Spectrum2DWidget::createIntensityDistribution_()
	{
		Histogram<UInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500.0);
		
		if (canvas()->getCurrentLayer().type==LayerData::DT_PEAK)
		{
			for (Spectrum2DCanvas::ExperimentType::ConstIterator spec_it = canvas()->getCurrentPeakData().begin(); spec_it != canvas()->getCurrentPeakData().end(); ++spec_it)
			{
				if (spec_it->getMSLevel()!=1)
				{
					continue;
				}
				for (Spectrum2DCanvas::ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
				{
					tmp.inc(peak_it->getIntensity());
				}
			}
		}
		else
		{
			for (Spectrum2DCanvas::FeatureMapType::ConstIterator it = canvas()->getCurrentLayer().features.begin(); it != canvas()->getCurrentLayer().features.end(); ++it)
			{
				tmp.inc(it->getIntensity());
			}
		}
		
		return tmp;
	}
	
	PreferencesDialogPage* Spectrum2DWidget::createPreferences(QWidget* parent)
	{
		PreferencesDialogPage* background = new Spectrum2DWidgetPDP(this,parent);
		return background;
	}

	void Spectrum2DWidget::hideProjections()
	{
		hide_button_->hide();
		projection_horz_->hide();
		projection_vert_->hide();
	}
	
	void Spectrum2DWidget::horizontalProjection(const MSExperiment<>& exp, Spectrum1DCanvas::DrawModes mode)
	{
		projection_horz_->setMainPreferences(prefs_);
		projection_horz_->mzToXAxis(true);
		projection_horz_->showLegend(false);
		projection_horz_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_horz_->canvas()->removeLayer(0);
		projection_horz_->canvas()->addLayer(exp);
		projection_horz_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_horz_->canvas()->setDrawMode(mode);
		projection_horz_->show();
		hide_button_->show();
	}
	
	void Spectrum2DWidget::verticalProjection(const MSExperiment<>& exp, Spectrum1DCanvas::DrawModes mode)
	{
		projection_vert_->setMainPreferences(prefs_);
		projection_vert_->mzToXAxis(false);
		projection_vert_->showLegend(false);
		projection_vert_->canvas()->setIntensityMode(SpectrumCanvas::IM_PERCENTAGE);
		projection_vert_->canvas()->removeLayer(0);
		projection_vert_->canvas()->addLayer(exp);
		projection_vert_->canvas()->setActionMode(SpectrumCanvas::AM_SELECT);
		projection_vert_->canvas()->setDrawMode(mode);
		projection_vert_->show();
		hide_button_->show();
	}

	const Spectrum1DWidget* Spectrum2DWidget::getHorizontalProjection() const
	{
		return projection_horz_;
	}

	const Spectrum1DWidget* Spectrum2DWidget::getVerticalProjection() const
	{
		return projection_vert_;
	}

	void Spectrum2DWidget::showGoToDialog()
	{
	  Spectrum2DGoToDialog goto_dialog(this);
	  const DRange<3>& area = canvas()->getDataRange();
	  goto_dialog.setMinRT(area.minY());
	  goto_dialog.setMaxRT(area.maxY());
	  goto_dialog.setMinMZ(area.minX());
	  goto_dialog.setMaxMZ(area.maxX());  
	  if(goto_dialog.exec())
	  {
	  	canvas()->setVisibleArea(SpectrumCanvas::AreaType( goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT()));
		}
	}

} //Namespace

