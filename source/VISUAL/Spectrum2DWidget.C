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
#include <QtGui/QGroupBox>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum2DWidget::Spectrum2DWidget(QWidget* parent)
		: SpectrumWidget(parent)
	{
		setCanvas_(new Spectrum2DCanvas(this),1,2);
		
		// add axes
		x_axis_->setLegend(String(RawDataPoint2D::shortDimensionName(RawDataPoint2D::MZ))+" ["+String(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::MZ))+"]");
		y_axis_->setLegend(String(RawDataPoint2D::shortDimensionName(RawDataPoint2D::RT))+" ["+String(RawDataPoint2D::shortDimensionUnit(RawDataPoint2D::RT))+"]");
		y_axis_->setMinimumWidth(50);
		
		addClient(canvas(), "Canvas", true);
		
		// add projetions
		grid_->setColumnStretch(2,3);
		grid_->setRowStretch(1,3);
		
		projection_vert_ = new 	Spectrum1DWidget(this);
		projection_vert_->hide();
		grid_->addWidget(projection_vert_,1,3,2,1);
		
		projection_horz_ = new Spectrum1DWidget(this);
		projection_horz_->hide();
		grid_->addWidget(projection_horz_,0,1,1,2);
		connect(canvas(), SIGNAL(showProjectionHorizontal(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)), this, SLOT(horizontalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)));
		connect(canvas(), SIGNAL(showProjectionVertical(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)), this, SLOT(verticalProjection(const MSExperiment<>&, Spectrum1DCanvas::DrawModes)));
		connect(canvas(), SIGNAL(showProjectionInfo(int,double)), this, SLOT(projectionInfo(int,double)));
		connect(canvas(), SIGNAL(showCurrentPeaksAs3D()), this, SIGNAL(showCurrentPeaksAs3D()));
		connect(canvas(), SIGNAL(showSpectrumAs1D(int)), this, SIGNAL(showSpectrumAs1D(int)));
		
		// add projections box
		projection_box_ = new QGroupBox("Projections",this);
		projection_box_->hide();
		grid_->addWidget(projection_box_, 0, 3);
		QGridLayout* box_grid = new QGridLayout(projection_box_); 
		
		QLabel* label = new QLabel("Peaks: ");
		box_grid->addWidget(label,0,0);
		projection_peaks_ = new QLabel("");
		box_grid->addWidget(projection_peaks_,0,1,1,2);
		
		label = new QLabel("Intensity: ");
		box_grid->addWidget(label,1,0);
		projection_sum_ = new QLabel("");
		box_grid->addWidget(projection_sum_,1,1,1,2);
		
		QPushButton* button = new QPushButton("Hide", projection_box_);
		connect(button, SIGNAL(clicked()), this, SLOT(hideProjections()));
		box_grid->addWidget(button,3,1);

		button = new QPushButton("Update", projection_box_);
		connect(button, SIGNAL(clicked()), canvas(), SLOT(showProjections()));
		box_grid->addWidget(button,3,2);
		
		box_grid->setRowStretch(2,2);
	}

	Spectrum2DWidget::~Spectrum2DWidget()
	{
		//cout << "DEST Spectrum2DWidget" << endl;
	}

	void Spectrum2DWidget::projectionInfo(int peaks, double intensity)
	{
		projection_peaks_->setText(QString::number(peaks));
		projection_sum_->setText(QString::number(intensity,'f',1));
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
		projection_box_->hide();
		projection_horz_->hide();
		projection_vert_->hide();
		grid_->setColumnStretch(3,0);
		grid_->setRowStretch(0,0);
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
		grid_->setColumnStretch(3,2);
		projection_horz_->show();
		projection_box_->show();
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
		grid_->setRowStretch(0,2);
		projection_vert_->show();
		projection_box_->show();
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

