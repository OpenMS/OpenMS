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
// $Id: Spectrum1DWidget.C,v 1.36 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// Qt
#include <qaction.h>

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DWidgetPDP.h>


using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	Spectrum1DWidget::Spectrum1DWidget(QWidget* parent, const char* name, WFlags f)
		: SpectrumWidget(parent, name, f)
	{
		//set the label mode for the axes  - side effect
		setCanvas(new Spectrum1DCanvas(this, "Spectrum1DCanvas"));
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)), this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)), this, SIGNAL(sendCursorStatus(double,double,double)));
		
		recalculateAxes();
		
		x_axis_->setLegend("m/z");
		y_axis_->setLegend("Intensity");
		addClient(canvas(),"Canvas",true);
	
		setMouseTracking(true);
	}
	
	Spectrum1DCanvas* Spectrum1DWidget::canvas() const
	{
		return static_cast<Spectrum1DCanvas*>(canvas_);
	}
	
	void Spectrum1DWidget::mouseMoveEvent( QMouseEvent* /*e*/)
	{
		
	}
	
	void Spectrum1DWidget::setVisibleArea(double x1, double x2)
	{
		canvas()->setVisibleArea(x1, x2);
	}
	
	void Spectrum1DWidget::recalculateAxes()
	{
		//set intensity axis to log scale if necessary
		if (canvas()->isMzToXAxis())
		{
			// y-axis is intensity axis
			x_axis_->setLogScale(false);
			y_axis_->setLogScale(canvas()->getIntensityMode() == SpectrumCanvas::IM_LOG);
		}
		else
		{
			// x-axis is intensity axis
			x_axis_->setLogScale(canvas()->getIntensityMode() == SpectrumCanvas::IM_LOG);
			y_axis_->setLogScale(false);
		}
		
		const SpectrumCanvas::AreaType& visible_area = canvas()->visible_area_;
		
		SpectrumCanvas::AreaType data_area;
		data_area.assign(canvas()->getDataRange());
	
		// recalculate gridlines
		double lx,hx,ly,hy;
	
		if (canvas()->isMzToXAxis())
		{
			lx = visible_area.minX();
			hx = visible_area.maxX();
			ly = visible_area.minY();
			hy = visible_area.maxY();
		}
		else
		{
			lx = visible_area.minY();
			hx = visible_area.maxY();
			ly = visible_area.minX();
			hy = visible_area.maxX();
		}

		if (canvas()->isMzToXAxis())  // y = intensity
		{
			x_axis_->setAxisBounds(lx, hx);
			if (canvas()->getIntensityMode() == SpectrumCanvas::IM_PERCENTAGE) // Adjust axis values in snap-to-max-intensity-mode
			{
				y_axis_->setAxisBounds(ly/canvas()->getSnapFactor(), hy/canvas()->getSnapFactor());
			}
			else
			{
				y_axis_->setAxisBounds(ly, hy);
			}
		}
		else  // x = intensity
		{
			y_axis_->setAxisBounds(ly, hy);
			if (canvas()->getIntensityMode() == SpectrumCanvas::IM_PERCENTAGE) // Adjust axis values in snap-to-max-intensity-mode
			{
				x_axis_->setAxisBounds(lx/canvas()->getSnapFactor(), hx/canvas()->getSnapFactor());
			}
			else
			{
				x_axis_->setAxisBounds(lx, hx);
			}
		}
	}
	
	void Spectrum1DWidget::intensityModeChange_()
	{
		//recalculate axes before and after (for grid update)
		recalculateAxes();
	 	canvas()->intensityModeChange_();
		recalculateAxes();
	}
	
	Histogram<UnsignedInt,float> Spectrum1DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500.0);
	
		for (Spectrum1DCanvas::ExperimentType::SpectrumType::ConstIterator it = canvas()->currentDataSet()[0].begin(); it != canvas()->currentDataSet()[0].end(); ++it)
		{
			tmp.inc(it->getIntensity());
		}
		return tmp;
	}
	
	// destructor
	Spectrum1DWidget::~Spectrum1DWidget()
	{
		
	}
	
	PreferencesDialogPage* Spectrum1DWidget::createPreferences(QWidget* parent)
	{
		PreferencesDialogPage* background = new Spectrum1DWidgetPDP(this, parent);
		return background;
	}

} //namespace


