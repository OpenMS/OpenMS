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
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DWidgetPDP.h>


using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum1DWidget::Spectrum1DWidget(QWidget* parent)
		: SpectrumWidget(parent)
	{
		//set the label mode for the axes  - side effect
		setCanvas_(new Spectrum1DCanvas(this));
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)), this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)), this, SIGNAL(sendCursorStatus(double,double,double)));
		
		x_axis_->setLegend("m/z");
		x_axis_->setAllowShortNumbers(false);
		y_axis_->setLegend("Intensity");
		y_axis_->setAllowShortNumbers(true);
		y_axis_->setMinimumWidth(50);
		addClient(canvas(),"Canvas",true);
	}
	
	Spectrum1DCanvas* Spectrum1DWidget::canvas()
	{
		return static_cast<Spectrum1DCanvas*>(canvas_);
	}
	
	void Spectrum1DWidget::recalculateAxes_()
	{
		//determine axes
		AxisWidget* mz_axis,* it_axis;
		if (canvas()->isMzToXAxis())
		{
			mz_axis = x_axis_;
			it_axis = y_axis_;
		}
		else
		{
			mz_axis = y_axis_;
			it_axis = x_axis_;
		}
		
		// recalculate gridlines
		mz_axis->setAxisBounds(canvas()->getVisibleArea().minX(), canvas()->getVisibleArea().maxX());
		switch(canvas()->getIntensityMode())
		{
			case SpectrumCanvas::IM_NONE:
				it_axis->setLogScale(false);
				it_axis->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
				break;
			case SpectrumCanvas::IM_PERCENTAGE:
				it_axis->setLogScale(false);
				it_axis->setAxisBounds(canvas()->getVisibleArea().minY() / canvas()->getDataRange().maxY() * 100.0, canvas()->getVisibleArea().maxY() / canvas()->getDataRange().maxY() * 100.0);
				break;
			case SpectrumCanvas::IM_LOG:
				it_axis->setLogScale(true);
				it_axis->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
				break;
			case SpectrumCanvas::IM_SNAP:
				it_axis->setLogScale(false);
				it_axis->setAxisBounds(canvas()->getVisibleArea().minY()/canvas()->getSnapFactor(), canvas()->getVisibleArea().maxY()/canvas()->getSnapFactor());
				break;
			default:
				throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
	}
	
	Histogram<UnsignedInt,float> Spectrum1DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500.0);
	
		for (Spectrum1DCanvas::ExperimentType::SpectrumType::ConstIterator it = canvas()->getCurrentPeakData()[0].begin(); it != canvas()->getCurrentPeakData()[0].end(); ++it)
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


