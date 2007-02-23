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

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum2DWidget::Spectrum2DWidget(QWidget* parent)
		: SpectrumWidget(parent)
	{
		setCanvas_(new Spectrum2DCanvas(this));
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)),
		        this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)),
		        this, SIGNAL(sendCursorStatus(double,double,double)));
	
		x_axis_->setLegend(DimensionDescription < LCMS_Tag >::dimension_unit_short[Spectrum2DCanvas::MZ]);
		y_axis_->setLegend(DimensionDescription < LCMS_Tag >::dimension_unit_short[Spectrum2DCanvas::RT]);
		y_axis_->setMinimumWidth(50);
		
		addClient(canvas(), "Canvas", true);
	
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
	
	Histogram<UnsignedInt,float> Spectrum2DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500.0);
		
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

} //Namespace

