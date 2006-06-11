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
// $Id: Spectrum2DWidget.C,v 1.28 2006/06/09 14:46:55 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum2DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum2DCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>

// QT
#include <qimage.h>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	
	Spectrum2DWidget::Spectrum2DWidget(QWidget* parent, const char* name, WFlags f)
		: SpectrumWidget(parent, name, f)
	{
		setCanvas(new Spectrum2DCanvas(this));
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)),
		        this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)),
		        this, SIGNAL(sendCursorStatus(double,double,double)));
	
		x_axis_->setLegend(Spectrum2DCanvas::DimensionDescription::dimension_unit_short[Spectrum2DCanvas::MZ]);
		y_axis_->setLegend(Spectrum2DCanvas::DimensionDescription::dimension_unit_short[Spectrum2DCanvas::RT]);
			
		addClient(canvas(), "Canvas", true);
	
	}
	
	Spectrum2DCanvas* Spectrum2DWidget::canvas() const
	{
		return static_cast<Spectrum2DCanvas*>(canvas_);
	}
	
	Spectrum2DWidget::~Spectrum2DWidget()
	{
	}
	
	void Spectrum2DWidget::showContours(bool on)
	{
		canvas()->showContours(on);
	}
	
	void Spectrum2DWidget::showColors(bool on)
	{
		canvas()->showColors(on);
	}
	
	void Spectrum2DWidget::showPoints(bool on)
	{
		canvas()->showPoints(on);
	}
	
	void Spectrum2DWidget::changeShowContours()
	{
		canvas()->changeShowContours();
	}
	
	void Spectrum2DWidget::changeShowColors()
	{
		canvas()->changeShowColors();
	}
	
	void Spectrum2DWidget::changeShowPoints()
	{
		canvas()->changeShowPoints();
	}
	
	bool Spectrum2DWidget::getShowContours()
	{
		return canvas()->getShowContours();
	}
	
	bool Spectrum2DWidget::getShowColors()
	{
		return canvas()->getShowColors();
	}
	
	bool Spectrum2DWidget::getShowPoints()
	{
		return canvas()->getShowPoints();
	}
	
	void Spectrum2DWidget::recalculateAxes()
	{
		const SpectrumCanvas::AreaType area = canvas()->visible_area_;
		
		if (canvas()->getMappingInfo()->isMzToXAxis())
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
	
	void Spectrum2DWidget::intensityModificationChange_()
	{
		//cout << "IM_CHANGE_WINDOW_2D"<<endl;
		canvas()->intensityModificationChange_();
	}
	
	
	void Spectrum2DWidget::legendModificationChange_()
	{
		y_axis_->showLegend(show_legend_);
		x_axis_->showLegend(show_legend_);
		update();
	}
	
	QImage Spectrum2DWidget::getImage(UnsignedInt width, UnsignedInt height, UnsignedInt flags)
	{
		return canvas()->getImage(width, height, flags);
	}
	
	Histogram<UnsignedInt,float> Spectrum2DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500);
		
		for (Spectrum2DCanvas::ExperimentType::ConstIterator spec_it = canvas()->currentDataSet().begin(); spec_it != canvas()->currentDataSet().end(); ++spec_it)
		{
			for (Spectrum2DCanvas::ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
			{
				tmp.inc(peak_it->getIntensity());
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

