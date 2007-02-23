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
// $Maintainer: Cornelia Friedle $
// --------------------------------------------------------------------------

#include<iostream.h>

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum3DWidgetPDP.h>
#include <OpenMS/VISUAL/Spectrum3DCanvas.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>

//QT
#include <QtGui/QPixmap>
#include <QtGui/QGridLayout>
#include <QtGui/QImage>

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum3DWidget::Spectrum3DWidget(QWidget* parent)
	  : SpectrumWidget(parent)		
	{
		setCanvas_(new Spectrum3DCanvas(this));
		
		x_axis_->hide();
		y_axis_->hide();	
		
		connect(canvas(), SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)),this, SIGNAL(sendStatusMessage(std::string, OpenMS::UnsignedInt)));
		
		connect(canvas(), SIGNAL(sendCursorStatus(double,double,double)),
		this, SIGNAL(sendCursorStatus(double,double,double)));
		addClient(canvas(),"Canvas",true);
	}
	
	Spectrum3DWidget::~Spectrum3DWidget()
	{
	
	}
	
	PreferencesDialogPage* Spectrum3DWidget::createPreferences(QWidget* parent)
	{
		PreferencesDialogPage* background = new Spectrum3DWidgetPDP(this, parent);
		return background;
	}
	
	
	void Spectrum3DWidget::recalculateAxes_()
	{
	}
	
	Histogram<UnsignedInt,float> Spectrum3DWidget::createIntensityDistribution_()
	{
		Histogram<UnsignedInt,float> tmp(canvas()->getCurrentMinIntensity(),canvas()->getCurrentMaxIntensity(),(canvas()->getCurrentMaxIntensity() - canvas()->getCurrentMinIntensity())/500.0);

		for (Spectrum3DCanvas::ExperimentType::ConstIterator spec_it = canvas()->getCurrentPeakData().begin(); spec_it != canvas()->getCurrentPeakData().end(); ++spec_it)
		{
			if (spec_it->getMSLevel()!=1)
			{
				continue;
			}
			for (Spectrum3DCanvas::ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
			{
				tmp.inc(peak_it->getIntensity());
			}
		}
		
		return tmp;
	}
	
	Spectrum3DCanvas * Spectrum3DWidget::canvas()
	{
	  return static_cast<Spectrum3DCanvas*>(canvas_);
	}
	
	QImage Spectrum3DWidget::getImage(UnsignedInt width, UnsignedInt height )
	{	
		QPixmap pix = canvas()->openglwidget()->renderPixmap(width,height,true);
		QImage img = pix.toImage();
		return img;
	}
	
	void Spectrum3DWidget::showLegend(int show)
	{
		legend_shown_ = (bool)show;
		canvas()->showLegend(legend_shown_);
	}

	bool Spectrum3DWidget::isLegendShown()  
	{
		return legend_shown_; 
	}
}//namespace

