// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream>

//OpenMS
#include <OpenMS/VISUAL/Spectrum3DWidget.h>
#include <OpenMS/VISUAL/Spectrum3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum2DGoToDialog.h>

//QT
#include <QtGui/QGridLayout>

using namespace std;

namespace OpenMS
{
	using namespace Internal;
	using namespace Math;
	
	Spectrum3DWidget::Spectrum3DWidget(const Param& preferences, QWidget* parent)
	  : SpectrumWidget(preferences, parent)		
	{
		setCanvas_(new Spectrum3DCanvas(preferences, this));
		
		x_axis_->hide();
		y_axis_->hide();
	}
	
	Spectrum3DWidget::~Spectrum3DWidget()
	{
	
	}
	
	void Spectrum3DWidget::recalculateAxes_()
	{
	}
	
	Histogram<> Spectrum3DWidget::createIntensityDistribution_() const
	{
		//initialize histogram
		DoubleReal min = canvas_->getCurrentMinIntensity();
		DoubleReal max = canvas_->getCurrentMaxIntensity();
		if (min==max)
		{
			min-=0.01;
			max+=0.01;
		}
		Histogram<> tmp(min,max,(max-min)/500.0);
		
		for (ExperimentType::ConstIterator spec_it = canvas_->getCurrentLayer().peaks.begin(); spec_it != canvas_->getCurrentLayer().peaks.end(); ++spec_it)
		{
			if (spec_it->getMSLevel()!=1) continue;
			for (ExperimentType::SpectrumType::ConstIterator peak_it = spec_it->begin(); peak_it != spec_it->end(); ++peak_it)
			{
				tmp.inc(peak_it->getIntensity());
			}
		}
		
		return tmp;
	}

	Histogram<> Spectrum3DWidget::createMetaDistribution_(const String& name) const
	{
		Histogram<> tmp;
		
		//determine min and max of the data
		Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
		for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().peaks.begin(); s_it!=canvas_->getCurrentLayer().peaks.end(); ++s_it)
		{
			if (s_it->getMSLevel()!=1) continue;
			//float arrays
			for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it=s_it->getFloatDataArrays().begin(); it!=s_it->getFloatDataArrays().end(); it++)
			{
				if (it->getName()==name)
				{
					for (Size i=0; i<it->size(); ++i)
					{
						if ((*it)[i]<min) min = (*it)[i];
						if ((*it)[i]>max) max = (*it)[i];
					}
					break;
				}
			}
			//integer arrays
			for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it=s_it->getIntegerDataArrays().begin(); it!=s_it->getIntegerDataArrays().end(); it++)
			{
				if (it->getName()==name)
				{
					for (Size i=0; i<it->size(); ++i)
					{
						if ((*it)[i]<min) min = (*it)[i];
						if ((*it)[i]>max) max = (*it)[i];
					}
					break;
				}
			}
		}
		if (min>=max) return tmp;
		
		//create histogram
		tmp.reset(min,max,(max-min)/500.0);
		for (ExperimentType::const_iterator s_it = canvas_->getCurrentLayer().peaks.begin(); s_it!=canvas_->getCurrentLayer().peaks.end(); ++s_it)
		{
			if (s_it->getMSLevel()!=1) continue;
			//float arrays
			for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it=s_it->getFloatDataArrays().begin(); it!=s_it->getFloatDataArrays().end(); it++)
			{
				if (it->getName()==name)
				{
					for (Size i=0; i<it->size(); ++i)
					{
						tmp.inc((*it)[i]);
					}
					break;
				}
			}
			//integer arrays
			for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it=s_it->getIntegerDataArrays().begin(); it!=s_it->getIntegerDataArrays().end(); it++)
			{
				if (it->getName()==name)
				{
					for (Size i=0; i<it->size(); ++i)
					{
						tmp.inc((*it)[i]);
					}
					break;
				}
			}
		}
		
		return tmp;
	}	

	void Spectrum3DWidget::showLegend(bool show)
	{
		canvas()->showLegend(show);
	}

	bool Spectrum3DWidget::isLegendShown() const
	{
		return static_cast<const Spectrum3DCanvas*>(canvas_)->isLegendShown();
	}

	void Spectrum3DWidget::showGoToDialog()
	{
		Spectrum2DGoToDialog goto_dialog(this);
		const DRange<3>& area = canvas()->getDataRange();
		goto_dialog.setRange(area.minY(),area.maxY(),area.minX(),area.maxX());
		goto_dialog.enableFeatureNumber(false);
		if(goto_dialog.exec())
		{
			canvas()->setVisibleArea(SpectrumCanvas::AreaType( goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT()));
		}
	}

}//namespace

