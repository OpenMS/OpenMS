// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <QtGui/QGridLayout>
#include <QtGui/QScrollBar>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	SpectrumWidget::SpectrumWidget(const Param& /*preferences*/, QWidget* parent)
		: QWidget(parent),
			canvas_(0)
	{
		setAttribute(Qt::WA_DeleteOnClose);
		setMinimumSize(250,250);
		grid_ = new QGridLayout(this);
		grid_->setSpacing(0);
		grid_->setMargin(1);
	}
	
	
	void SpectrumWidget::setCanvas_(SpectrumCanvas* canvas, UInt row, UInt col)
	{
		canvas_ = canvas;
		grid_->addWidget(canvas_, row, col);
		//axes
		y_axis_ = new AxisWidget(AxisWidget::LEFT, "",this);
		x_axis_ = new AxisWidget(AxisWidget::BOTTOM, "",this);
		grid_->addWidget(y_axis_,row,col-1);
		grid_->addWidget(x_axis_,row+1,col);
		connect(canvas_, SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes()));
		connect(canvas_, SIGNAL(recalculateAxes()), this, SLOT(updateAxes()));
		//scrollbars
		x_scrollbar_ = new QScrollBar(Qt::Horizontal, this);
		y_scrollbar_ = new QScrollBar(Qt::Vertical, this);
		y_scrollbar_->setInvertedAppearance(true);
		grid_->addWidget(y_scrollbar_,row,col-2);
		grid_->addWidget(x_scrollbar_,row+2,col);		
		x_scrollbar_->hide();
		y_scrollbar_->hide();
		connect(canvas_, SIGNAL(updateHScrollbar(float,float,float,float)), this, SLOT(updateHScrollbar(float,float,float,float)));
		connect(canvas_, SIGNAL(updateVScrollbar(float,float,float,float)), this, SLOT(updateVScrollbar(float,float,float,float)));
		connect(x_scrollbar_, SIGNAL(valueChanged(int)), canvas_, SLOT(horizontalScrollBarChange(int)));
		connect(y_scrollbar_, SIGNAL(valueChanged(int)), canvas_, SLOT(verticalScrollBarChange(int)));
		connect(canvas_, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)),this, SIGNAL(sendStatusMessage(std::string, OpenMS::UInt)));
		connect(canvas_, SIGNAL(sendCursorStatus(double,double,double)), this, SIGNAL(sendCursorStatus(double,double,double)));
		
		//swap axes if necessary
		updateAxes();
		
		canvas_->setSpectrumWidget(this);
	}
	
	SpectrumWidget::~SpectrumWidget()
	{
		//cout << "DEST SpectrumWidget" << endl;

		emit aboutToBeDestroyed(window_id);
	}
	
	Int SpectrumWidget::getActionMode() const 
	{
		return canvas_->getActionMode(); 
	}
	
	void SpectrumWidget::setActionMode(SpectrumCanvas::ActionModes mode)
	{
		if (getActionMode() != mode)
		{
			canvas_->setActionMode(mode);
			emit modesChanged(this);
		}
	}
	
	void SpectrumWidget::setIntensityMode(SpectrumCanvas::IntensityModes mode)
	{
		if (canvas_->getIntensityMode() != mode)
		{
			canvas_->setIntensityMode(mode);
			intensityModeChange_();
		}
	}
	
	void SpectrumWidget::showStatistics()
	{
		LayerStatisticsDialog lsd(this);
		lsd.exec();
	}
	
	void SpectrumWidget::showIntensityDistribution()
	{
		Histogram<UInt,float> dist = createIntensityDistribution_();
		HistogramDialog dw(dist);
		dw.setLegend("intensity");
		
		if (dw.exec() == QDialog::Accepted)
		{
			DataFilters filters;
			
			if (dw.getLeftSplitter()>dist.min())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getLeftSplitter();
				filter.field = DataFilters::INTENSITY;
				filter.op = DataFilters::GREATER_EQUAL;
				filters.add(filter);
			}
			
			if (dw.getRightSplitter()<dist.max())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getRightSplitter();
				filter.field = DataFilters::INTENSITY;
				filter.op = DataFilters::LESS_EQUAL;
				filters.add(filter);
			}
			
			canvas_->setFilters(filters);
			emit sendStatusMessage("Displayed intensity range: "+String(dw.getLeftSplitter())+" upto "+String(dw.getRightSplitter())+" m/z", 5000);
		}
	}
	
	void SpectrumWidget::showLegend(bool show)
	{
		y_axis_->showLegend(show);
		x_axis_->showLegend(show);
		update();
	}
	
	void SpectrumWidget::updateAxes()
	{
		//change axis lables if necessary
		if (canvas()->isMzToXAxis()==true && x_axis_->getLegend().size()>=2 && x_axis_->getLegend().prefix(2)=="RT"
		|| canvas()->isMzToXAxis()==false && y_axis_->getLegend().size()>=2 && y_axis_->getLegend().prefix(2)=="RT") 
		{
			std::string tmp = x_axis_->getLegend();
			x_axis_->setLegend(y_axis_->getLegend());
			y_axis_->setLegend(tmp);
		}
		recalculateAxes_();
	}

	void SpectrumWidget::intensityModeChange_()
	{
		
	}

	QImage SpectrumWidget::getImage(UInt width, UInt height)
	{		
		//hide scrollbars if necessary
		bool x_sc_on = x_scrollbar_->isVisible();
		bool y_sc_on = y_scrollbar_->isVisible();
		if (x_sc_on)
		{
			x_scrollbar_->hide();
		}
		if (y_sc_on)
		{
			y_scrollbar_->hide();
		}
		
		//store old background color
	 	QColor old_bg_color = palette().window().color();
	 	//set bg color to white
		QPalette new_palette;
		new_palette.setColor(backgroundRole(), Qt::white);
		setPalette(new_palette);

		//set white background
		y_axis_->setPalette(palette());
		x_axis_->setPalette(palette());
	
		//set pen width
		int pen = width/1024;
		y_axis_->setPenWidth(pen);
		x_axis_->setPenWidth(pen);
		
	 	//store size and resize
		int h = this->height();
		int w = this->width();

		grid_->activate();
		resize(width,height);		
				
		//take an image
		QImage image = QPixmap::grabWidget(this).toImage();
		
		//resore background colors
		new_palette.setColor(backgroundRole(), old_bg_color);
		setPalette(new_palette);
		y_axis_->setPalette(palette());
		x_axis_->setPalette(palette());
		
		//restore pen withs
		y_axis_->setPenWidth(0);
		x_axis_->setPenWidth(0);
		
		//show scrollbars again
		if (x_sc_on)
		{
			x_scrollbar_->show();
		}
		if (y_sc_on)
		{
			y_scrollbar_->show();
		}
		
		//restore size
		resize(w,h);
		
		return image;
	}

	bool SpectrumWidget::isLegendShown() const 
	{
		//Both are shown or hidden, so we simply return the label of the x-axis
		return x_axis_->isLegendShown(); 
	}

	void SpectrumWidget::hideAxes()
	{
		y_axis_->hide();
		x_axis_->hide();
	}

	void SpectrumWidget::updateHScrollbar(float min, float disp_min, float disp_max, float max)
	{
		if (min == disp_min && max == disp_max)
		{
			x_scrollbar_->hide();
		}
		else
		{
			//block signals as this causes repainting due to rounding (QScrollBar works with int ...)
			x_scrollbar_->blockSignals(true);
			x_scrollbar_->show();
			x_scrollbar_->setMinimum(static_cast<int>(min));
			x_scrollbar_->setMaximum(static_cast<int>(max-disp_max+disp_min));
			x_scrollbar_->setValue(static_cast<int>(disp_min));
			x_scrollbar_->setPageStep(static_cast<int>(disp_max-disp_min));
			x_scrollbar_->blockSignals(false);
		}
	}

	void SpectrumWidget::updateVScrollbar(float min, float disp_min, float disp_max, float max)
	{
		if (min == disp_min && max == disp_max)
		{
			y_scrollbar_->hide();
		}
		else
		{
			//block signals as this causes repainting due to rounding (QScrollBar works with int ...)
			y_scrollbar_->blockSignals(true);
			y_scrollbar_->show();
			y_scrollbar_->setMinimum(static_cast<int>(min));
			y_scrollbar_->setMaximum(static_cast<int>(max-disp_max+disp_min));
			y_scrollbar_->setValue(static_cast<int>(disp_min));
			y_scrollbar_->setPageStep(static_cast<int>(disp_max-disp_min));
			y_scrollbar_->blockSignals(false);
		}
	}

} //namespace OpenMS

