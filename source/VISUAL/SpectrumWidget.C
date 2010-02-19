// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <QtGui/QGridLayout>
#include <QtGui/QScrollBar>
#include <QtGui/QCloseEvent>
#include <QtGui/QMessageBox>
#include <QtGui/QFileDialog>
#include <QtCore/QMimeData>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	SpectrumWidget::SpectrumWidget(const Param& /*preferences*/, QWidget* parent)
		: QWidget(parent),
			canvas_(0)
	{
		setAttribute(Qt::WA_DeleteOnClose);
		grid_ = new QGridLayout(this);
		grid_->setSpacing(0);
		grid_->setMargin(1);

		setMinimumSize(250,250); //Canvas (200) + AxisWidget (30) + ScrollBar (20)
		setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::MinimumExpanding);
		
		setAcceptDrops(true);
	}
	
	void SpectrumWidget::setCanvas_(SpectrumCanvas* canvas, UInt row, UInt col)
	{
		canvas_ = canvas;
		setFocusProxy(canvas_);
		grid_->addWidget(canvas_, row, col);
		//axes
		y_axis_ = new AxisWidget(AxisWidget::LEFT, "",this);
		x_axis_ = new AxisWidget(AxisWidget::BOTTOM, "",this);
		grid_->addWidget(y_axis_,row,col-1);
		grid_->addWidget(x_axis_,row+1,col);
		connect(canvas_, SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes()));
		connect(canvas_, SIGNAL(recalculateAxes()), this, SLOT(updateAxes()));
		connect(canvas_, SIGNAL(changeLegendVisibility()), this, SLOT(changeLegendVisibility()));
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
		connect(canvas_, SIGNAL(sendCursorStatus(double,double)), this, SIGNAL(sendCursorStatus(double,double)));
		
		//swap axes if necessary
		updateAxes();
		
		canvas_->setSpectrumWidget(this);
	}
	
	SpectrumWidget::~SpectrumWidget()
	{
		emit aboutToBeDestroyed(window_id);
	}
	
	Int SpectrumWidget::getActionMode() const 
	{
		return canvas_->getActionMode(); 
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
		Histogram<> dist = createIntensityDistribution_();
		HistogramDialog dw(dist);
		dw.setLegend("intensity");
		dw.setLogMode(true);
		if (dw.exec() == QDialog::Accepted)
		{
			DataFilters filters;
			
			if (dw.getLeftSplitter()>dist.minBound())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getLeftSplitter();
				filter.field = DataFilters::INTENSITY;
				filter.op = DataFilters::GREATER_EQUAL;
				filters.add(filter);
			}
			
			if (dw.getRightSplitter()<dist.maxBound())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getRightSplitter();
				filter.field = DataFilters::INTENSITY;
				filter.op = DataFilters::LESS_EQUAL;
				filters.add(filter);
			}
			
			canvas_->setFilters(filters);
		}
	}

	void SpectrumWidget::showMetaDistribution(const String& name)
	{
		Histogram<> dist = createMetaDistribution_(name);
		HistogramDialog dw(dist);
		dw.setLegend(name);
		
		if (dw.exec() == QDialog::Accepted)
		{
			DataFilters filters;
			
			if (dw.getLeftSplitter()>dist.minBound())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getLeftSplitter();
				filter.field = DataFilters::META_DATA;
				filter.meta_name = name;
				filter.op = DataFilters::GREATER_EQUAL;
				filter.value_is_numerical = true;
				filters.add(filter);
			}
			
			if (dw.getRightSplitter()<dist.maxBound())
			{
				DataFilters::DataFilter filter;
				filter.value = dw.getRightSplitter();
				filter.field = DataFilters::META_DATA;
				filter.meta_name = name;
				filter.op = DataFilters::LESS_EQUAL;
				filter.value_is_numerical = true;
				filters.add(filter);
			}
			
			canvas_->setFilters(filters);
		}
	}
	
	void SpectrumWidget::showLegend(bool show)
	{
		y_axis_->showLegend(show);
		x_axis_->showLegend(show);
		update();
	}
	
	void SpectrumWidget::saveAsImage()
	{
		QString file_name = QFileDialog::getSaveFileName(this, "Save File", "", "Images (*.bmp *.png *.jpg *.gif)");
		bool x_visible = x_scrollbar_->isVisible();
		bool y_visible = y_scrollbar_->isVisible();
		x_scrollbar_->hide();
		y_scrollbar_->hide();
		QPixmap pixmap = QPixmap::grabWidget(this);
		x_scrollbar_->setVisible(x_visible);
		y_scrollbar_->setVisible(y_visible);
		pixmap.save(file_name);
	}
	
	void SpectrumWidget::updateAxes()
	{
		//change axis labels if necessary
		if ((canvas()->isMzToXAxis()==true && x_axis_->getLegend().size()>=2 && x_axis_->getLegend().prefix(2)=="RT")
		||( canvas()->isMzToXAxis()==false && y_axis_->getLegend().size()>=2 && y_axis_->getLegend().prefix(2)=="RT")) 
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

	void SpectrumWidget::changeLegendVisibility()
	{
		showLegend(!isLegendShown());
	}

	void SpectrumWidget::closeEvent(QCloseEvent* e)
	{
		for (UInt l=0; l<canvas()->getLayerCount(); ++l)
		{
			//modified => ask if it should be saved		
			const LayerData& layer = canvas()->getLayer(l);
			if (layer.modified)
			{
				QMessageBox::StandardButton result=QMessageBox::question(this,"Save?",(String("Do you want to save your changes to layer '")+ layer.name +  "'?").toQString(),QMessageBox::Ok|QMessageBox::Discard);
				if (result==QMessageBox::Ok)
				{
					canvas()->activateLayer(l);
					canvas()->saveCurrentLayer(false);
				}
			}
		}
		e->accept();
	}
	
	void SpectrumWidget::dragEnterEvent(QDragEnterEvent* event)
	{
		if (event->mimeData()->hasUrls())
		{
			event->acceptProposedAction();
		}
	}
	
	void SpectrumWidget::dragMoveEvent(QDragMoveEvent* event)
	{
		if (event->mimeData()->hasUrls())
		{
			event->acceptProposedAction();
		}
	}
	
	void SpectrumWidget::dropEvent(QDropEvent* event)
	{
		emit dropReceived(event->mimeData(), event->source(), window_id);
		event->acceptProposedAction();
	}

} //namespace OpenMS

