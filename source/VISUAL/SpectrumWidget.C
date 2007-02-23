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

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <QtGui/QGridLayout>
#include <QtGui/QScrollBar>

using namespace std;

namespace OpenMS
{
	
	SpectrumWidget::SpectrumWidget(QWidget* parent)
		: QWidget(parent),
			canvas_(0),
			spectrum_window_(0)
	{
		grid_ = new QGridLayout(this);
		
		//add axes
		y_axis_ = new AxisWidget(AxisWidget::LEFT, "",this);
		x_axis_ = new AxisWidget(AxisWidget::BOTTOM, "",this);
		grid_->addWidget(y_axis_,0,1);
		grid_->addWidget(x_axis_,1,2);

		//add scrollbars
		x_scrollbar_ = new QScrollBar(Qt::Horizontal, this);
		y_scrollbar_ = new QScrollBar(Qt::Vertical, this);
		grid_->addWidget(y_scrollbar_,0,0);
		grid_->addWidget(x_scrollbar_,2,2);		
		x_scrollbar_->hide();
		y_scrollbar_->hide();
	}
	
	
	void SpectrumWidget::setCanvas_(SpectrumCanvas* canvas)
	{
		canvas_ = canvas;
		grid_->addWidget(canvas_, 0, 2);

		//axes
		connect(canvas_, SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes()));
		connect(canvas_, SIGNAL(recalculateAxes()), this, SLOT(updateAxes()));
		//scrollbars
		connect(canvas_, SIGNAL(updateHScrollbar(float,float,float,float)), this, SLOT(updateHScrollbar(float,float,float,float)));
		connect(canvas_, SIGNAL(updateVScrollbar(float,float,float,float)), this, SLOT(updateVScrollbar(float,float,float,float)));
		connect(x_scrollbar_, SIGNAL(valueChanged(int)), canvas_, SLOT(horizontalScrollBarChange(int)));
		connect(y_scrollbar_, SIGNAL(valueChanged(int)), canvas_, SLOT(verticalScrollBarChange(int)));
		
		canvas_->setSpectrumWidget(this);
	}
	
	SpectrumWidget::~SpectrumWidget()
	{
		
	}
	
	SignedInt SpectrumWidget::getActionMode() const 
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
	
	void SpectrumWidget::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
		canvas()->setMainPreferences(prefs);
		
		x_axis_->showLegend(getPrefAsString("Preferences:Legend")=="Show");
		y_axis_->showLegend(getPrefAsString("Preferences:Legend")=="Show");
		
		if (! canvas()->isMzToXAxis())
		{
			// swap legend text
			std::string tmp = x_axis_->getLegend();
			// and update axis with swapped legend text
			x_axis_->setLegend(y_axis_->getLegend());
			y_axis_->setLegend(tmp);
		}
	}
	
	void SpectrumWidget::showIntensityDistribution()
	{
		HistogramDialog dw(createIntensityDistribution_());
		//cout << "showIntensityDistribution: " << canvas_->getCurrentLayer().min_int << " " << canvas_->getCurrentLayer().max_int << endl;
		dw.setLeftSplitter(canvas_->getCurrentLayer().min_int);
		dw.setRightSplitter(canvas_->getCurrentLayer().max_int);
		
		if (dw.exec() == QDialog::Accepted)
		{
			canvas_->setDispInt(dw.getLeftSplitter(), dw.getRightSplitter());
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
		recalculateAxes_();
	}
	
	void SpectrumWidget::mzToXAxis(bool mz_to_x_axis) 
	{
		// check if we have to swap axis
		if (mz_to_x_axis == canvas()->isMzToXAxis()) return;
	
		// update mapping info of cnavas
		canvas()->mzToXAxis(mz_to_x_axis);
	
		// swap legend text
		std::string tmp = x_axis_->getLegend();
		x_axis_->setLegend(y_axis_->getLegend());
		y_axis_->setLegend(tmp);
	
		recalculateAxes_();
	}

	void SpectrumWidget::intensityModeChange_()
	{
		
	}

	QImage SpectrumWidget::getImage(UnsignedInt width, UnsignedInt height)
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
			x_scrollbar_->show();
			x_scrollbar_->setMinimum(static_cast<int>(min));
			x_scrollbar_->setMaximum(static_cast<int>(max-disp_max+disp_min));
			x_scrollbar_->setValue(static_cast<int>(disp_min));
			x_scrollbar_->setPageStep(static_cast<int>(disp_max-disp_min));
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
			//cout << min << " " << disp_min << " " << disp_max << " " << max << endl;
			//cout << min << " " << max-disp_max+disp_min << " " << max-disp_max+min << endl << endl;
			y_scrollbar_->show();
			y_scrollbar_->setMinimum(static_cast<int>(min));
			y_scrollbar_->setMaximum(static_cast<int>(max-disp_max+disp_min));
			y_scrollbar_->setValue(static_cast<int>(max-disp_max+min));
			y_scrollbar_->setPageStep(static_cast<int>(disp_max-disp_min));
		}
	}

} //namespace OpenMS

