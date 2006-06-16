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
// $Id: SpectrumWidget.C,v 1.25 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>
#include <OpenMS/VISUAL/AxisWidget.h>

#include <qaction.h>
#include <qlayout.h>
#include <qimage.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	
	SpectrumWidget::SpectrumWidget(QWidget* parent, const char* name, WFlags f)
		: QWidget(parent,name,f),
			canvas_(0),
			old_max_intensity_(0.0),
			spectrum_window_(0)
	{
		grid_ = new QGridLayout(this, 2, 2);
		
		QVBoxLayout* vbox = new QVBoxLayout();
		QHBoxLayout* hbox = new QHBoxLayout();
		
		y_axis_ = new AxisWidget(AxisWidget::LEFT, "",this);
		x_axis_ = new AxisWidget(AxisWidget::TOP, "",this);
	
		vspacer_ = new QWidget(this);
		hspacer_ = new QWidget(this);
	
		addClient(y_axis_,"Y-Axis",true);
		addClient(x_axis_,"X-Axis",true);
	
		y_axis_->setPaletteBackgroundColor(backgroundColor());
		x_axis_->setPaletteBackgroundColor(backgroundColor());
		y_axis_->showLegend(show_legend_);
		x_axis_->showLegend(show_legend_);
	
		vspacer_->hide();
		hspacer_->hide();
		
//		vspacer_->setBackgroundColor(Qt::red); //For debugging
//		hspacer_->setBackgroundColor(Qt::red); //For debugging
//		x_axis_->setBackgroundColor(Qt::green); //For debugging
//		y_axis_->setBackgroundColor(Qt::green); //For debugging
		
		vspacer_->setSizePolicy(QSizePolicy(QSizePolicy::Minimum, QSizePolicy::Fixed));
		hspacer_->setSizePolicy(QSizePolicy(QSizePolicy::Fixed, QSizePolicy::Minimum));
		
		vbox->addWidget(y_axis_);
		vbox->addWidget(vspacer_);
		hbox->addWidget(x_axis_);
		hbox->addWidget(hspacer_);
		
		grid_->addLayout(vbox, 1, 0);
		grid_->addLayout(hbox, 0, 1);
	
		grid_->setRowStretch(0,1);
		grid_->setRowStretch(1,20);
		grid_->setColStretch(0,1);
		grid_->setColStretch(1,20);
	}
	
	
	void SpectrumWidget::setCanvas(SpectrumCanvas* canvas)
	{
		canvas_ = canvas;
		grid_->addWidget(canvas_, 1, 1);
	
		connect(canvas_, SIGNAL(contextMenu(QPoint)), this, SIGNAL(contextMenu(QPoint)));
		connect(canvas_, SIGNAL(visibleAreaChanged(DRange<2>)), this, SLOT(updateAxes_(DRange<2>)));
		
		canvas_->setSpectrumWidget(this);
	}
	
	SpectrumWidget::~SpectrumWidget()
	{
	}
	
	void SpectrumWidget::showGridLines(bool show)
	{
		canvas_->showGridLines(show);
	}
	
	SignedInt SpectrumWidget::getActionMode() const 
	{
		return canvas_->getActionMode(); 
	}
	
	void SpectrumWidget::setActionMode(QAction* a)
	{
		QString name;
		name=a->name();
		if(name=="setPickAction") actionModeSelect(); 
		else if(name=="setZoomAction") actionModeZoom(); 
		else if(name=="setTranslateAction") actionModeTranslate();
		else if(name=="setMeasureAction") actionModeMeasure();
		//todo else throw new std::runtime_error("unknown QAction");
	}
	
	// ### is this function really necessary?
	void SpectrumWidget::setActionMode(SignedInt mode)
	{
		if(mode==SpectrumCanvas::AM_SELECT) actionModeSelect(); 
		else if(mode==SpectrumCanvas::AM_ZOOM) actionModeZoom(); 
		else if(mode==SpectrumCanvas::AM_TRANSLATE) actionModeTranslate();
	}
	
	void SpectrumWidget::actionModeSelect()
	{
		if (getActionMode() != SpectrumCanvas::AM_SELECT)
		{
			canvas_->setActionMode(SpectrumCanvas::AM_SELECT);
			emit modesChanged(this);
		}
	}
	
	void SpectrumWidget::actionModeZoom()
	{
		if (getActionMode() != SpectrumCanvas::AM_ZOOM)
		{
			canvas_->setActionMode(SpectrumCanvas::AM_ZOOM);
			emit modesChanged(this);
		}
	}
	
	void SpectrumWidget::actionModeTranslate()
	{
		if (getActionMode() != SpectrumCanvas::AM_TRANSLATE)
		{
			canvas_->setActionMode(SpectrumCanvas::AM_TRANSLATE);
			emit modesChanged(this);
		}
	}
	
	void SpectrumWidget::actionModeMeasure()
	{
		if (getActionMode() != SpectrumCanvas::AM_MEASURE)
		{
			canvas_->setActionMode(SpectrumCanvas::AM_MEASURE);
			emit modesChanged(this);
		}
	}
	
	void SpectrumWidget::setIntensityModificationNone()
	{
		if (canvas_->getIntensityModification() != SpectrumCanvas::IM_NONE) {
			canvas_->setIntensityModification(SpectrumCanvas::IM_NONE);
			intensityModificationChange_();
		}
	}
	
	void SpectrumWidget::setIntensityModificationLog()
	{
		if (canvas_->getIntensityModification() != SpectrumCanvas::IM_LOG) {
			canvas_->setIntensityModification(SpectrumCanvas::IM_LOG);
			intensityModificationChange_();
		}
	}
	
	void SpectrumWidget::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
		canvas()->setMainPreferences(prefs);
	
		setMirroredXAxis(! canvas()->getMappingInfo().isXAxisAsc());
		setMirroredYAxis(! canvas()->getMappingInfo().isYAxisAsc());
	
		if (! canvas()->getMappingInfo().isMzToXAxis())
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
		//cout << "Min: " << canvas_->getMinDispInt() << "  Max: " << canvas_->getMaxDispInt() << endl;
		dw.setLeftSplitter(canvas_->getMinDispInt());
		dw.setRightSplitter(canvas_->getMaxDispInt());
		
		if (dw.exec() == QDialog::Accepted)
		{
			canvas_->setDispInt(dw.getLeftSplitter(), dw.getRightSplitter());
			emit sendStatusMessage("Displayed intensity range: "+String(dw.getLeftSplitter())+" upto "+String(dw.getRightSplitter())+" m/z", 5000);
		}
	}
	
	void SpectrumWidget::showLegend()
	{
		show_legend_ = true;
		legendModificationChange_();
	}
	
	void SpectrumWidget::showNoLegend()
	{
		show_legend_ = false;
		legendModificationChange_();
	}
	
	bool SpectrumWidget::isLogIntensity() const
	{
		return canvas_->getIntensityModification() == SpectrumCanvas::IM_LOG;
	}
	
	void SpectrumWidget::updateAxes_(SpectrumCanvas::AreaType /*area*/)
	{
		recalculateAxes();
		
		//show hide spaces when scollbars are shown/hidden
		if (canvas_->visibleWidth() < canvas_->contentsWidth())
		{
			vspacer_->setMinimumSize(y_axis_->width(), canvas_->horizontalScrollBar()->sizeHint().height());
			vspacer_->show();
		}
		else
		{
			vspacer_->hide();
		}
		
		if (canvas_->visibleHeight() < canvas_->contentsHeight())
		{
			hspacer_->setMinimumSize(canvas_->verticalScrollBar()->sizeHint().width(), x_axis_->height());
			hspacer_->show();
		}
		else
		{
			hspacer_->hide();
		}
	}
	
	void SpectrumWidget::setMirroredXAxis(bool b)
	{
		canvas()->setMirroredXAxis(b);
		x_axis_->setInverseOrientation(b);
		x_axis_->update();
	}
	
	void SpectrumWidget::setMirroredYAxis(bool  b )
	{
		canvas()->setMirroredYAxis(b);
		y_axis_->setInverseOrientation(b);
		y_axis_->update();
	}
	
	void SpectrumWidget::switchAxis(bool swapped_axes) 
	{
		// check if we have to swap axis
		if (swapped_axes != canvas()->getMappingInfo().isMzToXAxis()) return;
	
		// get current axis orientation
		bool mirrored_x_axis = x_axis_->hasInverseOrientation();
		bool mirrored_ordinate = y_axis_->hasInverseOrientation();
	
		// update mapping info
		canvas()->getMappingInfo().isMzToXAxis() ? canvas()->getMappingInfo().setMzToYAxis() : canvas()->getMappingInfo().setMzToXAxis();
	
		// swap information about inverse orientation
		x_axis_->setInverseOrientation(mirrored_ordinate);
		y_axis_->setInverseOrientation(mirrored_x_axis);
	
		// swap legend text
		std::string tmp = x_axis_->getLegend();
		// and update axis with swapped legend text
		x_axis_->setLegend(y_axis_->getLegend());
		y_axis_->setLegend(tmp);
	
		recalculateAxes();
		canvas()->update();
	}

	QImage SpectrumWidget::getImage(UnsignedInt width, UnsignedInt height)
	{
		//hide();
		this->setUpdatesEnabled(false);
		
		//hide scrollbars and spacers if necessary
		bool hor_sc_off_ = (canvas()->hScrollBarMode()==QScrollView::AlwaysOff);
		bool ver_sc_off_ = (canvas()->vScrollBarMode()==QScrollView::AlwaysOff);
		if (!hor_sc_off_)
		{
			canvas()->setHScrollBarMode(QScrollView::AlwaysOff);
			hspacer_->hide();
		}
		if (!ver_sc_off_)
		{
			canvas()->setVScrollBarMode(QScrollView::AlwaysOff);
			vspacer_->hide();
		}
		
		//store old background colors and size
	 	QColor c_a = x_axis_->paletteBackgroundColor();
	 	QColor c_o = y_axis_->paletteBackgroundColor();
	 	QColor c_t = paletteBackgroundColor();
		
		//set white background
		y_axis_->setPaletteBackgroundColor(white);
		x_axis_->setPaletteBackgroundColor(white);
		setBackgroundColor(white);
	
		//set pen width
		int pen = width/1024;
		canvas()->setPenWidth(pen);
		y_axis_->setPenWidth(pen);
		x_axis_->setPenWidth(pen);
		
	 	//store size and resize
		int h = this->height();
		int w = this->width();

		grid_->activate();
		resize(width,height);		
				
		//take an image
		QPixmap p = QPixmap::grabWidget(this);
		QImage image = p.convertToImage();
		
		//resore background colors
	 	y_axis_->setPaletteBackgroundColor(c_o);
	 	x_axis_->setPaletteBackgroundColor(c_a);
		setBackgroundColor(c_t);
		
		//restore pen withs
		canvas()->setPenWidth(0);
		y_axis_->setPenWidth(0);
		x_axis_->setPenWidth(0);
		
		//show scrollbars and spacers again
		if (!hor_sc_off_)
		{
			canvas()->setHScrollBarMode(QScrollView::Auto);
			hspacer_->show();
		}
		if (!ver_sc_off_)
		{
			canvas()->setVScrollBarMode(QScrollView::Auto);
			vspacer_->show();
		}
		
		//restore size
		resize(w,h);
		
		this->setUpdatesEnabled(true);
		//show();
		
		return image;
	}

} //namespace OpenMS

