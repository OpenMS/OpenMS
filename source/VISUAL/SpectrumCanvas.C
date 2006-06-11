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
// $Id: SpectrumCanvas.C,v 1.28 2006/06/08 14:29:19 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/SpectrumWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>

// custom cursor for translating
#include "ICONS/handopen.xpm"
#include "ICONS/handclosed.xpm"

// QT
#include <qtimer.h>
#include <qpainter.h>
#include <qpixmap.h>
#include <qbitmap.h>

// STLasdsd
#include <iostream>


using namespace std;
namespace OpenMS
{

	SpectrumCanvas::SpectrumCanvas(QWidget* parent, const char* name, WFlags f)
		: QScrollView(parent, name, f | Qt::WNoAutoErase | Qt::WStaticContents),
			buffer_(0),
			action_mode_(AM_ZOOM),
			intensity_modification_(IM_NONE),
			min_disp_ints_(),
			max_disp_ints_(),
			mapping_info_(),
			pen_width_(0),
			show_grid_(true),
			recalculate_(false),
			zoom_timeout_(false),
			current_data_(0),
			spectrum_widget_(0),
			datasets_()
	{
		setFrameShape(NoFrame);
		connect(this, SIGNAL(contentsMoving(int, int)), this, SLOT(move_(int, int)));
	
	  createCustomMouseCursors_();
	  
		// we need to initialize the painting buffer in the
		// constructor for maximum performance
		adjustBuffer_();
	}
	
	void SpectrumCanvas::createCustomMouseCursors_()
	{
		// create custom mouse cursor for translate action as Qt doesn't provide one
		QPixmap* pm1 = new QPixmap(handopen);
		pm1->setMask(pm1->createHeuristicMask());
		cursor_translate_ = QCursor(*pm1);
	
		QPixmap* pm2 = new QPixmap(handclosed);
		pm2->setMask(pm2->createHeuristicMask());
		cursor_translate_in_progress_ = QCursor(*pm2); 
	}
	
	void SpectrumCanvas::viewportResizeEvent(QResizeEvent* /* e */)
	{
		adjustBuffer_();
		updateScrollbars_();
		invalidate_();
	}
	
	void SpectrumCanvas::viewportPaintEvent(QPaintEvent* e)
	{
		// blit changed parts of backbuffer on widget
		QMemArray<QRect> rects = e->region().rects();
		for (int i = 0; i < (int)rects.size(); ++i)
		{
			bitBlt(viewport(), rects[i].topLeft(), buffer_, rects[i]);
		}
	}
	
	void SpectrumCanvas::adjustBuffer_()
	{
		if (buffer_) delete buffer_;
		
		buffer_ = new QPixmap(width(), height());
	}
	
	void SpectrumCanvas::setDispInt(double min, double max)
	{
		min_disp_ints_[current_data_] = min;
		max_disp_ints_[current_data_] = max;
		intensityDistributionChange_();
	}
	
	void SpectrumCanvas::showGridLines(bool show)
	{
		show_grid_ = show;
		invalidate_();
	}
	
	void SpectrumCanvas::intensityDistributionChange_()
	{
		recalculate_ = true;
		invalidate_();
	}
	
	void SpectrumCanvas::intensityModificationChange_()
	{
		recalculate_ = true;
		invalidate_();
	}
	
	void SpectrumCanvas::setMirroredXAxis(bool b)
	{
		if (b)
		{
			mapping_info_.setXAxisDesc();
		}
		else
		{
			mapping_info_.setXAxisAsc();
		}
		invalidate_();
	}
	
	void SpectrumCanvas::setMirroredYAxis(bool b)
	{
		if (b)
		{
			mapping_info_.setYAxisDesc();
		}
		else
		{
			mapping_info_.setYAxisAsc();
		}
		invalidate_();
	}
	
	SpectrumCanvas::PointType SpectrumCanvas::contextToChart_(QPoint p, int width, int height)
	{
		// transform point and context if necessary
		if (!mapping_info_.isXAxisAsc()) p.setX(width - p.x());
		if (mapping_info_.isYAxisAsc()) p.setY(height - p.y());
		
		if (!mapping_info_.isMzToXAxis())
		{
			int tmp = p.x();
			p.setX(p.y());
			p.setY(tmp);
			
			std::swap(width, height);
		}
		
		// calculate point4
		const float pixel_width = visible_area_.width() / static_cast<float>(width);
		const float pixel_height = visible_area_.height() / static_cast<float>(height);
		
		return PointType(visible_area_.minX() + pixel_width * static_cast<float>(p.x()),
		                    visible_area_.minY() + pixel_height * static_cast<float>(p.y()));
	}
	
	QPoint SpectrumCanvas::chartToContext_(const PointType& pos, int width, int height)
	{
//		cout << "Pos          -- X    : " << pos.X() << " Y: " << pos.Y() << endl;
//		cout << "Pixel        -- Width: " << width << " Height: " << height << endl;
//		cout << "Visible area -- Width: " << visible_area_.width() << " Height: " << visible_area_.height() << endl;
//		cout << "             -- X    : " << visible_area_.minX() << " " << visible_area_.maxX() << endl;
//		cout << "             -- Y    : " << visible_area_.minY() << " " << visible_area_.maxY() << endl;
		if (!mapping_info_.isMzToXAxis())
		{
			//cout << "SWAP" << endl;
			std::swap(width, height);
		}
		
		// calculate point for common case
		const float pixel_width = visible_area_.width() / static_cast<float>(width);
	
		const float pixel_height = visible_area_.height() / static_cast<float>(height);
		
		QPoint p(static_cast<int>((pos.X() - visible_area_.minX()) / pixel_width),
		         static_cast<int>((pos.Y() - visible_area_.minY()) / pixel_height));
		
		//cout << "Calculated   -- Width: " << p.x() << " Height: " << p.y() << endl;
		
		// transform point if necessary
		if (!mapping_info_.isMzToXAxis())
		{
			std::swap(width, height);
			
			int tmp = p.x();
			p.setX(p.y());
			p.setY(tmp);
		}
		
		if (mapping_info_.isYAxisAsc()) p.setY(height - p.y());
		if (!mapping_info_.isXAxisAsc()) p.setX(width - p.x());
		
		return p;
	}
	
	SpectrumCanvas::PointType SpectrumCanvas::widgetToChart_(const QPoint& pos)
	{
		return contextToChart_(pos, viewport()->width(), viewport()->height());
	}
	
	QPoint SpectrumCanvas::chartToWidget_(const PointType& pos)
	{
		return chartToContext_(pos, viewport()->width(), viewport()->height());
	}
	
	void SpectrumCanvas::changeVisibleArea_(const AreaType& new_area)
	{
		//store old zoom state
		if (zoom_timeout_)
		{
			zoom_stack_.push(visible_area_);
		}
		
		visible_area_ = new_area;
		
		updateScrollbars_();
		
		zoom_timeout_ = false;
		QTimer::singleShot(2000, this, SLOT(timeoutZoom_()));
	
		emit visibleAreaChanged(new_area);
		recalculate_ = true;
		invalidate_();
	}
	
	
	void SpectrumCanvas::updateScrollbars_()
	{
		const AreaType& area_ = getDataArea_();
		QPoint left_top = - chartToWidget_(PointType(area_.minX(), area_.minY()));
		QPoint size = left_top + chartToWidget_(PointType(area_.maxX(), area_.maxY()));
		
		// block contentsMoving signal, because calling the private
		// move_() slot falsely changes the visible area
		blockSignals(true);
		resizeContents(size.x(), size.y());
		setContentsPos(left_top.x(), left_top.y());
		blockSignals(false);
	}
	
	void SpectrumCanvas::move_(int x, int y)
	{
		PointType left_top = widgetToChart_(contentsToViewport(QPoint(x, y)));
		PointType right_bottom = left_top + (widgetToChart_(QPoint(viewport()->width(), viewport()->height()))
		                                      - widgetToChart_(QPoint(0, 0)));
		
	  visible_area_ = AreaType(left_top, right_bottom);
		emit visibleAreaChanged(visible_area_);
		recalculate_ = true;
		invalidate_();
	}
	
	void SpectrumCanvas::timeoutZoom_()
	{
		zoom_timeout_ = true;
	}
	
	
	void SpectrumCanvas::zoomBack_()
	{
		if (zoom_stack_.empty())
		{
			visible_area_ = getDataArea_();
		}
		else
		{
			visible_area_ = zoom_stack_.top();
			zoom_stack_.pop();
		}
		
		updateScrollbars_();
		emit visibleAreaChanged(visible_area_);
		recalculate_ = true;
		invalidate_();
	}
	
	void SpectrumCanvas::resetZoom()
	{
		changeVisibleArea_(getDataArea_());
	}
	
	void SpectrumCanvas::setVisibleArea(AreaType area)
	{
		if (visible_area_ != area)
		{
			changeVisibleArea_(area);
		}
	}
	
	
	void SpectrumCanvas::paintGridLines_(QPainter* p)
	{
		if (!show_grid_ || !spectrum_widget_) return;
	
		p->save();
	
		unsigned int xl, xh, yl, yh; //width/height of the diagram area, x, y coordinates of lo/hi x,y values
	
		if (mapping_info_.isXAxisAsc())
		{
			xl = 0;
			xh = viewport()->width();
		}
		else
		{
			xl = viewport()->width();
			xh = 0;
		}
	
		if (mapping_info_.isYAxisAsc())
		{
			yl = viewport()->height();
			yh = 0;
		}
		else
		{
			yl = 0;
			yh = viewport()->height();
		}
	
		// drawing of grid lines and associated text	
		QColor grid_line_color;
	
		for (unsigned int j = 0; j != spectrum_widget_->xAxis()->gridLines().size() ; j++) 
		{
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
				grid_line_color = QColor(90, 90, 90);
				break;
				case 1:	// style settings for small intervals
					grid_line_color = QColor(130, 130, 130);
				break;
				case 2: // style settings for smalles intervals
					grid_line_color = QColor(170, 170, 170);
				break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					grid_line_color = QColor(0, 0, 0);
				break;
			}
			p->setPen(QPen(grid_line_color,pen_width_));
			int x;
			for (std::vector<double>::const_iterator it = spectrum_widget_->xAxis()->gridLines()[j].begin(); it != spectrum_widget_->xAxis()->gridLines()[j].end(); it++) 
			{
				x = static_cast<int>(intervalTransformation(*it, spectrum_widget_->xAxis()->getAxisMinimum(), spectrum_widget_->xAxis()->getAxisMaximum(), xl, xh));
				p->drawLine(x, yl, x, yh);
			}
		}
		
		for (unsigned int j = 0; j != spectrum_widget_->yAxis()->gridLines().size() ; j++) 
		{
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					grid_line_color = QColor(0, 0, 0);
					break;
				case 1:	// style settings for small intervals
					grid_line_color = QColor(90, 90, 90);
					break;
				case 2: // style settings for smalles intervals
					grid_line_color = QColor(140, 140, 140);
					break;
				default:
					std::cout << "empty vertical grid line vector error!" << std::endl;
					grid_line_color = QColor(0, 0, 0);
				break;
			}
			p->setPen(QPen(grid_line_color,pen_width_));
			int y;
			for (std::vector<double>::const_iterator it = spectrum_widget_->yAxis()->gridLines()[j].begin(); it != spectrum_widget_->yAxis()->gridLines()[j].end(); it++) 
			{
				y = static_cast<int>(intervalTransformation(*it, spectrum_widget_->yAxis()->getAxisMinimum(), spectrum_widget_->yAxis()->getAxisMaximum(), yl, yh));
				p->drawLine(xl, y, xh, y);
			}
		}
		
		p->restore();
	}
	
	void SpectrumCanvas::setMainPreferences(const Param& prefs)
	{
		prefs_ = prefs;
	}
	
	UnsignedInt SpectrumCanvas::activeDataSetIndex() const
	{
		return current_data_;	
	}
	
	UnsignedInt SpectrumCanvas::getDataSetCount() const
	{
		return datasets_.size();
	}
	
	const String& SpectrumCanvas::getDataSetName(UnsignedInt index) const
	{
		if (index >= getDataSetCount())
		{
			return String::EMPTY;
		}
		return datasets_[index].getName();
	}
	
	bool SpectrumCanvas::isDataSetVisible(UnsignedInt index) const
	{
		return layer_visible_[index];
	}		

	SpectrumCanvas::ExperimentType& SpectrumCanvas::addEmptyDataSet()
	{
		datasets_.resize(getDataSetCount()+1);
		return datasets_[getDataSetCount()-1];
	}

	SpectrumCanvas::ExperimentType& SpectrumCanvas::getDataSet(UnsignedInt index) throw (Exception::IndexOverflow)
	{
		if (index >= getDataSetCount())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,index,getDataSetCount()-1);
		}
		return datasets_[index];
	}

	SpectrumCanvas::ExperimentType& SpectrumCanvas::currentDataSet() throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		return datasets_[current_data_];
	}

	const SpectrumCanvas::ExperimentType& SpectrumCanvas::getDataSet(UnsignedInt index) const throw (Exception::IndexOverflow)
	{
		if (index >= getDataSetCount())
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,index,getDataSetCount()-1);
		}
		return datasets_[index];
	}

	const SpectrumCanvas::ExperimentType& SpectrumCanvas::currentDataSet() const throw (Exception::IndexOverflow)
	{
		if (getDataSetCount()==0)
		{
			throw Exception::IndexOverflow(__FILE__, __LINE__, __PRETTY_FUNCTION__,0,0);
		}
		return datasets_[current_data_];
	}

	SignedInt SpectrumCanvas::addDataSet(const ExperimentType& in)
	{	
		datasets_.push_back(in);
		return finishAdding();
	}


	void SpectrumCanvas::changeVisibility(int i , bool b)
	{
		if (layer_visible_[i]!=b)
		{
			layer_visible_[i]=b;
			invalidate_();
		}
	}	

} //namespace


