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

// Qt
#include <qtooltip.h>
#include <qpixmap.h>

// STL
#include <iostream>


// OpenMS
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>

// ANSI C/C++
#include "stdio.h"
#include "stdlib.h"
#include "math.h"

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	AxisWidget::AxisWidget(UnsignedInt alignment, const char* legend, QWidget* parent, const char* name, WFlags f)
		: QWidget( parent, name, f | WRepaintNoErase),
		is_log_(false),
		show_legend_(false),
		alignment_(alignment),
		inverse_orientation_(false),
		margin_(0),
		legend_(legend),
		tick_level_(3),
		buffer_(0),
		pen_width_(0)
	{
		buffer_ = new QPixmap(1,1);
		if (!(alignment==RIGHT || alignment==LEFT || alignment==BOTTOM || alignment==TOP))
		{
			throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		setAxisBounds(0.0,100.0);  // min_, max_
	 	if (alignment==RIGHT || alignment==LEFT )	
	 	{
			setMinimumSize(30,100);	 
		}
		else
		{
			setMinimumSize(100, 30);
		}
		resize(minimumSize());
		setMouseTracking(true);
	}
	
	
	AxisWidget::~AxisWidget()
	{
		delete(buffer_);
	}
	
	void AxisWidget::resizeEvent(QResizeEvent* e)
	
	{
		buffer_->resize(e->size().width(),e->size().height());
		invalidate_();
	}
	
	void AxisWidget::resize(QSize s)
	{
		resize(s.width(), s.height());
	}
	
	void AxisWidget::resize(int w, int h)
	{
		bool isXAxis = (alignment_==BOTTOM || alignment_==TOP);
	
		int h_plus_m = (isXAxis)? h: h+margin_;
		int w_plus_m = (isXAxis)? w+margin_ : w;
	
		this->QWidget::resize(w_plus_m,h_plus_m);
		buffer_->resize(w_plus_m, h_plus_m);
		invalidate_();
	}
	
	void AxisWidget::paintEvent(QPaintEvent *)
	{
		bitBlt(this, 0, 0, buffer_);
	}
	
	void AxisWidget::invalidate_()
	{
		QColor text_color;
		int tick_size = 0;
	
		buffer_->fill(paletteBackgroundColor());
	 	painter_.begin(buffer_);
	
		bool isXAxis = (alignment_==BOTTOM || alignment_==TOP);
		int h = buffer_->height();
		int w = buffer_->width();
		w = (isXAxis)? w-margin_ : w;  // Remove margin to get the scale right
		h = (isXAxis)? h : h-margin_;
	
		// probe legend font size (use a third (fifth) of available space for legend)
		int legend_font_size = (isXAxis)? probeFont_(QString(legend_.c_str()),w,h/3.0) :  probeFont_(QString(legend_.c_str()),h,w/5.0);
	
		// probe text font size with last axis labels  
		QString probe = "";
		int index = 0;
		for (unsigned int i=0; i<3; i++)
		{
			QString s;
			if (grid_line_.size()>i)
			{
				s = QString("%1").arg(scale_(grid_line_[i][ grid_line_[i].size()-1])); //last axis label
			}
			if (s.length() > probe.length())
			{
				probe = s;
				index = i;
			}
		}
	
		double trans_dist = intervalTransformation(grid_line_dist_+min_, min_, max_, 0, (isXAxis)? w:h );
	
		int text_font_size;
		// use rest of space for axis label (i.e. 1/5 legend, 4/5 label) or use all available space if no legend is shown;
		if (show_legend_)
		{
			text_font_size = (isXAxis)? probeFont_(probe,trans_dist,0.8*h*2.0/3) : probeFont_(probe,0.8*w,trans_dist,index);
		}
		else
		{
			text_font_size = (isXAxis)? probeFont_(probe,trans_dist,0.8*h) : probeFont_(probe,w,trans_dist,index);
		}
		
		// length of one big tick equals a quarter of the text_font_size
		double grid_scaling = text_font_size/4;
		if (grid_scaling > w/5 ) grid_scaling = w/5;  // don't let font size get to big
		if (grid_scaling > h/5 ) grid_scaling = h/5;  
	
		for (unsigned int j = 0; j != grid_line_.size() ; j++) 
		{
	    if (is_log_ && j>0) break; // just draw text on big intervalls
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter_.setFont(QFont("courier", static_cast<UnsignedInt>(3*grid_scaling))); 
					tick_size = static_cast<UnsignedInt>( grid_scaling );
					text_color = QColor(0, 0, 0);				
				break;
				case 1:	// style settings for small intervals
					painter_.setFont(QFont("courier",static_cast<UnsignedInt>(2.4*grid_scaling))); 
					tick_size = static_cast<int>( grid_scaling*2/3 );
					text_color = QColor(20, 20, 20);
				break;
				case 2: // style settings for smalles intervals
					painter_.setFont(QFont("courier",static_cast<UnsignedInt>(2.0*grid_scaling))); 
					tick_size = static_cast<UnsignedInt>( grid_scaling/2);
					text_color = QColor(40, 40, 40);
				break;
				default:
					std::cerr << "empty grid line vector error! in Line: " << __LINE__ << " File: " << __FILE__ << std::endl;
					painter_.setFont(QFont("courier", static_cast<UnsignedInt>(3*grid_scaling))); 
					tick_size = static_cast<UnsignedInt>( grid_scaling);
					text_color = QColor(0, 0, 0);				
				break;
			}
	
			// drawing
			int x;
			UnsignedInt i_beg = (isXAxis)? 0 : h;
			UnsignedInt i_end = (isXAxis)? w : 0;
	
			for (std::vector<double>::const_iterator it = grid_line_[j].begin(); it != grid_line_[j].end(); it++) 
			{
	      if (inverse_orientation_)
	      {
					x = static_cast<int>(intervalTransformation(*it, min_, max_, i_end, i_beg)) + ((alignment_==LEFT || alignment_==RIGHT)?-1:1)*margin_;
				}
				else
				{
					x = static_cast<int>(intervalTransformation(*it, min_, max_, i_beg, i_end));
				}
	
				// small axis lines
				painter_.setPen(QPen(QPainter::black, pen_width_));
				switch (alignment_)	
				{
					case BOTTOM: painter_.drawLine(x, 0, x, tick_size); break;
					case TOP: painter_.drawLine(x, h, x,  h-tick_size); break;
					case LEFT: painter_.drawLine(w-tick_size, x+margin_, w, x+margin_); break;
					case RIGHT: painter_.drawLine(0, x+margin_, tick_size, x+margin_);	break;
				}			
					
				// values at axis lines
				QString text = QString("%1").arg(scale_(*it));
				painter_.setPen(QPen(text_color, pen_width_));
		
				// get bounding rectangle for text we want to layout
				QRect textbound = QFontMetrics(painter_.font()).boundingRect(text);
	
				// Calculate text position
				UnsignedInt posx = 0, posy = 0;
				switch (alignment_)	
				{
				  case BOTTOM: case TOP:	posx = x - static_cast<UnsignedInt>( textbound.width()/2 ); break;  // Center text around tick
			  	case LEFT: posx = w - static_cast<UnsignedInt>(1.5*grid_scaling) - textbound.width(); break;
				  case RIGHT: posx = static_cast<UnsignedInt>(1.5*grid_scaling); break; // Leave space for tick (longest one = grid_scaling)
				}	
	
				switch (alignment_)	
				{
				  case BOTTOM: posy = static_cast<UnsignedInt>(1.5*grid_scaling) + textbound.height() ; break;
				  case TOP:	posy = h-static_cast<UnsignedInt>(1.5*grid_scaling); break;
				  case LEFT: case RIGHT: posy = x+margin_ + static_cast<UnsignedInt>(textbound.height()/ 2);  break;
				}	
				painter_.drawText(posx, posy, text);
			}		
		}
	
		// style settings for legend
		painter_.setBrush(QPainter::NoBrush);
		painter_.setPen(QPen(QPainter::black,pen_width_));
		
		// text label at axis
		if (show_legend_ && legend_!="") 
		{	
			// style settings
			painter_.setFont(QFont("courier", legend_font_size));//static_cast<UnsignedInt>(2.4*grid_scaling)));
			painter_.setPen(QPen(QPainter::black,pen_width_));
	
			switch (alignment_)	
			{
	    		case BOTTOM:
						painter_.drawText(0, 0 ,  w, h, Qt::AlignBottom|Qt::AlignHCenter, legend_.c_str());	
						break;
				  case TOP:
						painter_.drawText(0, 0 ,  w, h, Qt::AlignTop|Qt::AlignHCenter, legend_.c_str());	
						break;
				  case LEFT: 
						painter_.rotate(270);
	 				painter_.drawText(-h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignTop, legend_.c_str());
						break;
				  case RIGHT:
						painter_.rotate(270);
						painter_.drawText(-h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignBottom, legend_.c_str());
						break;
			}	
		}
	
		painter_.end();
		update();
	}
	
	void AxisWidget::setAxisBounds(double min, double max)
	{
	  if (min == max)
	  {
	    return;
	  }
	  
		bool change=true;
		if (is_log_)
		{
			if (min_== linear2log(min) && max_ == linear2log(max)) change=false;
			min_ = linear2log(min);
			max_ = linear2log(max);
		}
		else
		{
			if (min_==min && max_==max) change=false;
			min_ = min; 
			max_ = max;
		}
	
		if (change)
	  {
	    if (is_log_)
	    {
				AxisTickCalculator::calcLogGridLines(min_,max_,grid_line_);
			} else
	    {	
				UnsignedInt max_num_big = (alignment_==BOTTOM || alignment_==TOP)? 7: 5;
				UnsignedInt max_num_small = (alignment_==BOTTOM || alignment_==TOP)? 5: 3;
				AxisTickCalculator::calcGridLines(min_, max_,tick_level_,grid_line_,max_num_big, max_num_small,grid_line_dist_);
	    }
			invalidate_();
		}
	
	}
	
	void AxisWidget::setLogScale(bool is_log)
	{
	  if (is_log_ != is_log)
	  {
	    is_log_ = is_log;
			
			if (is_log_) setAxisBounds(min_,max_);  // reset AxisBound - side effect: logarithmic transformation
			else setAxisBounds(log2linear(min_),log2linear(max_));  // reverse log transformation
	
		  invalidate_();
	  }
	}
	
	bool AxisWidget::isLogScale()
	{
		return is_log_;
	}
	
	UnsignedInt AxisWidget::margin()
	{
		return margin_;
	}
	
	void AxisWidget::setMargin(UnsignedInt margin)
	{
	  if (margin_ != margin){
	    margin_ = margin;
		  invalidate_();
	  }
	}
	
	
	void AxisWidget::showLegend(bool show_legend)
	{
		if (show_legend_ != show_legend)
		{
	  	show_legend_ = show_legend;
		  if (show_legend_)
		  { 
			  QToolTip::remove(this);
	  	}
	  	else
	  	{
		  	QToolTip::add(this,legend_.c_str());
	    }
	    invalidate_();
	  }
	}
	
	bool AxisWidget::isLegendShown() const
	{
		return show_legend_;
	}
	
	std::string AxisWidget::getLegend()
	{
		return legend_;
	}
	
	void AxisWidget::setInverseOrientation(bool inverse_orientation)
	{
		if (inverse_orientation_ != inverse_orientation){
	  	inverse_orientation_ = inverse_orientation;
	    invalidate_();
	  }
	}
	
	bool AxisWidget::hasInverseOrientation()
	{
		return inverse_orientation_;
	}
	
	void AxisWidget::setLegend(const std::string& legend)
	{
		if (legend_ != legend)
		{
			legend_ = legend;
			if (!show_legend_)
	  	{ 
			  QToolTip::remove(this);
	 	 		QToolTip::add(this,legend_.c_str());
   	 	}
  	}
	}
	
	QSize AxisWidget::sizeHint () const
	{
		return(QSize(width(),height()));
	}
	
	bool AxisWidget::setTickLevel(UnsignedInt level)
	{
		if (1<=level && level<=3){
			tick_level_ = level;
			return true;
		}else{
			return false;
		}
	}

	const AxisWidget::GridVector& AxisWidget::gridLines()
	{
		return grid_line_;
	}
	
	void AxisWidget::setGridLines(std::vector<double>& grid)
	{
		grid_line_.clear();
		grid_line_.push_back(grid);
	}
	
	void AxisWidget::mouseMoveEvent( QMouseEvent* /*e*/)
	{
		emit sendAxisMouseMovement();
	}

} //Namespace



