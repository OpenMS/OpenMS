// -*- mode: C++; tab-width: 2; -*-
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

// Qt
#include <QtGui/QPaintEvent>
#include <QtGui/QPainter>

// STL
#include <iostream>
#include <algorithm>

// OpenMS
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/VISUAL/AxisTickCalculator.h>

using namespace std;

namespace OpenMS
{
	using namespace Math;
	
	AxisWidget::AxisWidget(UInt alignment, const char* legend, QWidget* parent)
		: QWidget( parent),
		is_log_(false),
		show_legend_(true),
		alignment_(alignment),
		inverse_orientation_(false),
		margin_(0),
		legend_(legend),
		tick_level_(3),
		pen_width_(0),
		allow_short_numbers_(false)
	{		
		if (!(alignment==RIGHT || alignment==LEFT || alignment==BOTTOM || alignment==TOP))
		{
			throw Exception::OutOfRange(__FILE__, __LINE__, __PRETTY_FUNCTION__);
		}
		
		setAxisBounds(0.0,100.0);
		
	 	if (alignment==RIGHT || alignment==LEFT )	
	 	{
			setMinimumSize(30,100);
			setSizePolicy(QSizePolicy::Fixed,QSizePolicy::MinimumExpanding);
		}
		else
		{
			setMinimumSize(100, 30);
			setSizePolicy(QSizePolicy::MinimumExpanding,QSizePolicy::Fixed);
		}
		resize(minimumSize());
	}
	
	
	AxisWidget::~AxisWidget()
	{
	}
	
	void AxisWidget::paintEvent(QPaintEvent *)
	{
		const int legend_default_size = 12;
		const int tick1_default_size = 16;

		QColor text_color;
		int tick_size = 0;
		
		QPainter painter(this);
	 	
		bool isXAxis = (alignment_==BOTTOM || alignment_==TOP);
		int h = height();
		int w = width();
		w = (isXAxis)? w-margin_ : w;  // Remove margin to get the scale right
		h = (isXAxis)? h : h-margin_;
	
		// probe legend font size (use a third (fifth) of available space for legend)
		int legend_font_size(0);
		if (isXAxis)
		{
			legend_font_size = legend_default_size;
		}
		else
		{
			legend_font_size = legend_default_size;
		}
		
		// probe text font size with last axis labels  
		QString probe = "";
		int index = 0;
		for (unsigned int i=0; i<3; i++)
		{
			QString s;
			if (grid_line_.size() > i)
			{
				getShortenedNumber_(s, scale_(grid_line_[i][ grid_line_[i].size() - 1]));
			}
			if (s.length() > probe.length())
			{
				probe = s;
				index = i;
			}
		}
	
		int text_font_size(0);
		// use rest of space for axis label (i.e. 1/5 legend, 4/5 label) or use all available space if no legend is shown;
		if (show_legend_)
		{
			text_font_size = tick1_default_size;
		}
		else
		{
			text_font_size = tick1_default_size;
		}
		
		// length of one big tick equals a quarter of the text_font_size
		double grid_scaling = text_font_size/4;
		if (grid_scaling > w/5 ) grid_scaling = w/5;  // don't let font size get to big
		if (grid_scaling > h/5 ) grid_scaling = h/5;  
	
		for (unsigned int j = 0; j != min((unsigned int)grid_line_.size(), (unsigned int)2) ; j++) 
		{
	    if (is_log_ && j>0) break; // just draw text on big intervalls
			// style definitions
			switch(j)
			{
				case 0:	// style settings for big intervals 
					painter.setFont(QFont("courier", static_cast<UInt>(3*grid_scaling))); 
					tick_size = static_cast<UInt>( grid_scaling );
					text_color = QColor(0, 0, 0);				
				break;
				case 1:	// style settings for small intervals
					painter.setFont(QFont("courier",static_cast<UInt>(2.4*grid_scaling))); 
					tick_size = static_cast<int>( grid_scaling*2/3 );
					text_color = QColor(20, 20, 20);
				break;
				case 2: // style settings for smalles intervals
					painter.setFont(QFont("courier",static_cast<UInt>(2.0*grid_scaling))); 
					tick_size = static_cast<UInt>( grid_scaling/2);
					text_color = QColor(40, 40, 40);
				break;
				default:
					std::cerr << "empty grid line vector error! in Line: " << __LINE__ << " File: " << __FILE__ << std::endl;
					painter.setFont(QFont("courier", static_cast<UInt>(3*grid_scaling))); 
					tick_size = static_cast<UInt>( grid_scaling);
					text_color = QColor(0, 0, 0);				
				break;
			}
	
			// drawing
			int x;
			UInt i_beg = (isXAxis)? 0 : h;
			UInt i_end = (isXAxis)? w : 0;
	
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
				painter.setPen(QPen(Qt::black, pen_width_));
				switch (alignment_)	
				{
					case BOTTOM: painter.drawLine(x, 0, x, tick_size); break;
					case TOP: painter.drawLine(x, h, x,  h-tick_size); break;
					case LEFT: painter.drawLine(w-tick_size, x+margin_, w, x+margin_); break;
					case RIGHT: painter.drawLine(0, x+margin_, tick_size, x+margin_);	break;
				}			
					
				// values at axis lines
				QString text;
				getShortenedNumber_(text, scale_(*it));
				painter.setPen(QPen(text_color, pen_width_));
		
				// get bounding rectangle for text we want to layout
				QRect textbound = QFontMetrics(painter.font()).boundingRect(text);
	
				// Calculate text position
				UInt posx = 0, posy = 0;
				switch (alignment_)	
				{
				  case BOTTOM: case TOP:	posx = x - static_cast<UInt>( textbound.width()/2 ); break;  // Center text around tick
			  	case LEFT: posx = w - static_cast<UInt>(1.5*grid_scaling) - textbound.width(); break;
				  case RIGHT: posx = static_cast<UInt>(1.5*grid_scaling); break; // Leave space for tick (longest one = grid_scaling)
				}	
	
				switch (alignment_)	
				{
				  case BOTTOM: posy = static_cast<UInt>(1.5*grid_scaling) + static_cast<UInt>(textbound.height()/2.0); break;
				  case TOP:	posy = h-static_cast<UInt>(1.5*grid_scaling); break;
				  case LEFT: case RIGHT: posy = x + margin_ + static_cast<UInt>(textbound.height()/4.0);  break;
				}	
				painter.drawText(posx, posy, text);
			}		
		}
	
		// style settings for legend
		painter.setBrush(Qt::NoBrush);
		painter.setPen(QPen(Qt::black,pen_width_));
		
		// text label at axis
		if (show_legend_ && legend_!="") 
		{	
			// style settings
			painter.setFont(QFont("courier", legend_font_size));
			painter.setPen(QPen(Qt::black,pen_width_));
	
			switch (alignment_)	
			{
	    		case BOTTOM:
						painter.drawText(0, 0 ,  w, h, Qt::AlignBottom|Qt::AlignHCenter, legend_.c_str());	
						break;
				  case TOP:
						painter.drawText(0, 0 ,  w, h, Qt::AlignTop|Qt::AlignHCenter, legend_.c_str());	
						break;
				  case LEFT: 
						painter.rotate(270);
	 				painter.drawText(-h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignTop, legend_.c_str());
						break;
				  case RIGHT:
						painter.rotate(270);
						painter.drawText(-h, 0 ,h ,w, Qt::AlignHCenter|Qt::AlignBottom, legend_.c_str());
						break;
			}	
		}
	}
	
	
	void AxisWidget::setAxisBounds(double min, double max)
	{
	  if (min >= max)
	  {
	    return;
	  }
	  
		if (is_log_)
		{
			//abort if no change
			if (min_== linear2log(min) && max_ == linear2log(max)) return;
			min_ = linear2log(min);
			max_ = linear2log(max);
			AxisTickCalculator::calcLogGridLines(min_,max_,grid_line_);
		}
		else
		{
			//abort if no change
			if (min_==min && max_==max) return;
			min_ = min; 
			max_ = max;
			UInt max_num_big = (alignment_==BOTTOM || alignment_==TOP)? 7: 5;
			UInt max_num_small = (alignment_==BOTTOM || alignment_==TOP)? 5: 3;
			AxisTickCalculator::calcGridLines(min_, max_,tick_level_,grid_line_,max_num_big, max_num_small);
		}
		update();
	}
	
	void AxisWidget::setLogScale(bool is_log)
	{
	  if (is_log_ != is_log)
	  {
	    is_log_ = is_log;
			
			if (is_log_)
			{
				setAxisBounds(min_,max_);
			}
			else 
			{
				setAxisBounds(log2linear(min_),log2linear(max_));
			}
		  update();
	  }
	}
	
	bool AxisWidget::isLogScale()
	{
		return is_log_;
	}
	
	UInt AxisWidget::margin()
	{
		return margin_;
	}
	
	void AxisWidget::setMargin(UInt margin)
	{
	  if (margin_ != margin){
	    margin_ = margin;
		  update();
	  }
	}
	
	
	void AxisWidget::showLegend(bool show_legend)
	{
		if (show_legend_ != show_legend)
		{
	  	show_legend_ = show_legend;
		  if (show_legend_)
		  { 
			  setToolTip("");
	  	}
	  	else
	  	{
		  	setToolTip(legend_.c_str());
	    }
	    update();
	  }
	}
	
	bool AxisWidget::isLegendShown() const
	{
		return show_legend_;
	}
	
	const String& AxisWidget::getLegend()
	{
		return legend_;
	}
	
	void AxisWidget::setInverseOrientation(bool inverse_orientation)
	{
		if (inverse_orientation_ != inverse_orientation)
		{
	  	inverse_orientation_ = inverse_orientation;
	    update();
	  }
	}
	
	bool AxisWidget::hasInverseOrientation()
	{
		return inverse_orientation_;
	}
	
	void AxisWidget::setLegend(const String& legend)
	{
		if (legend_ != legend)
		{
			legend_ = legend;
			if (!show_legend_)
	  	{ 
			  setToolTip(legend_.c_str());
   	 	}
  	}
	}
	
	void AxisWidget::setTickLevel(UInt level)
	{
		if (level>=1 && level<=3)
		{
			tick_level_ = level;
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

	void AxisWidget::getShortenedNumber_(QString& short_num, double number)
	{
		if (!allow_short_numbers_)
		{
			short_num = QString("%1").arg(number);
			return;
		}
		if (number < 1000.0)
		{
			short_num = QString("%1").arg(number);
			return;
		}
		if (number < 1000000.0)
		{
			short_num = QString("%1k").arg(Math::round_decimal(number/1000.0, -1));
			return;
		}
		if (number < 1000000000.0)
		{
			short_num = QString("%1M").arg(number/1000000.0);
			return;
		}
		short_num = QString("%1G").arg(number/1000000000.0);
		return;
	}

	void AxisWidget::setAllowShortNumbers(bool short_nums)
	{
		allow_short_numbers_ = short_nums;
	}

} //Namespace



