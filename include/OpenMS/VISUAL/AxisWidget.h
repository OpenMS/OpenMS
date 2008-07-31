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

#ifndef OPENMS_VISUAL_AXISWIDGET_H
#define OPENMS_VISUAL_AXISWIDGET_H

// QT
#include <QtGui/QWidget>
class QPaintEvent;

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/MATH/MISC/MathFunctions.h>

namespace OpenMS
{
	/**
		@brief Widget that represents an axis of a graph.
		
		Additional to ticks and tick values a label e.g. the unit can be displayed.
		It supports both linear and logarithmic scale.
		
		@image html AxisWidget.png
		The above image shows a horizontal example axis.
		
		@ingroup Visual
	*/
	class AxisWidget 
		: public QWidget
	{
		Q_OBJECT
		
		public:
	
			///Vector of vector of doubles that defines the grid
			typedef std::vector<std::vector<double> > GridVector;
			
			/// Where the axis is placed
			static enum {TOP, BOTTOM, LEFT, RIGHT} ALIGNMENT_ENUM;
	
			/// constructor
			AxisWidget(UInt alignment, const char* legend="", QWidget* parent = 0);
			
			/// destructor
			virtual ~AxisWidget();
	
			/// sets the margin on the top/right side (default is 0)
			void setMargin(UInt size);
	
			/// returns the margin
			UInt margin();
	
			/// enable the display of the legend (default true)
			void showLegend(bool show_legend);
	
			/// returns true if legend is shown
			bool isLegendShown() const;
	
			/**
				@brief return constant reference to grid_line_
				
				gridlines are calculated by setAxisBounds or set by the user using setGridLines
			 */
			const GridVector& gridLines();
	
			/// set user-defined gridlines (instead of using setAxisBounds to calculate GridLines)
			void setGridLines(std::vector<double>&);
	
			/// sets the axis to logarithmic scale
			void setLogScale(bool is_log);
	
			/// returns true if the axis has logarithmic scale
			bool isLogScale();
	
			/// sets the legend text
			void setLegend(const String& legend);
	
			/// returns the actual legend text
			const String& getLegend();

			/// set true to display the axis label in inverse order (left to right or bottom to top)
			void setInverseOrientation(bool inverse_orientation);
			
			/// returns if the axis label is displayed in inverse order
			bool hasInverseOrientation();
	
			/// set true to allow for shortened numbers (with k/M/G units) on the axis label
			void setAllowShortNumbers(bool short_nums = true);
			
			/// returns the minimum value displayed on the axis
	    inline double getAxisMinimum() const
	    {
	    	return min_;
	    }

	    /// returns the maximum value displayed on the axis
	    inline double getAxisMaximum() const 
	    {
	    	return max_;
	    }
	    
	    /// sets the pen width for both text and lines and calls update()
			inline void setPenWidth(int p)
			{
				pen_width_ = p;
				update();
			}
	
		public slots:
		
			///sets min/max of the axis
			void setAxisBounds(double min, double max);
			
			/// set maximum number of tick levels (1 <= level <= 3)
			void setTickLevel(UInt level);
	
		protected:
			/// Vector that defines the position of the ticks/gridlines and the shown values on axis
			GridVector grid_line_;
	
			/// format of axis scale (linear or logarithmic)
			bool is_log_;
			/// display of legend enabled or not
			bool show_legend_;
			/// Position of the axis (right, left, top, down as defined in ALIGNMENT_ENUM)
			UInt alignment_;
			/// true if axis label are displayed in inverse order (left to right or bottom to top)
			bool inverse_orientation_;
			/// margin of axis
			UInt margin_;
			/// minimum value on the axis
			double min_;
			/// maximum value on the axis
			double max_;
			/// text/unit on axis
			String legend_;
			/// maximum number of tick levels (default=3)
			UInt tick_level_;
			/// drawing thicker lines (e.g. in printing) leads to better results
			UInt pen_width_;

			/// see QWidget
			void paintEvent( QPaintEvent * );
	
			/// Scale axis values to correct value (i.e. reverse log, unit conversion)
			inline double scale_(double x)
			{
				return (is_log_)? Math::round_decimal(pow(x,10),-8) : Math::round_decimal(x,-8);
			}
			
			/// sets @p short_num to a shortened string representation ("123.4 k/M/G") of @p number
			void getShortenedNumber_(QString& short_num, double number);

			/// true if k/M/G units can be used
			bool allow_short_numbers_;
 	};
} // namespace OpenMS


#endif

