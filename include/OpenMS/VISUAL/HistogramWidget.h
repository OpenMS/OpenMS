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
// $Id: HistogramWidget.h,v 1.7 2006/06/09 22:00:08 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_HISTOGRAMWIDGET_H
#define OPENMS_VISUAL_HISTOGRAMWIDGET_H

// QT
#include <qwidget.h>
#include <qpixmap.h>

// STL
#include <map>
#include <string>

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/AxisWidget.h>

namespace OpenMS
{
	/**
		@brief Widget which can visualize a histogram.
		
		It can also be used to define a left and right boundary inside the values.
		
		@todo splitters disappear sometimes (Marc)
		
		@ingroup Visual
	*/
	class HistogramWidget : public QWidget
	{
		Q_OBJECT
		public:

		HistogramWidget(const Histogram<UnsignedInt,float>& distribution, QWidget* parent = 0, const char* name = "HistogramWidget");
		virtual ~HistogramWidget();

		double getLeftSplitter();
		double getRightSplitter();

		/// set x-axis legend; default is 'x'
		void setLegend(const std::string& legend);

		public slots:
		void showSplitters(bool on);
		void setRightSplitter(double pos);
		void setLeftSplitter(double pos);

		protected:
		Histogram<UnsignedInt,float> dist_;
		bool show_splitters_;
		double left_splitter_;
		double right_splitter_;
		/// the splitter that is currently dragged (0=none, 1=left, 2=right)
		UnsignedInt moving_splitter_;
		AxisWidget *bottom_axis_;
		UnsignedInt margin_;
		/// internal buffer for the double buffering
		QPixmap* buffer_;
		/// repaints the contents to the buffer and calls update()
		void invalidate_();

		void paintEvent( QPaintEvent * );
		void mousePressEvent( QMouseEvent *);
		void mouseReleaseEvent( QMouseEvent *);
		void mouseMoveEvent( QMouseEvent *);
		void resizeEvent( QResizeEvent *);
	};
} // namespace OpenMS


#endif

