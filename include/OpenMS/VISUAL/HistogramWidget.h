// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_HISTOGRAMWIDGET_H
#define OPENMS_VISUAL_HISTOGRAMWIDGET_H

// QT
#include <QtGui/QWidget>
#include <QtGui/QPixmap>
class QPaintEvent;
class QResizeEvent;
class QMouseEvent;

//OpenMS
#include <OpenMS/MATH/STATISTICS/Histogram.h>

namespace OpenMS
{
	class AxisWidget;
	class String;
	
	/**
		@brief Widget which can visualize a histogram.
		
		@image html HistogramWidget.png
		
		It can also be used to define a left and right boundary inside the values.
    It supports normal and log scaling via the context menu.
		
		@ingroup Visual
	*/
	class OPENMS_GUI_DLLAPI HistogramWidget
		: public QWidget
	{
		Q_OBJECT
		
		public:
			/// Constructor
			HistogramWidget(const Math::Histogram<>& distribution, QWidget* parent = 0);

			/// Destructor
			virtual ~HistogramWidget();
			
			/// Returns the value f the lower splitter
			DoubleReal getLeftSplitter();

			/// Returns the value of the upper splitter
			DoubleReal getRightSplitter();

      /// Set axis legends
			void setLegend(const String& legend);

		public slots:
			/// Shows the splitters if @p on is true. Hides them otherwise.
			void showSplitters(bool on);

			/// Sets the value of the right splitter
			void setRightSplitter(DoubleReal pos);

			/// Sets the value of the left splitter
			void setLeftSplitter(DoubleReal pos);

			/// Enables/disables log mode
			void setLogMode(bool log_mode);

		protected:
      /// The histogram to display
			Math::Histogram<> dist_;

			/// Flag that indicates if splitters are shown
			bool show_splitters_;

      /// Value of the right splitter
			DoubleReal left_splitter_;

      /// Value of the right splitter
			DoubleReal right_splitter_;

      /// The splitter that is currently dragged (0=none, 1=left, 2=right)
			UInt moving_splitter_;

      /// X axis
			AxisWidget *bottom_axis_;

      /// Margin around plot
			UInt margin_;

      /// Internal buffer for the double buffering
			QPixmap buffer_;

      /// Flag that indicates the current mode
			bool log_mode_;

      /// Repaints the contents to the buffer and calls update()
			void invalidate_();

			///@name reimplemented Qt events
			//@{
			void paintEvent(QPaintEvent*);
			void mousePressEvent(QMouseEvent*);
			void mouseReleaseEvent(QMouseEvent*);
			void mouseMoveEvent(QMouseEvent*);
			void resizeEvent(QResizeEvent*);
			//@}

		protected slots:

      /// Context menu event
			void showContextMenu(const QPoint& pos);
	};
} // namespace OpenMS


#endif
