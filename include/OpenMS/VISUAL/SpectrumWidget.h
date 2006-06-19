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
// $Id: SpectrumWidget.h,v 1.26 2006/06/08 15:51:32 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUMWIDGET_H
#define OPENMS_VISUAL_SPECTRUMWIDGET_H

//OpenMS
#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/VISUAL/AxisWidget.h>

//STL

// QT
class QPopupMenu;
class QAction;
class QGridLayout;

namespace OpenMS
{
	class AxisWidget;
	class SpectrumCanvas;
	class SpectrumWindow;
	
	/**
		@brief Base class for spectrum widgets
		
		This class is the base class for widgets which contain a
		widget derived from SpectrumCanvas. For each class derived
		from SpectrumCanvas, there must be another class derived
		from this class.
		
		A SpectrumWidget instance has two instances of AxisWidget
		as children, which can be accessed by xAxis() and
		yAxis().
		
		@todo document, remove unneeded methods (Marc)
		
		@ingroup spectrum_widgets
	*/
	class SpectrumWidget : public QWidget, public PreferencesManager
	{
		Q_OBJECT
	
	public:
		/**
			@brief Returns image of the diagram.
			
			Creates an image of the diagram with the specified size and possibly with special setting (e.g. for printing or saving).
			
			@param width The image's width.
			@param height The image's height.
			@param flags Image creation flags.
			@returns The created image.
		*/
		virtual QImage getImage(UnsignedInt width, UnsignedInt height);
		
		/**
			@brief Returns pointer to canvas object
			
			Returns a pointer to the canvas object. The canvas object
			is set with the setCanvas() method. This is usually done
			in the constructor.
			
			@note Don't delete this object
			
			@return the canvas widget
		*/
		inline SpectrumCanvas* canvas() { return canvas_; }
		
		/**
			@brief Returns pointer to x-axis widget
			
			Returns a pointer to the x-axis axis widget.
			
			@note Don't delete this object
			
			@return the x-axis widget
		*/
		inline AxisWidget* xAxis() 
		{ 
			return x_axis_; 
		}
		
		/**
			@brief Returns pointer to y-axis widget
			
			Returns a pointer to the y-axis axis widget.
			
			@note Don't delete this object
			
			@return the y-axis widget
		*/
		inline AxisWidget* yAxis() 
		{ 
			return y_axis_; 
		}
		
		/**
			@brief Sets the spectrum window
			
			Sets the pointer to the SpectrumWindow instance which
			contains this widget. This is neccessary for integration
			with the TOPPView application. You can get the pointer
			to the SpectrumWindow instance with the getSpectrumWindow()
			method.
			
			@param window The SpectrumWindow instance, can be 0
		*/
		inline void setSpectrumWindow(SpectrumWindow* window) { spectrum_window_ = window; }
		
		/**
			@brief Returns pointer to spectrum window
			
			Returns the pointer to the SpectrumWindow instance
			containing this widget, as set with setSpectrumWindow().
			
			@note Don't delete this object
			@note The return value can be 0
			
			@return the Spectrum
		*/
		inline SpectrumWindow* getSpectrumWindow() const { return spectrum_window_; }
		
		///Get the mouse action mode
		SignedInt getActionMode() const;
		
		///Set the main Param object
		void setMainPreferences(const Param& prefs);

		///PreferencesManager
		virtual PreferencesDialogPage* createPreferences(QWidget* parent)=0;
		
		bool isLogIntensity() const;
		inline bool getShowLegend() const { return show_legend_; }
		
	signals:
		///signals that draw or display mode changed (e.g. used to update the tool bar)
		void modesChanged(QWidget*);
		/// Displays a status message. See SpectrumMDIWindow::showStatusMessage .
		void sendStatusMessage(std::string, OpenMS::UnsignedInt);
		void sendCursorStatus(double,double,double);
		void contextMenu(QPoint pos);
		
	public slots:
		void actionModeSelect();
		void actionModeZoom();
		void actionModeTranslate();
		void actionModeMeasure();
		void setActionMode(QAction* a);
		void setActionMode(OpenMS::SignedInt mode);
		void setIntensityModificationNone();
		void setIntensityModificationLog();
		void showIntensityDistribution();
		void showLegend();
		void showNoLegend();
		void setMirroredXAxis(bool b);
		void setMirroredYAxis(bool b);

		/**
			@brief Sets whether grid lines are shown or not.
			
			@param show Boolean variable deciding whether or not to show the grid lines.
		*/
		void showGridLines(bool show);

		virtual void switchAxis(bool swapped_axes);
		
	protected:
		/// Default constructor
		SpectrumWidget(QWidget* parent = 0, const char* name="SpectrumWidget", WFlags f=0);
		/// Destructor
		~SpectrumWidget();
		
		/// Adds the canvas to the layout and connects some signals/slots
		void setCanvas(SpectrumCanvas* canvas);
  	
  	/// Switches between log/normal intensities
  	virtual void intensityModificationChange_() = 0;
  	
  	/// Shows/hides the axis units
		virtual void legendModificationChange_() = 0;

		/// creates the intensity distribution of the widget
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_() = 0;
		
		/// recalculates the Axis ticks
		virtual void recalculateAxes() = 0;
		
		/// Canvas widget
		SpectrumCanvas* canvas_;
		
		///Main layout
		QGridLayout* grid_;
		
		/// Vertical axis
		AxisWidget* y_axis_;
		/// Horizontal axis
		AxisWidget* x_axis_;
		/// Spacer for the horizontal axis
		QWidget* hspacer_;
		/// Spacer for the vertical axis
		QWidget* vspacer_;
		
		/// Flag that indicates if the axis legend is shown
		bool show_legend_;
		
		///for storing the old maximum when the intensities are transformed
		double old_max_intensity_;
	
		
		SpectrumWindow* spectrum_window_;	
	
	private slots:
		/// updates the axes, when the visible area changes
		void updateAxes_(DRange<2> area);
	};
}

#endif
