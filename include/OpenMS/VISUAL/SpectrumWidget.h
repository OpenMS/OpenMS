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

#ifndef OPENMS_VISUAL_SPECTRUMWIDGET_H
#define OPENMS_VISUAL_SPECTRUMWIDGET_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QGridLayout;
class QScrollBar;

namespace OpenMS
{

	class AxisWidget;
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
		
		@ingroup spectrum_widgets
	*/
	class SpectrumWidget 
		: public QWidget,
			public PreferencesManager
	{
		Q_OBJECT
	
	public:
		/**
			@brief Returns image of the diagram.
			
			Creates an image of the diagram with the specified size and possibly with special setting (e.g. for printing or saving).
			
			@param width The image's width.
			@param height The image's height.
			@returns The created image.
		*/
		virtual QImage getImage(UnsignedInt width, UnsignedInt height);
		
		/**
			@brief Returns a pointer to canvas object
			
			Returns a pointer to the canvas object. The canvas object
			is set with the setCanvas_() method. This is usually done
			in the constructor.
			
			@note Don't delete this object
			
			@return the canvas widget
		*/
		SpectrumCanvas* canvas()
		{ 
			return canvas_; 
		}
		
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

		/// SpectrumWidgetActionModes
		void setActionMode(SpectrumCanvas::ActionModes mode);
		
		///Set the main Param object
		void setMainPreferences(const Param& prefs);

		// Docu in base class
		virtual PreferencesDialogPage* createPreferences(QWidget* parent)=0;
		
		/// Returns if the axis labels are shown
		bool isLegendShown() const;

		/// Returns if the axis labels are shown
		void showLegend(bool show);

		/// Sets the intensity mode of the SpectrumCanvas
		void setIntensityMode(SpectrumCanvas::IntensityModes mode);

		/// Hides x-axis and y-axis
		void hideAxes();
		
	signals:
		/// Signals that draw mode or display mode changed (e.g. used to update the tool bar)
		void modesChanged(QWidget*);
		/// Displays a status message. See TOPPViewBase::showStatusMessage .
		void sendStatusMessage(std::string, OpenMS::UnsignedInt);
		/// Displays peak information in the status bar (m/z, RT, intensity)
		void sendCursorStatus(double,double,double);
		
	public slots:
		/// Shows the intensity distribution of the data
		void showIntensityDistribution();
		/// Sets mapping of m/z values to x-axis or y-axis
		virtual void mzToXAxis(bool mz_to_x_axis);
		/// Updates the axes by calling recalculateAxes_();
		void updateAxes();
		/**
			@brief Updates the horizontal scrollbar
			
			@param min The overall minimum of the range
			@param disp_min The displayed minimum 
			@param disp_max The displayed maximum
			@param max The overall maximum of the range
		*/
		void updateHScrollbar(float min, float disp_min, float disp_max, float max);
		/**
			@brief Updates the vertical scrollbar
			
			@param min The overall minimum of the range
			@param disp_min The displayed minimum 
			@param disp_max The displayed maximum
			@param max The overall maximum of the range
		*/
		void updateVScrollbar(float min, float disp_min, float disp_max, float max);
	protected:
		/// Default constructor
		SpectrumWidget(QWidget* parent = 0);
		/// Destructor
		~SpectrumWidget();
		/// Adds the canvas to the layout and connects some signals/slots
		void setCanvas_(SpectrumCanvas* canvas);
  	/// Switch between different intensitiy modes
  	virtual void intensityModeChange_();
		/// creates the intensity distribution of the widget
		virtual Math::Histogram<UnsignedInt,float> createIntensityDistribution_() = 0;
		/// recalculates the Axis ticks
		virtual void recalculateAxes_() = 0;
		
		/// Pointer to the canvas widget
		SpectrumCanvas* canvas_;
		/// Pointer to the main window widget
		SpectrumWindow* spectrum_window_;	
		///Main layout
		QGridLayout* grid_;
		/// Vertical axis
		AxisWidget* y_axis_;
		/// Horizontal axis
		AxisWidget* x_axis_;
		/// Horizontal scrollbar
		QScrollBar* x_scrollbar_;
		/// Vertical scrollbar
		QScrollBar* y_scrollbar_;
	};
}

#endif
