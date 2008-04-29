// -*- Mode: C++; tab-width: 2; -*-
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

#ifndef OPENMS_VISUAL_SPECTRUMWIDGET_H
#define OPENMS_VISUAL_SPECTRUMWIDGET_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QGridLayout;
class QScrollBar;

namespace OpenMS
{

	class AxisWidget;
	
	/**
		@brief Base class for spectrum widgets
		
		This class is the base class for the different MDI window
		types in the TOPPView application. For each type of spectrum
		view (such as 1D view, 2D view etc.), there must exist a
		corresponding class derived from this class.
		
		To integrate a new spectrum view (i.e. classes derived from
		SpectrumWidget and SpectrumCanvas) into the TOPPView application,
		a class must be derived from this class which holds an
		instance of the SpectrumCanvas class as a child widget.
		
		This Widget also provides axis widgets and scrollbars.
		
		@ingroup SpectrumWidgets
	*/
	class SpectrumWidget 
		: public QWidget
	{
		Q_OBJECT

		public:
			/// Default constructor
			SpectrumWidget(const Param& preferences, QWidget* parent = 0);
			/// Destructor
			virtual ~SpectrumWidget();
			
			/**
				@brief Returns a pointer to canvas object
				
				The canvas object is set with the setCanvas_() method. 
				This is usually done in the constructor.
			*/
			SpectrumCanvas* canvas()
			{ 
				return canvas_; 
			}
			
			///Returns a pointer to the x-axis axis widget.
			inline AxisWidget* xAxis() 
			{ 
				return x_axis_; 
			}
			
			///Returns a pointer to the y-axis axis widget.
			inline AxisWidget* yAxis() 
			{ 
				return y_axis_; 
			}
			
			///Get the mouse action mode
			Int getActionMode() const;
	
			/// SpectrumWidgetActionModes
			void setActionMode(SpectrumCanvas::ActionModes mode);
			
			/// Returns if the axis labels are shown
			virtual bool isLegendShown() const;
	
			/// Shows/hides axis labels
			virtual void showLegend(bool show);
	
			/// Sets the intensity mode of the SpectrumCanvas
			void setIntensityMode(SpectrumCanvas::IntensityModes mode);
	
			/// Hides x-axis and y-axis
			void hideAxes();
			
			/// Widget id used as identifier
			Int window_id;
			
		signals:
			/// Signals that draw mode or display mode changed (e.g. used to update the tool bar)
			void modesChanged(QWidget*);
			/// Displays a status message. See TOPPViewBase::showStatusMessage .
			void sendStatusMessage(std::string, OpenMS::UInt);
			/// Displays peak information in the status bar (m/z, RT, intensity)
			void sendCursorStatus(double mz=-1.0, double intens=-1.0, double rt=-1.0);
		  /// Message about the destruction of this widget
		  void aboutToBeDestroyed(int window_id);
		  /// Shows the main preferences dialog
		  void openPreferences();
			  
		public slots:
			/// Shows statistics about the data (count,min,max,avg of Intensity, Charge, Quality and meta data)
			void showStatistics();
			/// Shows the intensity distribution of the data
			void showIntensityDistribution();
			/// Updates the axes by setting the right labels and calling recalculateAxes_();
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
			/// Shows a goto dialog
			virtual void showGoToDialog() = 0;
			/// Toggles the axis legend visibility
			void changeLegendVisibility();
			
		protected:
			/**
				@brief Adds the canvas, axes and scrollbars to the layout
			
				@p row and @p col define the position of the canvas. 
				Axes and scrollbars are added to the left and bottom of the canvas.
			*/
			void setCanvas_(SpectrumCanvas* canvas, UInt row=0, UInt col=2);
	  	/// Switch between different intensitiy modes
	  	virtual void intensityModeChange_();
			/// creates the intensity distribution of the widget
			virtual Math::Histogram<UInt,float> createIntensityDistribution_() = 0;
			/// recalculates the Axis ticks
			virtual void recalculateAxes_() = 0;
			
			/// Pointer to the canvas widget
			SpectrumCanvas* canvas_;
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
