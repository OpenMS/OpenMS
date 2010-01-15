// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_SPECTRUMWIDGET_H
#define OPENMS_VISUAL_SPECTRUMWIDGET_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/MATH/STATISTICS/Histogram.h>
#include <OpenMS/VISUAL/SpectrumCanvas.h>

class QGridLayout;
class QScrollBar;
class QCloseEvent;
class QMimeData;

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
		
		@todo Add support to store the displayed data as SVG image (HiWi)
	*/
	class OPENMS_DLLAPI SpectrumWidget 
		: public QWidget
	{
		Q_OBJECT

		public:
			/**	@name Type definitions */
			//@{
			
			/// Main data type (experiment)
			typedef LayerData::ExperimentType ExperimentType;
			/// Main data type (features)
			typedef LayerData::FeatureMapType FeatureMapType;
			/// Spectrum type
			typedef ExperimentType::SpectrumType SpectrumType;
			//@}
			
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
	
			/// Returns if the axis labels are shown
			virtual bool isLegendShown() const;
	
			/// Shows/hides axis labels
			virtual void showLegend(bool show);
	
			/// Sets the intensity mode of the SpectrumCanvas
			void setIntensityMode(SpectrumCanvas::IntensityModes mode);
	
			/// Hides x-axis and y-axis
			virtual void hideAxes();
			
			/// Saves the widget's content as image file
			void saveAsImage();
			
			/// Widget id used as identifier
			Int window_id;
			
		signals:
			/// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
			void sendStatusMessage(std::string, OpenMS::UInt);
			/// Emitted when the cursor position changes (for displaying e.g. in status bar)
			void sendCursorStatus(double mz=-1.0, double rt=-1.0);
		  /// Message about the destruction of this widget
		  void aboutToBeDestroyed(int window_id);
		  /// Shows the main preferences dialog
		  void openPreferences();
		  /// Signal that is emitted, when a drag-and-drop action ends on this widget
			void dropReceived(const QMimeData* data, QWidget* source, int id);			
			  
		public slots:
			/// Shows statistics about the data (count, min, max, avg of intensity, charge, quality and meta data)
			void showStatistics();
			/// Shows the intensity distribution of the current layer
			void showIntensityDistribution();
			/// Shows the meta data distribution of vlaue @p name of the current layer
			void showMetaDistribution(const String& name);
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
			/// @name Reimplemented Qt events
			//@{ 
			void closeEvent(QCloseEvent* e);
			//@}
				
			/**
				@brief Adds the canvas, axes and scrollbars to the layout
			
				@p row and @p col define the position of the canvas. 
				Axes and scrollbars are added to the left and bottom of the canvas.
			*/
			void setCanvas_(SpectrumCanvas* canvas, UInt row=0, UInt col=2);
	  	/// Switch between different intensitiy modes
	  	virtual void intensityModeChange_();
			/// creates the intensity distribution of the current layer
			virtual Math::Histogram<> createIntensityDistribution_() const = 0;
			/// creates the meta data distribution of value @p name of the current layer
			virtual Math::Histogram<> createMetaDistribution_(const String& name) const = 0;
			/// recalculates the Axis ticks
			virtual void recalculateAxes_() = 0;
			
			///@name reimplemented Qt events
      //@{
     	void dragEnterEvent(QDragEnterEvent* event);
			void dragMoveEvent(QDragMoveEvent* event);
			void dropEvent(QDropEvent* event);
			//@}
			
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
