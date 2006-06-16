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
// $Id: SpectrumCanvas.h,v 1.35 2006/06/08 14:29:18 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_SPECTRUMCANVAS_H
#define OPENMS_VISUAL_SPECTRUMCANVAS_H

//OpenMS
#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/MappingInfo.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>

//STL
#include <stack>
#include <vector>

//QT
#include <qscrollview.h>
#include <qcursor.h>

class QPainter;
class QPixmap;

namespace OpenMS
{
	class AxisWidget;
	class SpectrumWidget;
	
	/**
		@brief Base class for visualization canvas classes
		
		This class is the base class for the spectrum data views. It
		derives from QScrollView, so scrollbars are provided. The viewing
		area is also managed by this class.
		
		It also provides commonly used constants such as ActionModes or
		IntensityModifications.
		
		To provide additional spectrum views, you can derive from this class.
		You should also create a subclass from SpectrumWidget which encloses
		your class derived from SpectrumCanvas. To integrate you class into
		TOPPView, you also need to derive from SpectrumWindow.
		
		@ingroup spectrum_widgets
	*/
	class SpectrumCanvas : public QScrollView, public PreferencesManager
	{
		Q_OBJECT
		
	public:
		/**	@name Type definitions */
		//@{
		
		/// Main data type (experiment)
		typedef MSExperiment<> ExperimentType;
		/// Spectrum type
		typedef ExperimentType::SpectrumType SpectrumType;
		/// Spectrum iterator type (iterates over peaks)
		typedef SpectrumType::Iterator SpectrumIteratorType;
		/// Peak type
		typedef SpectrumType::PeakType PeakType;

    ///Type of the Points
		typedef DPosition<2> PointType;
		///Types of Ranges/Areas
		typedef DRange<2> AreaType;
		
		
		///All action modes for the mouse
		enum ActionModes 
		{
			AM_SELECT,		///< select a peaks
			AM_ZOOM,			///< zoom in / out
			AM_TRANSLATE,	///< move the visible area
			AM_MEASURE		///< measure distance between peaks
		};
		
		///All modifications of intensity
		enum IntensityModifications
		{
			IM_NONE,		///f(x)=x
			IM_LOG			///f(x)=ln(x)
		};
		
		//@}

		/**
			brief Constructor.
			@param parent the parent QWidget
			@param name the widget's name
			@param f Window flags
		*/
		SpectrumCanvas(QWidget* parent = 0, const char* name="SpectrumCanvas", WFlags f=0);
		
		/**
			@brief Sets the spectrum widget.
			
			Sets the enclosing spectrum widget. Call this from your
			SpectrumWidget derived class.
			@param widget the spectrum widget
		*/
		inline void setSpectrumWidget(SpectrumWidget* widget) 
		{ 
			spectrum_widget_ = widget; 
		}
		
		/**
			@brief Returns the spectrum widget.
			
			Returns the enclosing spectrum widget
			@return the spectrum widget
		*/
		inline SpectrumWidget* getSpectrumWidget() const 
		{ 
			return spectrum_widget_; 
		}
		
		/**
			@brief Returns the action mode
			
			Returns the current action mode of type ActionModes
			@return the current action mode
		*/
		inline SignedInt getActionMode() const 
		{ 
			return action_mode_;
		}
		
		/**
			@brief Sets the action mode
			
			Sets the action mode for the left mouse button, e.g. zoom, translate etc.
			@param mode the new action mode.
			@see ActionModes
		*/
		inline void setActionMode(ActionModes mode) 
		{ 
			action_mode_ = mode; 
		}
		
		/**
			@brief Returns the intensity modification
			
			Returns the current intensity modification of type IntensityModifications
			
			@return the current intensity modification
		*/
		inline SignedInt getIntensityModification() const 
		{ 
			return intensity_modification_; 
		}
		
		/**
			@brief Sets the intensity modification
			
			Sets the intensity modification
			
			@param mod the new intensity modification
		*/
		inline void setIntensityModification(IntensityModifications mod) 
		{
			intensity_modification_ = mod;
			intensityModificationChange_();
		}
		
		/**
			@brief Returns the grid mode
			
			Returns if the grid is currently shown
			@return @c true if the grid is visible, @c false otherwise
		*/
		inline bool getGridMode() const 
		{ 
			return show_grid_; 
		}
		
		/**
			@brief Returns the minimum displayed intensity for the current layer
			
			Returns the minimum intensity a peak needs to be shown
			@return the minimum displayed intensity
		*/
		inline double getMinDispInt() const 
		{ 
			return min_disp_ints_[current_data_]; 
		}
		
		/**
			@brief Returns the maximum displayed intensity for the current layer
			
			Returns the maximum intensity a peak can have to be shown
			@return the maximum displayed intensity
		*/
		inline double getMaxDispInt() const 
		{ 
			return max_disp_ints_[current_data_];
		}
		
		/**
			@brief Returns the visible area
			
			Returns the current visible area.
			@return the visible area
		*/
		inline AreaType getVisibleArea() 
		{ 
			return visible_area_;
		}
		
		/**
			@brief Sets the minimum and maximum displayed intensities
			
			Sets the range of intensities. To be visible, a peak's intensity must be inside this range.
			@param min the minimum displayed intensity
			@param max the maximum displayed intensity
		*/
		virtual void setDispInt(double min, double max);
		
		/**
			@brief Returns the mapping info
			
			Returns a reference to the mapping info object
			@return the mapping info
			@see MappingInfo
		*/
		inline MappingInfo& getMappingInfo() 
		{ 
			return mapping_info_; 
		}
		
		/**
			@brief Sets the pen width
			
			Sets the pen width and repaints the widget. This is useful for printing.
			@param p the new pen width
		*/
		inline void setPenWidth(int p)
		{
			pen_width_ = p;
			invalidate_();
		}
		
		/**
			@brief Sets the preferences object
			
			Sets the preferences object for this spectrum view
			@param prefs the preferences object
		*/
		virtual void setMainPreferences(const Param& prefs);		

		/** 
			@name Dataset handling methods
			
			@note see changeVisibility(int,bool) as well.
		*/
		//@{	
		
		/// Returns the number of datasets
		UnsignedInt getDataSetCount() const;
		/// Returns the @p index'th dataset (mutable)
		ExperimentType& getDataSet(UnsignedInt index) throw (Exception::IndexOverflow);
		/// Returns the @p index'th dataset (not mutable)
		const ExperimentType& getDataSet(UnsignedInt index) const throw (Exception::IndexOverflow);
		/// Returns the active dataset (mutable)
		ExperimentType& currentDataSet() throw (Exception::IndexOverflow);
		/// Returns the active dataset (not mutable)
		const ExperimentType& currentDataSet() const throw (Exception::IndexOverflow);
		/// Returns the name associated with dataset @p index
		const String& getDataSetName(UnsignedInt index) const;		
		/// Returns true if dataset @p index is visible. false otherwise
		bool isDataSetVisible(UnsignedInt index) const;			
		/// Returns the index of the active dataset
		UnsignedInt activeDataSetIndex() const;
		///change the active spectrum (the one that is used for selecting and so on)
		virtual void activateDataSet(int data_set)=0;
		///removes the dataset with index @p data_set
		virtual void removeDataSet(int data_set)=0;
		/**
    	@brief Adds another dataset to fill afterwards
    	
    	Call finishAdding() after you filled the dataset!
    	
    	@return reference to the new dataset
    */
    ExperimentType& addEmptyDataSet();
		/**
			@brief Finish adding data after call to addEmtpyDataSet()
		
			You can use this method instead of addDataSet (add by copy).
			First call addEmptyDataSet(),then fill returned reference and finally call finishAdding().
		
			@return the index of the new dataset
		*/
		virtual SignedInt finishAdding() = 0;
		/**
			@brief Add another dataset by copy
		
			@return the index of the new dataset. -1 if no new dataset was created.
		*/
		SignedInt addDataSet(const ExperimentType&);
		
		//@}
		
		/// Returns the minimum intensity of the active spectrum
		inline double getCurrentMinIntensity() const 
		{ 
			return currentDataSet().getMinInt(); 
		}

		/// Returns the maximum intensity of the active spectrum
		inline double getCurrentMaxIntensity() const 
		{ 
			return currentDataSet().getMaxInt(); 
		}
	
		
	public slots:
		/**
			@brief change the visibility of a spectrum
		
			@param i the index of the dataset
			@param b true if spectrum is supposed to be visible
		*/
		void changeVisibility(int i, bool b);

		/**
			@brief Whether or not to show grid lines
			
			Sets whether grid lines are shown or not.
			@param show Boolean variable deciding whether or not to show the grid lines.
		*/
		void showGridLines(bool show);
		
		/**
			@brief Zooms fully out.
			
			Sets the visible area to the initial value, such that all data is shown.
		*/
		void resetZoom();
		
		/**
			@brief Sets the visible area.
			
			Sets the visible area to a new value. Note that it does not emit visibleAreaChanged()
			@param area the new visible area
		*/
		void setVisibleArea(AreaType area);
		
		/**
			@brief Whether or not the X axis should be mirrored
			
			Sets whether the X axis should be the other way around
			@param b if @c true, the X axis is mirrored, if @c falseit is shown normally
		*/
		void setMirroredXAxis(bool b);
		
		/**
			@brief Whether or not the Y axis should be mirrored
			
			Sets whether the Y axis should be the other way around
			@param b if @c true, the Y axis is mirrored, if @c falseit is shown normally
		*/
		void setMirroredYAxis(bool b);
		
	signals:
		/// Signal emitted whenever a new Layer is activated within the current window
		void layerActivated(QWidget* w);
		
		/**
			@brief Change of the visible area
			
			Signal emitted whenever the visible area changes.
			@param area The new visible area.
		*/
		void visibleAreaChanged(DRange<2> area); //Do not change this to AreaType! QT needs the exact type...
				
		/// Emited when the cursor position changes (for displaying in status bar)
		void sendCursorStatus(double pos=-1.0, double intens=-1.0, double rt=-1.0);
		
		/**
			@brief Context menu request
			
			This signal is emitted when a context menu has to be shown, e.g. after the right mouse button is hit
			@param pos The global position where the context menu
			should pop up
		*/
		void contextMenu(QPoint pos);

		/// Displays a status message. See SpectrumMDIWindow::showStatusMessage .
		void sendStatusMessage(std::string, OpenMS::UnsignedInt);

	private slots:
		/// Updates the scroll bars
		void move_(int x, int y);
		/// Sets zoom_timeout_ to true. Connected to a timer
		void timeoutZoom_();
		
	protected:
		
		/**
			@brief QT resize event of the widget
			
		*/
		virtual void viewportResizeEvent(QResizeEvent* e);
		
		/**
			@brief QT repaint event of the widget
			
		*/
		virtual void viewportPaintEvent(QPaintEvent* e);
		
		/**
			@brief Change of the intensity distribution
			
			This function is called whenever the intensity distribution changes. Reimplement if you need to react on such changes.
		*/
		virtual void intensityDistributionChange_();

		/**
			@brief Change of the intensity modifiaction
			
			This function is called whenever the intensity modification changes. Reimplement if you need to react on such changes.
		*/
		virtual void intensityModificationChange_();
		
		/**
			@brief Invalidates the contents of the back buffer and repaints.
			
			Repaints the content into buffer_ after a data or view change (e.g. zoom, translate, displayed intesity). You need to
			reimplement this method and you need to draw all contents into buffer_ rather than to paint on the widget directly.
			
			@see recalculate_
		*/
		virtual void invalidate_() = 0;
		
		/**
			@brief Sets the visible area
			
			Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
			If parts of the area are outside of the data area, the new area will be adjusted.
			@param new_area The new visible area.
		*/
		virtual void changeVisibleArea_(const AreaType& new_area);
		
		/**
			@brief Go back in zoom history
			
			Pops the topmost item from the zoom stack and resets the visible area to the area of that item
		*/
		void zoomBack_();
		
		/**
			@brief Updates the scroll bars
			
			Updates the scrollbars after a change of the visible area.
		*/
		void updateScrollbars_();
		
		/**
			@brief Returns the data area
			
			Returns the area which encloses all data.
		*/
		virtual const AreaType& getDataArea_() = 0;
		
		/**
			@brief Convert pixel to chart coordinates
			
			Translates pixel coordinates to chart coordinates with respect toa given window size
			@param p the pixel coordinates
			@param width the window width
			@param height the window height
			@return chart coordinates
		*/
		PointType contextToChart_(QPoint p, int width, int height);
		
		/**
			@brief Convert chart to pixel coordinates
			
			Translates chart coordinates to pixel coordinates with respect to a given window size
			@param pos the chart coordinates
			@param width the window width
			@param height the window height
			@return pixel coordinates
		*/
		QPoint chartToContext_(const PointType& pos, int width, int height);
		
		/**
			@brief Convert widget to chart coordinates
			
			Translates widget coordinates to chart coordinates.
			@param pos the widget coordinates
			@return chart coordinates
		*/
		PointType widgetToChart_(const QPoint& pos);
		
		/**
			@brief Convert chart to widget coordinates
			
			Translates chart coordinates to widget coordinates.
			@param pos the chart coordinates
			@return widget coordinates
		*/
		QPoint chartToWidget_(const PointType& pos);
		
		/**
			@brief Paints grid lines
			
			Helper function to paint grid lines for the spectrum view
			
			@param p the QPainter to paint the grid lines on
		*/
		void paintGridLines_(QPainter* p);
		
		/// Buffer pixmap for fast reblitting of damaged content
		QPixmap* buffer_;
		
		/// Stores the current action mode (Pick, Zoom, Translate)
		ActionModes action_mode_;
		
		/// Stores the used intensity modification function
		IntensityModifications intensity_modification_;
		
		/// Stores the minimum displayed intensities for all layers
		std::vector<float> min_disp_ints_;
		
		/// Stores the maximum displayed intensities for all layers
		std::vector<float> max_disp_ints_;
		
		/// Stores the mapping info
		MappingInfo mapping_info_;
		
		/// Stores the pen width. Drawing thicker lines (e.g. in printing) leads to better results
		UnsignedInt pen_width_;
		
		/// Stores the currently visible area.
		AreaType visible_area_;
		
		/// Stores whether or not to show a grid.
		bool show_grid_;
		
		/// The zoom stack. This is dealt with in the changeVisibleArea_() and zoomBack_() functions.
		std::stack<AreaType> zoom_stack_;

		/// Whether to recalculate the data before drawing. This is used to optimize redrawing in invalidate_().
		bool recalculate_;
		
		/**
			@brief Creates mouse cursors
			
			Creates custom cursors for translate action
		*/
		void createCustomMouseCursors_();

		/// The cursor used in the @c translate action mode
		QCursor cursor_translate_;

		/// The cursor used in while the view is dragged in the @c translate action mode
		QCursor cursor_translate_in_progress_;

		/// Stores if the current area was visible for >= 2 secs. If so, the area gets added to the zoom stack later on
		bool zoom_timeout_;

		/// Stores the index of the currently marked dataset in the layerbar.
		UnsignedInt current_data_;

		/// Stores for each layer, if it is shown
		std::vector<bool> layer_visible_;

		/// Changes the size of the paint buffer to the currently required size
		void adjustBuffer_();
		
		/// Back-pointer to the enclosing spectrum widget
		SpectrumWidget* spectrum_widget_;
		
		/// Array of datasets
		std::vector<ExperimentType > datasets_;
	};
}

#endif
