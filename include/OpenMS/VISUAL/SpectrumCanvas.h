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


#ifndef OPENMS_VISUAL_SPECTRUMCANVAS_H
#define OPENMS_VISUAL_SPECTRUMCANVAS_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/PreferencesManager.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/VISUAL/LayerData.h>

//QT
#include <QtGui/QWidget>
#include <QtGui/QRubberBand>

//STL
#include <stack>
#include <vector>

namespace OpenMS
{
	class SpectrumWidget;
	class DataReducer;
	
	/**
		@brief Base class for visualization canvas classes
		
		This class is the base class for the spectrum data views. It
		derives from QScrollView, so scrollbars are provided. The viewing
		area is also managed by this class.
		
		It also provides commonly used constants such as ActionModes or IntensityModes.
		
		To provide additional spectrum views, you can derive from this class.
		You should also create a subclass from SpectrumWidget which encloses
		your class derived from SpectrumCanvas. To integrate you class into
		TOPPView, you also need to derive from SpectrumWindow.
		
		@ingroup spectrum_widgets
	*/
	class SpectrumCanvas 
		: public QWidget, 
			public PreferencesManager
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
		/// Spectrum iterator type (iterates over peaks)
		typedef SpectrumType::Iterator SpectrumIteratorType;
		/// Peak type
		typedef SpectrumType::PeakType PeakType;
		/// Feature type
		typedef FeatureMapType::FeatureType FeatureType;

    ///Type of the Points
		typedef DPosition<2> PointType;
		///Types of Ranges/Areas
		typedef DRange<2> AreaType;
		
		
		///Mouse action modes
		enum ActionModes 
		{
			AM_SELECT,		///< select a peaks
			AM_ZOOM,			///< zoom in / out
			AM_TRANSLATE,	///< move the visible area
			AM_MEASURE		///< measure distance between peaks
		};
		
		///Display modes of intensity
		enum IntensityModes
		{
			IM_NONE,		    ///< Normal mode: f(x)=x
			IM_LOG,			    ///< Log mode: f(x)=ln(x)
			IM_PERCENTAGE,  ///< Shows intensities normalized by layer maximum: f(x)=x/max(x)*100
			IM_SNAP         ///< Shows the maximum displayed intensity as if it was the overall maximum intensity
		};
		
		//@}

		/// Default constructor
		SpectrumCanvas(QWidget* parent = 0);
		
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
			
			@see actionModeChange_()
		*/
		inline void setActionMode(ActionModes mode) 
		{ 
			action_mode_ = mode;
			if (mode == AM_ZOOM)
			{
				setCursor(Qt::CrossCursor);
			}
			else if (mode == AM_SELECT)
			{
			  setCursor(Qt::ArrowCursor);
			}
			else if (mode == AM_TRANSLATE)
			{
				setCursor(cursor_translate_);	
			}
			else if (mode == AM_MEASURE)
			{
				setCursor(Qt::ArrowCursor);
			}
			actionModeChange_();
		}
		
		/**
			@brief Returns the intensity mode
			
			Returns the current intensity mode of type IntensityModes
			
			@return the current intensity mode
		*/
		inline SignedInt getIntensityMode() const 
		{ 
			return intensity_mode_; 
		}
		
		/**
			@brief Sets the intensity mode
			
			Sets the intensity mode
			
			@param mod the new intensity mode
			
			@see intensityModeChange_()
		*/
		inline void setIntensityMode(IntensityModes mod) 
		{
			intensity_mode_ = mod;
			intensityModeChange_();
		}
		
		/**
			@brief Returns if the grid is currently shown
			
			@return @c true if the grid is visible, @c false otherwise
		*/
		inline bool gridLinesShown() const 
		{ 
			return show_grid_; 
		}
		/// returns the layer data with index @p index
		inline const LayerData& getLayer(UnsignedInt index) const
		{
			OPENMS_PRECONDITION(index < layers_.size(), "SpectrumCanvas::getLayer() index overflow");
			return layers_[index];
		}
		
		/// returns the layer data of the active layer
		inline const LayerData& getCurrentLayer() const
		{
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::getLayer() index overflow");
			return layers_[current_layer_];
		}
		
		/**
			@brief Returns the currently visible area
			
			Dimension 0 is the m/z dimension.
			Dimension 1 is the RT dimension (not used in 1D).
		*/
		inline const AreaType& getVisibleArea() 
		{ 
			return visible_area_;
		}
		
		/**
			@brief Sets the minimum and maximum displayed intensities
			
			Sets the range of intensities. To be visible, a peak's intensity must be inside this range.
			@param min the minimum displayed intensity
			@param max the maximum displayed intensity
		*/
		virtual void setDispInt(float min, float max);
		
		/// Returns the mapping of m/z to axes
		bool isMzToXAxis();
		
		/// Sets the mapping of m/z to axes
		void mzToXAxis(bool mz_to_x_axis);
		
		/**
			@brief Sets the preferences object
			
			Sets the preferences object for this spectrum view
			@param prefs the preferences object
		*/
		virtual void setMainPreferences(const Param& prefs);		

		/** 
			@name Dataset handling methods
			
			@see changeVisibility
		*/
		//@{
		/// Returns the number of layers
		inline UnsignedInt getLayerCount() const
		{
			return layers_.size();
		}
		/// Returns the peak data (reduced or normal) of the @p index'th layer (not mutable)
		inline const ExperimentType& getPeakData(UnsignedInt index) const 
		{
			if(show_reduced_)
			{
				return getLayer(index).reduced;
			}
			return getLayer(index).peaks;
		}
		/// Returns the peak data (reduced or normal) of the active layer (not mutable)
		inline const ExperimentType& getCurrentPeakData() const
		{
			if(show_reduced_)
			{
				return getLayer(current_layer_).reduced;
			}
			return getLayer(current_layer_).peaks;
		}
			
		/// Returns the index of the active layer
		UnsignedInt activeLayerIndex() const;
		///change the active layer (the one that is used for selecting and so on)
		virtual void activateLayer(int layer_index)=0;
		///removes the layer with index @p layer_index
		virtual void removeLayer(int layer_index)=0;
		/**
    	@brief Adds another peak layer to fill afterwards
    	
    	Call finishAdding(float) after you filled the layer!
    	
    	@return reference to the new layer
    */
    ExperimentType& addEmptyPeakLayer();
		/**
			@brief Finish adding data after call to addEmptyPeakLayer()
		
			You can use this method instead of addLayer (add by copy).
			First call addEmptyPeakLayer(),then fill returned reference and finally call finishAdding(float).
		
			@return the index of the new layer
		*/
		virtual SignedInt finishAdding(float low_intensity_cutoff = 0) = 0;
		/**
			@brief Add a peak data layer (data is copied)
		
			@return the index of the new layer. -1 if no new layer was created.
		*/
		SignedInt addLayer(const ExperimentType&);

		/**
			@brief Add a feature data layer (data is copied)
			
			@param pairs Flag that indicates that a feature pair file was read.
			
			@return the index of the new layer. -1 if no new layer was created.
		*/
		SignedInt addLayer(const FeatureMapType& map, bool pairs);
		
		//@}
		
		/// Returns the minimum intensity of the active layer
		inline double getCurrentMinIntensity() const 
		{ 
			if (getCurrentLayer().type==LayerData::DT_PEAK)
			{
				return getCurrentPeakData().getMinInt(); 
			}
			else
			{
				return getCurrentLayer().features.getMinInt(); 
			}
		}

		/// Returns the maximum intensity of the active layer
		inline double getCurrentMaxIntensity() const 
		{ 
			if (getCurrentLayer().type==LayerData::DT_PEAK)
			{
				return getCurrentPeakData().getMaxInt(); 
			}
			else
			{
				return getCurrentLayer().features.getMaxInt(); 
			}
		}

		/// Returns the minimum intensity of the layer with index @p index
		inline double getMinIntensity(UnsignedInt index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK)
			{
				return getCurrentPeakData().getMinInt(); 
			}
			else
			{
				return getLayer(index).features.getMinInt(); 
			}
		}

		/// Returns the name of the layer with index @p index
		inline void setCurrentLayerName(const String& name) 
		{ 
		  getCurrentLayer_().name = name; 
		}

		/// Returns the maximum intensity of the active layer
		inline double getMaxIntensity(UnsignedInt index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK)
			{
				return getPeakData(index).getMaxInt(); 
			}
			else
			{
				return getLayer(index).features.getMaxInt(); 
			}
		}

		/**
			@brief Returns the area which encloses all data points.
			
			The order domensions is dependent on the derived class.
		*/
		const DRange<3>& getDataRange();	

		/**
			@brief Returns repaints the whole widget.
			
			Call this method after you changed the settings through the pulic interface of PreferencesManager
			in order to notify the widget of the changes.
		*/
		virtual void repaintAll();	

		/**
			@brief Returns the intensity scaling factor for 'snap to maximum intensity mode'.
			
			@see snap_factor_
		*/
		double getSnapFactor();
		
	public slots:
		/**
			@brief change the visibility of a layer
		
			@param i the index of the layer
			@param b true if layer is supposed to be visible
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
			@brief Notifies the canvas that the horizontal scollbar has been moved.
		
			Reimplement this slot to react on scrollbar events.
		*/
		virtual void horizontalScrollBarChange(int value);

		/**
			@brief Notifies the canvas that the vertical scollbar has been moved.
		
			Reimplement this slot to react on scrollbar events.
		*/
		virtual void verticalScrollBarChange(int value);
		
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

		/// Displays a status message. See TOPPViewBase::showStatusMessage .
		void sendStatusMessage(std::string, OpenMS::UnsignedInt);
			
		/// Forces recalculation of axis ticks in the connected widget.
		void recalculateAxes();
		
		/// Triggers the update of the vertical scrollbar
		void updateVScrollbar(float,float,float,float);

		/// Triggers the update of the horizontal scrollbar
		void updateHScrollbar(float,float,float,float);
	protected:

		inline LayerData& getLayer_(UnsignedInt index)
		{
			OPENMS_PRECONDITION(index < layers_.size(), "SpectrumCanvas::getLayer() index overflow");
			return layers_[index];
		}

		inline LayerData& getCurrentLayer_()
		{
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::getLayer() index overflow");
			return getLayer_(current_layer_);
		}

		/// Returns the @p index'th layer (mutable)
		inline ExperimentType& getPeakData_(UnsignedInt index)
		{
			if(show_reduced_)
			{
				return getLayer_(index).reduced;
			}
			return getLayer_(index).peaks;
		}
		
		/// Returns the active layer (mutable)
		inline ExperimentType& currentPeakData_()
		{
			if(show_reduced_)
			{
				return getCurrentLayer_().reduced;
			}
			return getCurrentLayer_().peaks;
		}
	
		/// reimplemented QT event
		void resizeEvent(QResizeEvent* e);
		
		/**
			@brief Change of the intensity distribution
			
			This function is called whenever the intensity distribution changes. Reimplement if you need to react on such changes.
		*/
		virtual void intensityDistributionChange_();

		///This function is called whenever the intensity mode changes. Reimplement if you need to react on such changes.
		virtual void intensityModeChange_();

		///This function is called whenever the action mode changes. Reimplement if you need to react on such changes.
		virtual void actionModeChange_();

		///This function is called whenever the action mode changes. Reimplement if you need to react on such changes.
		virtual void axisMappingChange_();
		
		/**
			@brief Sets the visible area
			
			Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
			If parts of the area are outside of the data area, the new area will be adjusted.
			
			@param new_area The new visible area.
			@param add_to_stack If the new area is to add to the zoom_stack_
		*/
		virtual void changeVisibleArea_(const AreaType& new_area, bool add_to_stack = false);

		/**
			@brief REcalculates the intensity scaling factor for 'snap to maximum intensity mode'.
			
			@see snap_factor_
		*/
		virtual void recalculateSnapFactor_();
		
		/**
			@brief Go back in zoom history
			
			Pops the topmost item from the zoom stack and resets the visible area to the area of that item
		*/
		void zoomBack_();
		
		/**
			@brief Updates the scroll bars
			
			Updates the scrollbars after a change of the visible area.
		*/
		virtual void updateScrollbars_();
		
		/**
			@brief Convert widget to chart coordinates
			
			Translates widget coordinates to chart coordinates.
			@param x the widget coordinate x
			@param y the widget coordinate y
			@return chart coordinates
		*/
		inline PointType widgetToData_(float x, float y)
		{
			if (!isMzToXAxis())
			{
				return PointType(
								visible_area_.minX() + (height() - y) / height()  * visible_area_.width(),
								visible_area_.minY() + x  / width() * visible_area_.height() 
								);
			}
			else
			{
				return PointType(
								visible_area_.minX() + x / width() * visible_area_.width(),
								visible_area_.minY() + (height() - y) / height() * visible_area_.height() 
								);			
			}		
		}

		/// Calls widgetToData_(float, float) with x and y position of @p pos
		inline PointType widgetToData_(const QPoint& pos)
		{
			return widgetToData_(pos.x(), pos.y());
		}
				
		/**
			@brief Convert chart to widget coordinates
			
			Translates chart coordinates to widget coordinates.
			@param x the chart coordinate x
			@param y the chart coordinate y
			@param point returned widget coordinates
		*/
		inline void dataToWidget_(float x, float y, QPoint& point)
		{
			if (!isMzToXAxis())
			{
				point.setX( static_cast<int>((y - visible_area_.minY()) / visible_area_.height() * width()));
				point.setY(height() - static_cast<int>((x - visible_area_.minX()) / visible_area_.width() * height()));
			}
			else
			{
				point.setX( static_cast<int>((x - visible_area_.minX()) / visible_area_.width() * width()));
				point.setY( height() - static_cast<int>((y - visible_area_.minY()) / visible_area_.height() * height()));
			}
		}
		
		/// Calls dataToWidget_(float x, float y, QPoint& point) with x and y position of @p pos
		inline void dataToWidget_(const PointType& pos, QPoint& point)
		{
			dataToWidget_(pos.X(), pos.Y(),point);
		}
		
		/**
			@brief Paints grid lines
			
			Helper function to paint grid lines
			
			@param p the QPainter to paint the grid lines on
		*/
		void paintGridLines_(QPainter& painter);
		
		/// Buffer that stores the actual peak information
		QPixmap buffer_;
		
		/// Stores the current action mode (Pick, Zoom, Translate)
		ActionModes action_mode_;
		
		/// Stores the used intensity mode function
		IntensityModes intensity_mode_;
		
		/// Layer data
		std::vector< LayerData > layers_;
		
		/// Stores the mapping of m/z
		bool mz_to_x_axis_;
		
		/// Stores the currently visible area.
		AreaType visible_area_;

		/**
			@brief Updates data and intensity range with the values of layer @p layer_index
			
			@param layer_index layer index
			@param mz_dim Index of m/z in overall_data_range_
			@param rt_dim Index of RT in overall_data_range_			
			@param it_dim Index of intensity in overall_data_range_	
			
			@see overall_data_range_
			
			@note Make sure the updateRanges() of the layers has been called before this method is called
		*/
		void updateRanges_(UnsignedInt layer_index, UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim);
		
		/**
			@brief Recalculates the data range.
			
			This method resets overall_data_range_ and calls updateRanges_(UnsignedInt,UnsignedInt,UnsignedInt,UnsignedInt)
			for all layers.
	
			@param mz_dim Index of m/z in overall_data_range_
			@param rt_dim Index of RT in overall_data_range_		
			@param it_dim Index of intensity in overall_data_range_	
			
			@see overall_data_range_
		*/
		void recalculateRanges_(UnsignedInt mz_dim, UnsignedInt rt_dim, UnsignedInt it_dim);
		
		/// Stores the data range (m/z, RT and intensity) of all layers
		DRange<3> overall_data_range_;
		
		/// Stores whether or not to show a grid.
		bool show_grid_;
		
		/// Flag for reduced displaying of peak layers
		bool show_reduced_;
		
		/// The zoom stack. This is dealt with in the changeVisibleArea_() and zoomBack_() functions.
		std::stack<AreaType> zoom_stack_;

		/**
			@brief Updates the diplayed data
			
			The default implementation calls QQidget::update().
			
			This method is reimplemented in the 3D view to update the OpenGL widget.
			
			@param caller_name Name of the calling function (use __PRETTY_FUNCTION__).
		*/
		virtual void update_(const char* caller_name);

		/// Whether to recalculate the data in the buffer when repainting
		bool update_buffer_;

		/// The cursor used in the @c translate action mode
		QCursor cursor_translate_;

		/// The cursor used in while the view is dragged in the @c translate action mode
		QCursor cursor_translate_in_progress_;

		/// Stores the index of the currently active layer.
		UnsignedInt current_layer_;

		/// Changes the size of the paint buffer to the currently required size
		void adjustBuffer_();
		
		/// Back-pointer to the enclosing spectrum widget
		SpectrumWidget* spectrum_widget_;

		/// pointer to the used datareducer implementation
		DataReducer * datareducer_;
		
		/// start position of mouse actions
		QPoint last_mouse_pos_;
		
		/**
			@brief Intensity scaling factor for relative scale with multiple layers.
			
			In this mode all layers are scaled to the same maximum.
		*/
		double percentage_factor_;
		
		/**
			@brief Intensity scaling factor for 'snap to maximum intensity mode'.
			
			In this mode the highest currently visible intensisty is treated like the maximum overall intensity.
		*/
		double snap_factor_;
		
		/// Rubber band for selected area
		QRubberBand rubber_band_;
	};
}

#endif
