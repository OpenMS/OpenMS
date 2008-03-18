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


#ifndef OPENMS_VISUAL_SPECTRUMCANVAS_H
#define OPENMS_VISUAL_SPECTRUMCANVAS_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

//QT
#include <QtGui/QWidget>
#include <QtGui/QRubberBand>

//STL
#include <stack>
#include <vector>

namespace OpenMS
{
	class SpectrumWidget;
	
	/**
		@brief Base class for visualization canvas classes
		
		This class is the base class for the spectrum data views.
		
		It also provides commonly used constants such as ActionModes or IntensityModes.
		
		To provide additional spectrum views, you can derive from this class.
		You should also create a subclass from SpectrumWidget which encloses
		your class derived from SpectrumCanvas. To integrate your class into
		TOPPView, you also need to derive a class from SpectrumWidget.

		@ingroup SpectrumWidgets
	*/
	class SpectrumCanvas 
		: public QWidget,
			public DefaultParamHandler
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
			AM_SELECT,		///< select + measure
			AM_ZOOM 			///< zoom + translate
		};
		
		///Display modes of intensity
		enum IntensityModes
		{
			IM_NONE,		    ///< Normal mode: f(x)=x
			IM_PERCENTAGE,  ///< Shows intensities normalized by layer maximum: f(x)=x/max(x)*100
			IM_SNAP         ///< Shows the maximum displayed intensity as if it was the overall maximum intensity
		};
		
		//@}

		/// Default constructor
		SpectrumCanvas(const Param& preferences, QWidget* parent = 0);

		/// Destructor
		~SpectrumCanvas();
		
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
		inline Int getActionMode() const 
		{ 
			return action_mode_;
		}
		
		/**
			@brief Sets the action mode
			
			Sets the action mode for the left mouse button, e.g. zoom, translate etc.
			@param mode the new action mode.
		*/
		inline void setActionMode(ActionModes mode) 
		{ 
			action_mode_ = mode;
			switch (mode)
			{
				case AM_ZOOM:
					setCursor(Qt::CrossCursor);
					break;
				default:
					setCursor(Qt::ArrowCursor);
			}
		}
		
		/**
			@brief Returns the intensity mode
			
			Returns the current intensity mode of type IntensityModes
			
			@return the current intensity mode
		*/
		inline Int getIntensityMode() const 
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
		inline const LayerData& getLayer(UInt index) const
		{
			OPENMS_PRECONDITION(index < layers_.size(), "SpectrumCanvas::getLayer(index) index overflow");
			return layers_[index];
		}
		
		/// returns the layer data of the active layer
		inline const LayerData& getCurrentLayer() const
		{
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::getCurrentLayer() index overflow");
			return layers_[current_layer_];
		}

		/// returns a layer flag of the current layer
		bool getLayerFlag(LayerData::Flags f) const
		{
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::getLayerFlag() index overflow");
			switch(f)
			{
				case LayerData::F_HULLS:
					return layers_[current_layer_].f1;
				case LayerData::F_NUMBERS:
					return layers_[current_layer_].f2;
				case LayerData::F_HULL:
					return layers_[current_layer_].f3;
				case LayerData::P_PRECURSORS:
					return layers_[current_layer_].f1;
				case LayerData::P_PROJECTIONS:
					return layers_[current_layer_].f2;
			}
			std::cout << "Error: SpectrumCanvas::getLayerFlag -- unknown flag '" << f << "'!" << std::endl;
			return false;
		}

		/// sets a layer flag of the current layer
		void setLayerFlag(LayerData::Flags f, bool value)
		{
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::setLayerFlag() index overflow");
			switch(f)
			{
				case LayerData::F_HULLS:
					layers_[current_layer_].f1 = value;
					break;
				case LayerData::F_NUMBERS:
					layers_[current_layer_].f2 = value;
					break;
				case LayerData::F_HULL:
					layers_[current_layer_].f3 = value;
					break;
				case LayerData::P_PRECURSORS:
					layers_[current_layer_].f1 = value;
					break;
				case LayerData::P_PROJECTIONS:
					layers_[current_layer_].f2 = value;
					break;
			}
			update_buffer_ = true;
			update();
		}

		/// returns a layer flag of the layer @p layer
		bool getLayerFlag(UInt layer, LayerData::Flags f) const
		{
			OPENMS_PRECONDITION(layer < layers_.size(), "SpectrumCanvas::getLayerFlag() index overflow");
			switch(f)
			{
				case LayerData::F_HULLS:
					return layers_[layer].f1;
				case LayerData::F_NUMBERS:
					return layers_[layer].f2;
				case LayerData::F_HULL:
					return layers_[layer].f3;
				case LayerData::P_PRECURSORS:
					return layers_[layer].f1;
				case LayerData::P_PROJECTIONS:
					return layers_[layer].f2;
			}
			std::cout << "Error: SpectrumCanvas::getLayerFlag -- unknown flag '" << f << "'!" << std::endl;
			return false;
		}

		/// sets a layer flag of the layer @p layer
		void setLayerFlag(UInt layer, LayerData::Flags f, bool value)
		{
			OPENMS_PRECONDITION(layer < layers_.size(), "SpectrumCanvas::setLayerFlag() index overflow");
			switch(f)
			{
				case LayerData::F_HULLS:
					layers_[layer].f1 = value;
					break;
				case LayerData::F_NUMBERS:
					layers_[layer].f2 = value;
					break;
				case LayerData::F_HULL:
					layers_[layer].f3 = value;
					break;
				case LayerData::P_PRECURSORS:
					layers_[layer].f1 = value;
					break;
				case LayerData::P_PROJECTIONS:
					layers_[layer].f2 = value;
					break;
			}
			update_buffer_ = true;
			update();
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
			@brief Sets the filters applied to the data before drawing (for the current layer)
		*/
		virtual void setFilters(const DataFilters& filters);
		
		/// Returns the mapping of m/z to axes
		bool isMzToXAxis();
		
		/// Sets the mapping of m/z to axes
		void mzToXAxis(bool mz_to_x_axis);

		/** 
			@name Dataset handling methods
			
			@see changeVisibility
		*/
		//@{
		/// Returns the number of layers
		inline UInt getLayerCount() const
		{
			return layers_.size();
		}
			
		/// Returns the index of the active layer
		UInt activeLayerIndex() const;
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
		virtual Int finishAdding() = 0;
		/**
			@brief Add a peak data layer (data is copied)
		
			@return the index of the new layer. -1 if no new layer was created.
		*/
		Int addLayer(const ExperimentType&);

		/**
			@brief Add a feature data layer (data is copied)
			
			@param pairs Flag that indicates that a feature pair file was read.
			@param map Feature map
			
			@return the index of the new layer. -1 if no new layer was created.
		*/
		Int addLayer(const FeatureMapType& map, bool pairs);
		
		//@}
		
		/// Returns the minimum intensity of the active layer
		inline double getCurrentMinIntensity() const 
		{ 
			if (getCurrentLayer().type==LayerData::DT_PEAK)
			{
				return getCurrentLayer().peaks.getMinInt(); 
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
				return getCurrentLayer().peaks.getMaxInt(); 
			}
			else
			{
				return getCurrentLayer().features.getMaxInt(); 
			}
		}

		/// Returns the minimum intensity of the layer with index @p index
		inline double getMinIntensity(UInt index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK)
			{
				return getCurrentLayer().peaks.getMinInt(); 
			}
			else
			{
				return getLayer(index).features.getMinInt(); 
			}
		}

		/// Sets the @p name of layer @p i
		void setLayerName(UInt i, const String& name);

		/// Sets the parameters of the current layer
		inline void setCurrentLayerParameters(const Param& param) 
		{ 
		  getCurrentLayer_().param = param;
		  currentLayerParamtersChanged_();
		}

		/// Returns the maximum intensity of the layer with index @p index
		inline double getMaxIntensity(UInt index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK)
			{
				return getLayer(index).peaks.getMaxInt(); 
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
			@brief Returns the intensity scaling factor for 'snap to maximum intensity mode'.
			
			@see snap_factor_
		*/
		double getSnapFactor();
		
		/// Shows the preferences dialog of the active layer
		virtual void showCurrentLayerPreferences() = 0;
		
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
		void sendCursorStatus(double mz=-1.0, double intens=-1.0, double rt=-1.0);

		/// Displays a status message. See TOPPViewBase::showStatusMessage .
		void sendStatusMessage(std::string, OpenMS::UInt);
			
		/// Forces recalculation of axis ticks in the connected widget.
		void recalculateAxes();
		
		/// Triggers the update of the vertical scrollbar
		void updateVScrollbar(float,float,float,float);

		/// Triggers the update of the horizontal scrollbar
		void updateHScrollbar(float,float,float,float);
		
		/// Toggle axis legend visibility change
		void changeLegendVisibility();

	protected:
		inline LayerData& getLayer_(UInt index)
		{
			OPENMS_PRECONDITION(index < layers_.size(), "SpectrumCanvas::getLayer_(index) index overflow");
			return layers_[index];
		}

		inline LayerData& getCurrentLayer_()
		{
			return getLayer_(current_layer_);
		}
		
		/// Returns the active layer (mutable)
		inline ExperimentType& currentPeakData_()
		{
			return getCurrentLayer_().peaks;
		}
	
		///reimplemented QT event
		void resizeEvent(QResizeEvent* e);
		
		/**
			@brief Change of layer parameters
			
			This method is called whenever the paramters of the current layer change. Reimplement if you need to react on such changes.
		*/
		virtual void currentLayerParamtersChanged_();

		///This method is called whenever the intensity mode changes. Reimplement if you need to react on such changes.
		virtual void intensityModeChange_();
		
		/**
			@brief Saves the current layer data.
			
			@param visible If true, only the visible data is stored. Otherwise the whole data is stored.
		*/
		virtual void saveCurrentLayer(bool visible)=0;
		
		
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
		
		///Helper function to paint grid lines
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
			@brief Recalculates the data range.
			
			A small margin is added to each side of the range in order to display all data.
	
			@param mz_dim Int of m/z in overall_data_range_
			@param rt_dim Int of RT in overall_data_range_		
			@param it_dim Int of intensity in overall_data_range_	
			
			@see overall_data_range_
		*/
		void recalculateRanges_(UInt mz_dim, UInt rt_dim, UInt it_dim);
		
		/// Stores the data range (m/z, RT and intensity) of all layers
		DRange<3> overall_data_range_;
		
		/// Stores whether or not to show a grid.
		bool show_grid_;
		
		/// The zoom stack. This is dealt with in the changeVisibleArea_() and zoomBack_() functions.
		std::stack<AreaType> zoom_stack_;

		/**
			@brief Updates the diplayed data
			
			The default implementation calls QWidget::update().
			
			This method is reimplemented in the 3D view to update the OpenGL widget.
			
			@param caller_name Name of the calling function (use __PRETTY_FUNCTION__).
		*/
		virtual void update_(const char* caller_name);

		/// Whether to recalculate the data in the buffer when repainting
		bool update_buffer_;

		/// Stores the index of the currently active layer.
		UInt current_layer_;

		/// Changes the size of the paint buffer to the currently required size
		void adjustBuffer_();
		
		/// Back-pointer to the enclosing spectrum widget
		SpectrumWidget* spectrum_widget_;

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
