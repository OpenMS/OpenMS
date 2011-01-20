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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_SPECTRUMCANVAS_H
#define OPENMS_VISUAL_SPECTRUMCANVAS_H

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>

//QT
#include <QtGui/QWidget>
#include <QtGui/QRubberBand>
class QWheelEvent;
class QKeyEvent;
class QMouseEvent;
class QFocusEvent;
class QMenu;

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
		
		All derived classes should follow these interface conventions:
		- Translate mode
		  - Activated by default
		  - Arrow keys can be used to translate without entering translate mode
		- Zoom mode
		  - Activated using the CTRL key
		  - Zoom stack traversal with CTRL+/CTRL- or mouses wheel
		  - Pressing the @em Backspace key resets the zoom (and stack)
    - Measure mode
      - Activated using the SHIFT key
    
    @improvement Add log mode (Hiwi)

    @todo Allow reordering the layer list by drag-and-drop (Hiwi, Johannes)

		@htmlinclude OpenMS_SpectrumCanvas.parameters

		@ingroup SpectrumWidgets
	*/
	class OPENMS_DLLAPI SpectrumCanvas 
		: public QWidget,
			public DefaultParamHandler
	{
		Q_OBJECT
	
	public:
		/**	@name Type definitions */
		//@{
		
    /// Main data type (experiment)
    typedef LayerData::ExperimentType ExperimentType;
    /// Main managed data type (experiment)
    typedef LayerData::ExperimentSharedPtrType ExperimentSharedPtrType;
		/// Main data type (features)
		typedef LayerData::FeatureMapType FeatureMapType;
    /// Main managed data type (features)
    typedef LayerData::FeatureMapSharedPtrType FeatureMapSharedPtrType;
		/// Main data type (consensus features)
		typedef LayerData::ConsensusMapType ConsensusMapType;
    /// Main managed data type (consensus features)
    typedef LayerData::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

		/// Spectrum type
		typedef ExperimentType::SpectrumType SpectrumType;
		/// Spectrum iterator type (iterates over peaks)
    typedef SpectrumType::ConstIterator SpectrumConstIteratorType;
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
			AM_TRANSLATE, ///< translate
			AM_ZOOM, 			///< zoom
			AM_MEASURE    ///< measure
		};
		
		///Display modes of intensity
		enum IntensityModes
		{
			IM_NONE,		    ///< Normal mode: f(x)=x
			IM_PERCENTAGE,  ///< Shows intensities normalized by layer maximum: f(x)=x/max(x)*100
      IM_SNAP,        ///< Shows the maximum displayed intensity as if it was the overall maximum intensity
      IM_LOG          ///< Logarithmic mode
		};
		
		//@}

		/// Default constructor
		SpectrumCanvas(const Param& preferences, QWidget* parent = 0);

		/// Destructor
		virtual ~SpectrumCanvas();
		
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
		inline const LayerData& getLayer(Size index) const
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

    /// returns the layer data of the active layer
    inline LayerData& getCurrentLayer()
    {
      OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::getCurrentLayer() index overflow");
      return layers_[current_layer_];
    }

		/// returns a layer flag of the current layer
		bool getLayerFlag(LayerData::Flags f) const
		{
			return getLayerFlag(current_layer_, f);
		}

		/// sets a layer flag of the current layer
		void setLayerFlag(LayerData::Flags f, bool value)
		{
			setLayerFlag(current_layer_, f, value);
		}

		/// returns a layer flag of the layer @p layer
		bool getLayerFlag(Size layer, LayerData::Flags f) const
		{
			OPENMS_PRECONDITION(layer < layers_.size(), "SpectrumCanvas::getLayerFlag() index overflow");
			return layers_[layer].flags.test(f);
		}
		
		/// sets a layer flag of the layer @p layer
		void setLayerFlag(Size layer, LayerData::Flags f, bool value)
		{
			//abort if there are no layers
			if (layers_.empty()) return;
			OPENMS_PRECONDITION(layer < layers_.size(), "SpectrumCanvas::setLayerFlag() index overflow");
			
			layers_[layer].flags.set(f, value);
			update_buffer_ = true;
			update();
		}

		inline void setLabel(LayerData::LabelType label)
		{
			//abort if there are no layers
			if (layers_.empty()) return;
			
			OPENMS_PRECONDITION(current_layer_ < layers_.size(), "SpectrumCanvas::setLabel() index overflow");
			layers_[current_layer_].label = label;
			
			update_buffer_ = true;
			update();
		}
		
		/**
			@brief Returns the currently visible area
			
			@see visible_area_
		*/
		inline const AreaType& getVisibleArea() const
		{ 
			return visible_area_;
		}
		
		/**
			@brief Sets the filters applied to the data before drawing (for the current layer)
		*/
		virtual void setFilters(const DataFilters& filters);
		
		/// Returns the mapping of m/z to axes
		inline bool isMzToXAxis()
		{ 
			return mz_to_x_axis_; 
		}
		
		/// Sets the mapping of m/z to axes
		void mzToXAxis(bool mz_to_x_axis);

		/** 
			@name Dataset handling methods
			
			@see changeVisibility
		*/
		//@{
		/// Returns the number of layers
		inline Size getLayerCount() const
		{
			return layers_.size();
		}
			
		/// Returns the index of the active layer
		Size activeLayerIndex() const;
		///change the active layer (the one that is used for selecting and so on)
		virtual void activateLayer(Size layer_index)=0;
		///removes the layer with index @p layer_index
		virtual void removeLayer(Size layer_index)=0;
		/**
			@brief Add a peak data layer
				
			If chromatograms are present, a chromatogram layer is shown. Otherwise a peak layer is shown. Make sure to remove chromatograms from peak data and vice versa.
		
      @param map Shared Pointer to input map. It can be performed in constant time and does not double the required memory.
			@param filename This @em absolute filename is used to monitor changes in the file and reload the data

			@return If a new layer was created
		*/
    bool addLayer(ExperimentSharedPtrType map, const String& filename="");

		/**
			@brief Add a feature data layer
			
      @param map Shared Pointer to input map. It can be performed in constant time and does not double the required memory.
      @param filename This @em absolute filename is used to monitor changes in the file and reload the data
			
			@return If a new layer was created
		*/
    bool addLayer(FeatureMapSharedPtrType map, const String& filename="");
		
		/**
			@brief Add a consensus feature data layer
			
      @param map Shared Pointer to input map. It can be performed in constant time and does not double the required memory.
      @param filename This @em absolute filename is used to monitor changes in the file and reload the data
			
			@return If a new layer was created
		*/
    bool addLayer(ConsensusMapSharedPtrType map, const String& filename="");
		//@}
		
		/**
			@brief Add an identification data layer
			
			@param peptides Input list of peptides, which has to be mutable and will be empty after adding. Swapping is used to insert the data. It can be performed in constant time and does not double the required memory. 
			@param filename This @em absolute filename is used to monitor changes in the file and reload the data
			
			@return If a new layer was created
		*/
		bool addLayer(std::vector<PeptideIdentification>& peptides, 
									const String& filename="");
		
		/// Returns the minimum intensity of the active layer
		inline Real getCurrentMinIntensity() const 
		{ 
			if (getCurrentLayer().type==LayerData::DT_PEAK || getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
			{
        return getCurrentLayer().getPeakData()->getMinInt();
			}
			else if (getCurrentLayer().type==LayerData::DT_FEATURE)
			{
        return getCurrentLayer().getFeatureMap()->getMinInt();
			}
			else
			{
        return getCurrentLayer().getConsensusMap()->getMinInt();
			}
		}

		/// Returns the maximum intensity of the active layer
		inline Real getCurrentMaxIntensity() const 
		{ 
			if (getCurrentLayer().type==LayerData::DT_PEAK || getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
			{
        return getCurrentLayer().getPeakData()->getMaxInt();
			}
			else if (getCurrentLayer().type==LayerData::DT_FEATURE)
			{
        return getCurrentLayer().getFeatureMap()->getMaxInt();
			}
			else
			{
        return getCurrentLayer().getConsensusMap()->getMaxInt();
			}
		}

		/// Returns the minimum intensity of the layer with index @p index
		inline Real getMinIntensity(Size index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK || getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
			{
        return getLayer(index).getPeakData()->getMinInt();
			}
			else if (getLayer(index).type==LayerData::DT_FEATURE)
			{
        return getLayer(index).getFeatureMap()->getMinInt();
			}
			else
			{
        return getLayer(index).getConsensusMap()->getMinInt();
			}
		}

		/// Returns the maximum intensity of the layer with index @p index
		inline Real getMaxIntensity(Size index) const 
		{ 
			if (getLayer(index).type==LayerData::DT_PEAK || getCurrentLayer().type==LayerData::DT_CHROMATOGRAM)
			{
        return getLayer(index).getPeakData()->getMaxInt();
			}
			else if (getLayer(index).type==LayerData::DT_FEATURE)
			{
        return getLayer(index).getFeatureMap()->getMaxInt();
			}
			else
			{
        return getLayer(index).getConsensusMap()->getMaxInt();
			}
		}

		/// Sets the @p name of layer @p i
		void setLayerName(Size i, const String& name);

    /// Gets the name of layer @p i
    String getLayerName(Size i);

		/// Sets the parameters of the current layer
		inline void setCurrentLayerParameters(const Param& param) 
		{ 
		  getCurrentLayer_().param = param;
		  emit preferencesChange();
		}

		/**
			@brief Returns the area which encloses all data points.
			
			@see overall_data_range_
		*/
		const DRange<3>& getDataRange();	

		/**
			@brief Returns the first intensity scaling factor for 'snap to maximum intensity mode'.
			
			@see snap_factors_
		*/
		DoubleReal getSnapFactor();
		
		/// Returns the percentage factor
		DoubleReal getPercentageFactor();
		
		/// Shows the preferences dialog of the active layer
		virtual void showCurrentLayerPreferences() = 0;

		/**
			@brief Shows a dialog with the meta data
			
			@param modifiable indicates if the data can be modified.
			@param index If given, the meta data of the corresponding element (spectrum, feature, consensus feature) is shown instead of the layer meta data.
		*/
		virtual void showMetaData(bool modifiable=false, Int index = -1);

		/**
			@brief Saves the current layer data.
			
			@param visible If true, only the visible data is stored. Otherwise the whole data is stored.
		*/
		virtual void saveCurrentLayer(bool visible)=0;
		
	public slots:
		
		/**
			@brief change the visibility of a layer
		
			@param i the index of the layer
			@param b true if layer is supposed to be visible
		*/
		void changeVisibility(Size i, bool b);

		/**
			@brief change if the defined data filters are used
		
			@param i the index of the layer
			@param b true if layer is supposed to be visible
		*/
		void changeLayerFilterState(Size i, bool b);

		/**
			@brief Whether or not to show grid lines
			
			Sets whether grid lines are shown or not.
			@param show Boolean variable deciding whether or not to show the grid lines.
		*/
		void showGridLines(bool show);
		
		/**
			@brief Zooms fully out and resets the zoom stack
			
			Sets the visible area to the initial value, such that all data is shown.
			
			@param repaint If @em true a repaint is forced. Otherwise only the new area is set.
		*/
		void resetZoom(bool repaint = true);
		
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
		
		///Sets the additional context menu. If not 0, this menu is added to the context menu of the canvas
		void setAdditionalContextMenu(QMenu* menu);
		
		/**
			@brief Fills the handed over @p map with the visible peaks of the current layer. 
			
			Takes zoom area and data filters into account.
			
			If the current layer is not a peak layer, @p map is cleared only.
		*/
		void getVisiblePeakData(ExperimentType& map) const;


		/**
			@brief Fills the handed over @p map with the visible features of the current layer. 
			
			Takes zoom area and data filters into account.
			
			If the current layer is not a feature layer, @p map is cleared only.
		*/
		void getVisibleFeatureData(FeatureMapType& map) const;

		/**
			@brief Fills the handed over @p map with the visible consensus features of the current layer. 
			
			Takes zoom area and data filters into account.
			
			If the current layer is not a consensus feature layer, @p map is cleared only.
		*/
		void getVisibleConsensusData(ConsensusMapType& map) const;
		
		/**
			@brief Fills the handed over @p peptides with the visible peptide identifications of the current layer. 
			
			Takes zoom area into account.
			
			If the current layer is not an identification data layer, @p peptides is cleared only.
		*/
		void getVisibleIdentifications(std::vector<PeptideIdentification>& peptides) const;

    ///Updates layer @p i when the data in the corresponding file changes
    virtual void updateLayer(Size i) = 0;

	signals:

		/// Signal emitted whenever the modification status of a layer changes (editing and storing)
		void layerModficationChange(Size layer, bool modified);
		
		/// Signal emitted whenever a new layer is activated within the current window
		void layerActivated(QWidget* w);
		
		/**
			@brief Change of the visible area
			
			Signal emitted whenever the visible area changes.
			@param area The new visible area.
		*/
		void visibleAreaChanged(DRange<2> area); //Do not change this to AreaType! QT needs the exact type...
				
		/// Emitted when the cursor position changes (for displaying e.g. in status bar)
		void sendCursorStatus(double mz=-1.0, double rt=-1.0);

		/// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
		void sendStatusMessage(std::string message, OpenMS::UInt time);
			
		/// Forces recalculation of axis ticks in the connected widget.
		void recalculateAxes();
		
		/// Triggers the update of the vertical scrollbar
		void updateVScrollbar(float,float,float,float);

		/// Triggers the update of the horizontal scrollbar
		void updateHScrollbar(float,float,float,float);
		
		/// Toggle axis legend visibility change
		void changeLegendVisibility();
		
		/// Emitted when the action mode changes
		void actionModeChange();
		
		/// Emitted when the layer preferences have changed
		void preferencesChange();
		
	protected slots:
	
		///Updates the cursor accoring to the current action mode
		void updateCursor_();

	protected:

		/// Draws several lines of text to the upper right corner of the widget
		void drawText_(QPainter& painter, QStringList text);

		/// Returns the m/z value of an identification depending on the m/z source of the layer (precursor mass/theoretical peptide mass)
		DoubleReal getIdentificationMZ_(const Size layer_index, 
																		const PeptideIdentification& peptide) const;
	
		///Method that is called when a new layer has been added
		virtual bool finishAdding_() = 0;
		
		///Returns the layer with index @p index
		inline LayerData& getLayer_(Size index)
		{
			OPENMS_PRECONDITION(index < layers_.size(), "SpectrumCanvas::getLayer_(index) index overflow");
			return layers_[index];
		}

		///Returns the currently active layer
		inline LayerData& getCurrentLayer_()
		{
			return getLayer_(current_layer_);
		}
		
    /// Returns the currently active layer (mutable)
    inline ExperimentSharedPtrType currentPeakData_()
    {
            return getCurrentLayer_().getPeakData();
    }
	
		///@name reimplemented QT events
		//@{
		void resizeEvent(QResizeEvent* e);
		void wheelEvent(QWheelEvent* e);
		void keyPressEvent(QKeyEvent* e);
		void keyReleaseEvent(QKeyEvent* e);
		void focusOutEvent(QFocusEvent* e);
		void leaveEvent(QEvent* e);
		void enterEvent(QEvent* e);
		//@}
		
		///This method is called whenever the intensity mode changes. Reimplement if you need to react on such changes.
		virtual void intensityModeChange_();
				
		/**
			@brief Sets the visible area
			
			Changes the visible area, adjustes the zoom stack and notifies interested clients about the change. 
			If parts of the area are outside of the data area, the new area will be adjusted.
			
			@param new_area The new visible area.
			@param repaint If @em true, a complete repaint is forced.
			@param add_to_stack If @em true the new area is to add to the zoom_stack_.
		*/
		virtual void changeVisibleArea_(const AreaType& new_area, bool repaint = true, bool add_to_stack = false);

		/**
			@brief REcalculates the intensity scaling factor for 'snap to maximum intensity mode'.
			
			@see snap_factors_
		*/
		virtual void recalculateSnapFactor_();
		
		///@name Zoom stack methods
		//@{
		///Go backward in zoom history
		void zoomBack_();
		///Go forward in zoom history
		virtual void zoomForward_();
		/// Add a visible area to the zoom stack
		void zoomAdd_(const AreaType& area);
		/// Clears the zoom stack and invalidates the current zoom position. After calling this, a valid zoom position has to be added immediately.
		void zoomClear_();
		//@}
		
		///@name Translation methods, which are called when cursor buttons are pressed
		//@{
		/// Translation bound to the 'Left' key
		virtual void translateLeft_();
		/// Translation bound to the 'Rightt' key
		virtual void translateRight_();
		/// Translation bound to the 'Up' key
		virtual void translateForward_();
		/// Translation bound to the 'Down' key
		virtual void translateBackward_();
		//@}
		
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


        if (intensity_mode_ != SpectrumCanvas::IM_LOG)
        {
        point.setX( int((y - visible_area_.minY()) / visible_area_.height() * width()));
        } else  // IM_LOG
        {
          point.setX( int(
              std::log10((y - visible_area_.minY())+1) / std::log10(visible_area_.height()+1) * width())
                      );
        }

				point.setY(height() - int((x - visible_area_.minX()) / visible_area_.width() * height()));
			}
			else
			{
        point.setX( int((x - visible_area_.minX()) / visible_area_.width() * width()));

        if (intensity_mode_ != SpectrumCanvas::IM_LOG)
        {
          point.setY( height() - int((y - visible_area_.minY()) / visible_area_.height() * height()));
        } else  // IM_LOG
        {
          point.setY( height() - int(
              std::log10((y-visible_area_.minY())+1)/std::log10(visible_area_.height()+1)*height()
              ));
        }
			}
		}
		
		///Helper function to paint grid lines
		virtual void paintGridLines_(QPainter& painter);
		
		/// Buffer that stores the actual peak information
		QImage buffer_;
		
		/// Stores the current action mode (Pick, Zoom, Translate)
		ActionModes action_mode_;
		
		/// Stores the used intensity mode function
		IntensityModes intensity_mode_;
		
		/// Layer data
		std::vector< LayerData > layers_;
		
		/// Stores the mapping of m/z
		bool mz_to_x_axis_;
		
		/**
			@brief Stores the currently visible area.
			
			Dimension 0 is the m/z dimension.@n
			Dimension 1 is the RT dimension (2D and 3D view) or the intensity dimension (1D view).
		*/
		AreaType visible_area_;
		
		/**
			@brief Recalculates the overall_data_range_
			
			A small margin is added to each side of the range in order to display all data.
	
			@param mz_dim Int of m/z in overall_data_range_
			@param rt_dim Int of RT in overall_data_range_		
			@param it_dim Int of intensity in overall_data_range_	
		*/
		void recalculateRanges_(UInt mz_dim, UInt rt_dim, UInt it_dim);
		
		/**
			@brief Stores the data range (m/z, RT and intensity) of all layers
			
			Dimension 0 is the m/z dimension.@n
			Dimension 1 is the RT dimension (2D and 3D view) or the intensity dimension (1D view).@n
			Dimension 2 is the intensity dimension (2D and 3D view) or the RT dimension (1D view).
		*/
		DRange<3> overall_data_range_;
		
		/// Stores whether or not to show a grid.
		bool show_grid_;
		
		/// The zoom stack.
		std::vector<AreaType> zoom_stack_;
		/// The current position in the zoom stack
		std::vector<AreaType>::iterator zoom_pos_;

		/**
			@brief Updates the diplayed data
			
			The default implementation calls QWidget::update().
			
			This method is reimplemented in the 3D view to update the OpenGL widget.
			
			@param caller_name Name of the calling function (use __PRETTY_FUNCTION__).
		*/
		virtual void update_(const char* caller_name);
		
		///Takes all actions necessary when the modification status of a layer changes (signals etc.)
    void modificationStatus_(Size layer_index, bool modified);
		
		/// Whether to recalculate the data in the buffer when repainting
		bool update_buffer_;

		/// Stores the index of the currently active layer.
		Size current_layer_;

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
		std::vector<DoubleReal> snap_factors_;
		
		/// Rubber band for selected area
		QRubberBand rubber_band_;		
		
    /// External context menu extension
		QMenu* context_add_;
		
    /// Flag that determines if timimg data is printed to the command line
		bool show_timing_;

		/// selected peak
		PeakIndex selected_peak_;
		/// start peak of measuring mode
    PeakIndex measurement_start_;

    /// Data processing setter for peak maps
		template<typename PeakType>
		void addDataProcessing_(MSExperiment<PeakType>& map, DataProcessing::ProcessingAction action) const
		{
			std::set<DataProcessing::ProcessingAction> actions;
			actions.insert(action);
			
			DataProcessing p;
			//actions
			p.setProcessingActions(actions);
			//software
			p.getSoftware().setName("SpectrumCanvas");
			//version
			p.getSoftware().setVersion(VersionInfo::getVersion());
			//time
			p.setCompletionTime(DateTime::now());
			
			for (Size i=0; i<map.size(); ++i)
			{
				map[i].getDataProcessing().push_back(p);          
			}
		}

	};
}

#endif
