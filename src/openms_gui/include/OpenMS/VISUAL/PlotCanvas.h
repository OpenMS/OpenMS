// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

//OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/DATASTRUCTURES/DRange.h>
#include <OpenMS/KERNEL/DimMapper.h>
#include <OpenMS/VISUAL/LayerDataBase.h>
#include <OpenMS/VISUAL/MISC/CommonDefs.h>

//QT
#include <QtWidgets>
#include <QRubberBand>

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
  class PlotWidget;
  class LayerDataChrom;
  class LayerDataPeak;
  class LayerDataFeature;
  class LayerDataConsensus;

  using LayerDataBaseUPtr = std::unique_ptr<LayerDataBase>;
  using LayerDataChromUPtr = std::unique_ptr<LayerDataChrom>;
  using LayerDataPeakUPtr = std::unique_ptr<LayerDataPeak>;
  using LayerDataFeatureUPtr = std::unique_ptr<LayerDataFeature>;
  using LayerDataConsensusUPtr = std::unique_ptr<LayerDataConsensus>;

  /**
    A class to manage a stack of layers as shown in the layer widget in TOPPView.
    The order of layers is automatically determined based on LayerDataBase::type (in short: peak data below, ID data on top).

  */
  class LayerStack
  {
    public:
      /// adds a new layer and makes it the current layer
      /// @param new_layer Takes ownership of the layer!
      void addLayer(LayerDataBaseUPtr new_layer);

      const LayerDataBase& getLayer(const Size index) const;

      LayerDataBase& getLayer(const Size index);

      const LayerDataBase& getCurrentLayer() const;

      LayerDataBase& getCurrentLayer();

      /// throws Exception::IndexOverflow unless @p index is smaller than getLayerCount()
      void setCurrentLayer(Size index);
      
      Size getCurrentLayerIndex() const;
      
      bool empty() const;

      Size getLayerCount() const;

      void removeLayer(Size layer_index);

      void removeCurrentLayer();
  
  protected:
      std::vector<LayerDataBaseUPtr> layers_;
  private:
      Size current_layer_ = -1;
  };

  /**
      @brief Base class for visualization canvas classes

      This class is the base class for the spectrum data views which are used
      for 1D, 2D and 3D visualization of data. In TOPPView, each PlotCanvas
      is paired with an enclosing PlotWidget (see also the
      getPlotWidget() function that provides a back-reference).  To provide
      additional spectrum views, you can derive from this class and you should
      also create a subclass from PlotWidget which encloses your class
      derived from PlotCanvas. A spectrum canvas can display multiple data
      layers at the same time (see layers_ member variable).

      The actual data to be displayed is stored as a vector of LayerDataBase
      objects which hold the actual data.  It also stores information about the
      commonly used constants such as ActionModes or IntensityModes.

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

      @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI PlotCanvas : public QWidget, public DefaultParamHandler
  {
    Q_OBJECT

  public:
    /**@name Type definitions */
    //@{

    /// Main data type (experiment)
    typedef LayerDataBase::ExperimentType ExperimentType;
    /// Main managed data type (experiment)
    typedef LayerDataBase::ExperimentSharedPtrType ExperimentSharedPtrType;
    typedef LayerDataBase::ConstExperimentSharedPtrType ConstExperimentSharedPtrType;
    typedef LayerDataBase::ODExperimentSharedPtrType ODExperimentSharedPtrType;
    /// Main data type (features)
    typedef LayerDataBase::FeatureMapType FeatureMapType;
    /// Main managed data type (features)
    typedef LayerDataBase::FeatureMapSharedPtrType FeatureMapSharedPtrType;
    /// Main data type (consensus features)
    typedef LayerDataBase::ConsensusMapType ConsensusMapType;
    /// Main managed data type (consensus features)
    typedef LayerDataBase::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

    /// Spectrum type
    typedef ExperimentType::SpectrumType SpectrumType;
    /// Spectrum iterator type (iterates over peaks)
    typedef SpectrumType::ConstIterator SpectrumConstIteratorType;
    /// Peak type
    typedef SpectrumType::PeakType PeakType;
    /// a generic range for the most common units
    using RangeType = RangeAllType;

    /// The range of data shown on the X and Y axis (unit depends on runtime config)
    using AreaXYType = Area<2>::AreaXYType;
    /// The visible range of data on X and Y axis as shown on plot axis (not necessarily the range of actual data, e.g. no data to show).
    using VisibleArea = Area<2>;

    /// A generic range of data on X and Y axis as shown on plot axis
    using GenericArea = Area<2>;

    /// The number of pixels on the axis. The lower point of the area will be zero. The maxima will reflect the number of pixels
    /// in either dimension. Note that the unit of the Axis etc is in PIXELS (-> *not* in seconds, m/z or whatever)
    using PixelArea = Area<2>;

    using UnitRange = RangeAllType;

    using PointOnAxis = DimMapper<2>::Point;

    /// Mouse action modes
    enum ActionModes
    {
      AM_TRANSLATE, ///< translate
      AM_ZOOM,      ///< zoom
      AM_MEASURE    ///< measure
    };

    /// Display modes of intensity
    enum IntensityModes
    {
      IM_NONE,       ///< Normal mode: f(x)=x; y-axis shows the maximum of all layers, no scaling
      IM_PERCENTAGE, ///< Shows intensities normalized by each layer's maximum: f(x)=x/max(x)*100
      IM_SNAP,       ///< Shows the maximum displayed intensity (across all layers) as if it was the overall maximum intensity
      IM_LOG         ///< Logarithmic version of normal mode
    };

    //@}

    /// Default constructor
    PlotCanvas(const Param& preferences, QWidget* parent = nullptr);

    /// Destructor
    ~PlotCanvas() override;

    /**
        @brief Sets the spectrum widget.

        Sets the enclosing spectrum widget. Call this from your
        PlotWidget derived class.
        @param widget the spectrum widget
    */
    inline void setPlotWidget(PlotWidget* widget)
    {
      spectrum_widget_ = widget;
    }

    /**
        @brief Returns the spectrum widget.

        Returns the enclosing spectrum widget
        @return the spectrum widget
    */
    inline PlotWidget* getPlotWidget() const
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
    inline IntensityModes getIntensityMode() const
    {
      return intensity_mode_;
    }

    /**
        @brief Sets the intensity mode

        Sets the intensity mode and calls intensityModeChange_()

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
    const LayerDataBase& getLayer(Size index) const
    {
      return layers_.getLayer(index);
    }
    /// returns the layer data with index @p index
    LayerDataBase& getLayer(Size index)
    {
      return layers_.getLayer(index);
    }

    /// returns the layer data of the active layer
    const LayerDataBase& getCurrentLayer() const
    {
      return layers_.getCurrentLayer();
    }
    /// returns the layer data of the active layer
    LayerDataBase& getCurrentLayer()
    {
      return layers_.getCurrentLayer();
    }

    /// returns the index of the active layer
    inline Size getCurrentLayerIndex() const
    {
      return layers_.getCurrentLayerIndex();
    }

    /// returns a layer flag of the current layer
    bool getLayerFlag(LayerDataBase::Flags f) const
    {
      return getLayerFlag(layers_.getCurrentLayerIndex(), f);
    }

    /// sets a layer flag of the current layer
    void setLayerFlag(LayerDataBase::Flags f, bool value)
    {
      setLayerFlag(layers_.getCurrentLayerIndex(), f, value);
    }

    /// returns a layer flag of the layer @p layer
    bool getLayerFlag(Size layer, LayerDataBase::Flags f) const
    {
      return layers_.getLayer(layer).flags.test(f);
    }

    /// sets a layer flag of the layer @p layer
    void setLayerFlag(Size layer, LayerDataBase::Flags f, bool value)
    {
      // abort if there are no layers
      if (layers_.empty())
        return;

      layers_.getLayer(layer).flags.set(f, value);
      update_buffer_ = true;
      update();
    }

    inline void setLabel(LayerDataBase::LabelType label)
    {
      // abort if there are no layers
      if (layers_.empty())
        return;
      layers_.getCurrentLayer().label = label;
      update_buffer_ = true;
      update();
    }

    /**
        @brief Returns the currently visible area. This is the authority which determines the X and Y axis' scale.

        @see visible_area_
    */
    const VisibleArea& getVisibleArea() const
    {
      return visible_area_;
    }

    /// Given a 2D axis coordinate, is it in the currently visible area? (useful to avoid plotting stuff outside the visible area)
    /// Note: The input @p p must have unit coordinates (i.e. the result of widgetToData_), not pixel coordinates.
    bool isVisible(const PointOnAxis& p) const
    {
      return visible_area_.getAreaXY().encloses(p);
    }

    /// Get the number of pixels of the current canvas (this is independent of the current visible area and zoom level).
    /// It's just the size of the canvas.
    PixelArea getPixelRange() const
    {
      int X_pixel_count = buffer_.width();
      int Y_pixel_count = buffer_.height();
      PixelArea area(&unit_mapper_);
      area.setArea(AreaXYType(0, 0, X_pixel_count, Y_pixel_count));
      return area;
    }

    /**
        @brief Sets the filters applied to the data before drawing (for the current layer)
    */
    virtual void setFilters(const DataFilters& filters);

    /**
        @name Dataset handling methods

        @see changeVisibility
    */
    //@{
    /// Returns the number of layers
    inline Size getLayerCount() const
    {
      return layers_.getLayerCount();
    }

    /// change the active layer (the one that is used for selecting and so on)
    virtual void activateLayer(Size layer_index) = 0;
    /// removes the layer with index @p layer_index
    virtual void removeLayer(Size layer_index) = 0;

    /// removes all layers by calling removeLayer() for all layer indices (from highest to lowest)
    void removeLayers()
    {
      for (Size i = getLayerCount(); i > 0; --i)
      {
        removeLayer(i - 1);
      }
      visible_area_.clear(); // reset visible area
    }

    /// Add an already constructed layer (e.g. for projections)
    bool addLayer(std::unique_ptr<LayerData1DBase> layer);

    /**
      @brief Add a peak data layer

      @param map Shared pointer to input map. It can be performed in constant time and does not double the required memory.
      @param od_map Shared pointer to on disk data which potentially caches some data to save memory (the map can be empty, but do not pass nullptr).
      @param filename This @em absolute filename is used to monitor changes in the file and reload the data
      @param caption The caption of the layer (shown in the layer window)
      @param use_noise_cutoff Add a noise filter which removes low-intensity peaks

      @return If a new layer was created
    */
    bool addPeakLayer(const ExperimentSharedPtrType& map,
                      ODExperimentSharedPtrType od_map,
                      const String& filename = "",
                      const String& caption = "",
                      const bool use_noise_cutoff = false);

    /**
      @brief Add a chrom data layer

      @param map Shared pointer to input map. It can be performed in constant time and does not double the required memory.
      @param od_map Shared pointer to on disk data which potentially caches some data to save memory (the map can be empty, but do not pass nullptr).
      @param filename This @em absolute filename is used to monitor changes in the file and reload the data
      @param caption The caption of the layer (shown in the layer window)

      @return If a new layer was created
    */
    bool addChromLayer(const ExperimentSharedPtrType& map, ODExperimentSharedPtrType od_map, const String& filename = "", const String& caption = "");


    /**
        @brief Add a feature data layer

        @param map Shared Pointer to input map. It can be performed in constant time and does not double the required memory.
        @param filename This @em absolute filename is used to monitor changes in the file and reload the data
        @param caption The caption of the layer (shown in the layer window)

        @return If a new layer was created
    */
    bool addLayer(FeatureMapSharedPtrType map, const String& filename = "", const String& caption = "");

    /**
        @brief Add a consensus feature data layer

        @param map Shared Pointer to input map. It can be performed in constant time and does not double the required memory.
        @param filename This @em absolute filename is used to monitor changes in the file and reload the data
        @param caption The caption of the layer (shown in the layer window)

        @return If a new layer was created
    */
    bool addLayer(ConsensusMapSharedPtrType map, const String& filename = "", const String& caption = "");
    //@}

    /**
        @brief Add an identification data layer

        @param peptides Input list of peptides, which has to be mutable and will be empty after adding. 
               Swapping is used to insert the data. It can be performed in constant time and does not double
               the required memory.
        @param filename This @em absolute filename is used to monitor changes in the file and reload the data
        @param caption The caption of the layer (shown in the layer window)

        @return If a new layer was created
    */
    bool addLayer(std::vector<PeptideIdentification>& peptides, const String& filename = "", const String& caption = "");

    /// Returns the minimum intensity of the active layer
    inline float getCurrentMinIntensity() const
    {
      return layers_.getCurrentLayer().getMinIntensity();
    }

    /// Returns the maximum intensity of the active layer
    inline float getCurrentMaxIntensity() const
    {
      return layers_.getCurrentLayer().getMaxIntensity();
    }

    /// Returns the minimum intensity of the layer with index @p index
    inline float getMinIntensity(Size index) const
    {
      return getLayer(index).getMinIntensity();
    }

    /// Returns the maximum intensity of the layer with index @p index
    inline float getMaxIntensity(Size index) const
    {
      return getLayer(index).getMaxIntensity();
    }

    /// Sets the @p name of layer @p i
    void setLayerName(Size i, const String& name);

    /// Gets the name of layer @p i
    String getLayerName(Size i);

    /// Sets the parameters of the current layer
    inline void setCurrentLayerParameters(const Param& param)
    {
      getCurrentLayer().param = param;
      emit preferencesChange();
    }

    /**
        @brief Returns the area which encloses all data points of all layers.

        @see overall_data_range_
    */
    virtual const RangeType& getDataRange() const;

    /**
        @brief Returns the first intensity scaling factor for 'snap to maximum intensity mode' (for the currently visible data range).

        @see snap_factors_
    */
    double getSnapFactor();

    /// Returns the percentage factor
    double getPercentageFactor() const;

    /// Shows the preferences dialog of the active layer
    virtual void showCurrentLayerPreferences() = 0;

    /**
        @brief Shows a dialog with the meta data

        @param modifiable indicates if the data can be modified.
        @param index If given, the meta data of the corresponding element (spectrum, feature, consensus feature) is shown instead of the layer meta data.
    */
    virtual void showMetaData(bool modifiable = false, Int index = -1);

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
    virtual void resetZoom(bool repaint = true);

    /**
        @brief Sets the visible area.

        Sets the visible area to a new value and emits visibleAreaChanged() if the area is different from the old one.

        @param area the new visible area
    */
    void setVisibleArea(const VisibleArea& area);

    /**
        @brief Sets the visible area.

        Sets the visible area to a new value and emits visibleAreaChanged() if the area is different from the old one.

        @param area the new visible area
    */
    void setVisibleArea(const RangeAllType& area);

    /**
        @brief Sets the visible area.

        Sets the visible area to a new value and emits visibleAreaChanged() if the area is different from the old one.

        @param area the new visible area
    */
    void setVisibleArea(const AreaXYType& area);

    /**
     * @brief Set only the visible area for the x axis; other axes are untouched.
     * @param min 
     * @param max 
    */
    void setVisibleAreaX(double min, double max);

    /**
     * @brief Set only the visible area for the y axis; other axes are untouched.
     * @param min
     * @param max
     */
    void setVisibleAreaY(double min, double max);

    /**
        @brief Saves the current layer data.

        @param visible If true, only the visible data is stored. Otherwise the whole data is stored.
    */
    void saveCurrentLayer(bool visible);

    /**
        @brief Notifies the canvas that the horizontal scrollbar has been moved.

        Reimplement this slot to react on scrollbar events.
    */
    virtual void horizontalScrollBarChange(int value);

    /**
        @brief Notifies the canvas that the vertical scrollbar has been moved.

        Reimplement this slot to react on scrollbar events.
    */
    virtual void verticalScrollBarChange(int value);

    /// Sets the additional context menu. If not 0, this menu is added to the context menu of the canvas
    void setAdditionalContextMenu(QMenu * menu);

    /// Updates layer @p i when the data in the corresponding file changes
    virtual void updateLayer(Size i) = 0;


    /**
     * \brief Get the Area in pixel coordinates of the current canvas for X and Y axis.
     * \return
     */
    AreaXYType canvasPixelArea() const
    {
      return AreaXYType({0, 0}, {(float)width(), (float)height()});
    }

    /**
     * \brief Get Mapper to translate between values for axis (X/Y) and units (m/z, RT, intensity, ...)
     * \return The translation from axis to units
     */
    const DimMapper<2>& getMapper() const;

    /**
     * \brief Set a new mapper for the canvas.
     * \param mapper The new mapper for translating between units and axis
     */
    void setMapper(const DimMapper<2>& mapper);

  signals:
    /// Signal emitted whenever the modification status of a layer changes (editing and storing)
    void layerModficationChange(Size layer, bool modified);

    /// Signal emitted whenever a new layer is activated within the current window
    void layerActivated(QWidget * w);

    /// Signal emitted whenever the zoom changed
    void layerZoomChanged(QWidget * w);

    /**
        @brief Change of the visible area

        Signal emitted whenever the visible area changes.
        @param area The new visible area.
    */
    void visibleAreaChanged(const VisibleArea& area);

    /// Emitted when the cursor position changes (for displaying e.g. in status bar)
    void sendCursorStatus(const String& x_value, const String& y_value);

    /// Emits a status message that should be displayed for @p time ms. If @p time is 0 the message should be displayed until the next message is emitted.
    void sendStatusMessage(std::string message, OpenMS::UInt time);

    /// Forces recalculation of axis ticks in the connected widget.
    void recalculateAxes();

    /// Triggers the update of the vertical scrollbar
    void updateVScrollbar(float f_min, float disp_min, float disp_max, float f_max);

    /// Triggers the update of the horizontal scrollbar
    void updateHScrollbar(float f_min, float disp_min, float disp_max, float f_max);

    /// Toggle axis legend visibility change
    void changeLegendVisibility();

    /// Emitted when the action mode changes
    void actionModeChange();

    /// Emitted when the layer preferences have changed
    void preferencesChange();

protected slots:

    ///Updates the cursor according to the current action mode
    void updateCursor_();

protected:

    /// Draws several lines of text to the upper right corner of the widget
    void drawText_(QPainter & painter, const QStringList& text);

    /// Returns the m/z value of an identification depending on the m/z source of the layer (precursor mass/theoretical peptide mass)
    double getIdentificationMZ_(const Size layer_index,
                                    const PeptideIdentification & peptide) const;

    /// Method that is called when a new layer has been added
    virtual bool finishAdding_() = 0;

    /// remove already added layer which did not pass final checks in finishAdding_()
    /// @param error_message Optional error message to show as messagebox
    void popIncompleteLayer_(const QString& error_message = "");

    ///@name reimplemented QT events
    //@{
    void resizeEvent(QResizeEvent * e) override;
    void wheelEvent(QWheelEvent * e) override;
    void keyPressEvent(QKeyEvent * e) override;
    void keyReleaseEvent(QKeyEvent * e) override;
    void focusOutEvent(QFocusEvent * e) override;
    void leaveEvent(QEvent * e) override;
    void enterEvent(QEnterEvent * e) override;
    //@}

    /// This method is called whenever the intensity mode changes. Reimplement if you need to react on such changes.
    virtual void intensityModeChange_();

    /// Call this whenever the DimMapper receives new dimensions; will update the axes and scrollbars
    void dimensionsChanged_();

    /**
        @brief Sets the visible area

        Changes the visible area, adjusts the zoom stack and notifies interested clients about the change.
        If the area is outside the overall data range, the new area is pushed back into the overall range.

        @param new_area The new visible area.
        @param repaint If @em true, a complete repaint is forced.
        @param add_to_stack If @em true the new area is to add to the zoom_stack_.
    */
    virtual void changeVisibleArea_(VisibleArea new_area, bool repaint = true, bool add_to_stack = false);

    /**
        @brief Recalculates the intensity scaling factor for 'snap to maximum intensity mode'.

        @see snap_factors_
    */
    virtual void recalculateSnapFactor_();

    ///@name Zoom stack methods
    //@{
    /// Zooms such that screen point x, y would still point to the same data point
    virtual void zoom_(int x, int y, bool zoom_in);
    ///Go backward in zoom history
    void zoomBack_();
    ///Go forward in zoom history
    virtual void zoomForward_();
    /// Add a visible area to the zoom stack
    void zoomAdd_(const VisibleArea& area);
    /// Clears the zoom stack and invalidates the current zoom position. After calling this, a valid zoom position has to be added immediately.
    void zoomClear_();
    //@}

    ///@name Translation methods, which are called when cursor buttons are pressed
    //@{
    /// Translation bound to the 'Left' key
    virtual void translateLeft_(Qt::KeyboardModifiers m);
    /// Translation bound to the 'Right' key
    virtual void translateRight_(Qt::KeyboardModifiers m);
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
        @brief Convert widget (pixel) to chart (unit) coordinates

        Translates widget coordinates to chart coordinates.

        @param x the widget coordinate x
        @param y the widget coordinate y
        @return chart coordinates
    */
    inline PointXYType widgetToData_(double x, double y)
    {
      const auto& xy = visible_area_.getAreaXY();
      return PointXYType(
                xy.minX() + x / width() * xy.width(),
                xy.minY() + (height() - y) / height() * xy.height()
                );
    }

    /// Calls widgetToData_ with x and y position of @p pos
    inline PointXYType widgetToData_(const QPoint& pos)
    {
      return widgetToData_(pos.x(), pos.y());
    }

    /// Helper function to paint grid lines
    virtual void paintGridLines_(QPainter & painter);

    /// Buffer that stores the actual peak information
    QImage buffer_;

    /// Mapper for X and Y axis
    DimMapper<2> unit_mapper_;

    /// Stores the current action mode (Pick, Zoom, Translate)
    ActionModes action_mode_ = AM_TRANSLATE;

    /// Stores the used intensity mode function
    IntensityModes intensity_mode_ = IM_NONE;

    /// Layer data
    LayerStack layers_;

    /**
        @brief Stores the currently visible area in data units (e.g. seconds, m/z, intensity etc) and axis (X,Y) area.

        This is always (and only) the data shown in the widget.
        In 1D, the gravity axis may get some headroom on the y-axis, see Plot1DCanvas::correctGravityAxisOfVisibleArea_()
    */
    VisibleArea visible_area_;

    /**
        @brief Recalculates the overall_data_range_

        A small margin is added to each side of the range in order to display all data.
    */
    virtual void recalculateRanges_();

    /**
        @brief Stores the data range (m/z, RT and intensity) of all layers
    */
    RangeType overall_data_range_;

    /// Stores whether or not to show a grid.
    bool show_grid_ = true;

    /// The zoom stack.
    std::vector<VisibleArea> zoom_stack_;
    /// The current position in the zoom stack
    std::vector<VisibleArea>::iterator zoom_pos_ = zoom_stack_.end();

    /**
        @brief Updates the displayed data

        The default implementation calls QWidget::update().

        This method is reimplemented in the 3D view to update the OpenGL widget.

        @param caller_name Name of the calling function (use OPENMS_PRETTY_FUNCTION).
    */
    virtual void update_(const char* caller_name);

    /// Takes all actions necessary when the modification status of a layer changes (signals etc.)
    void modificationStatus_(Size layer_index, bool modified);

    /// Whether to recalculate the data in the buffer when repainting
    bool update_buffer_ = false;

    /// Back-pointer to the enclosing spectrum widget
    PlotWidget* spectrum_widget_ = nullptr;

    /// start position of mouse actions
    QPoint last_mouse_pos_;

    /**
        @brief Intensity scaling factor for relative scale with multiple layers.

        In this mode all layers are scaled to the same maximum.
        FIXME: this factor changes, depending on the layer which is currently plotted! Ouch!
    */
    double percentage_factor_ = 1.0;

    /**
        @brief Intensity scaling factor for 'snap to maximum intensity mode'.

        In this mode the highest currently visible intensity is treated like the maximum overall intensity.

        Only used in 2D mode.
    */
    std::vector<double> snap_factors_;

    /// Rubber band for selected area
    QRubberBand rubber_band_;

    /// External context menu extension
    QMenu* context_add_ = nullptr;

    /// Flag that determines if timing data is printed to the command line
    bool show_timing_ = false;

    /// selected peak
    PeakIndex selected_peak_;
    /// start peak of measuring mode
    PeakIndex measurement_start_;

    /// Data processing setter for peak maps
    void addDataProcessing_(PeakMap & map, DataProcessing::ProcessingAction action) const
    {
      std::set<DataProcessing::ProcessingAction> actions;
      actions.insert(action);

      DataProcessingPtr p = boost::shared_ptr<DataProcessing>(new DataProcessing);
      //actions
      p->setProcessingActions(actions);
      //software
      p->getSoftware().setName("PlotCanvas");
      //version
      p->getSoftware().setVersion(VersionInfo::getVersion());
      //time
      p->setCompletionTime(DateTime::now());

      for (Size i = 0; i < map.size(); ++i)
      {
        map[i].getDataProcessing().push_back(p);
      }
    }

  };
}
