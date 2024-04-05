// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Timo Sachsenberg, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// OpenMS
#include <OpenMS/VISUAL/LayerData1DBase.h>
#include <OpenMS/VISUAL/LayerDataChrom.h>
#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/Painter1DBase.h>

// QT
#include <QTextDocument>
#include <QPoint>

// STL
#include <vector>
#include <utility>

// QT
class QAction;

namespace OpenMS
{
  class Annotation1DItem;
  
  /**
   * \brief Manipulates X or Y component of points in the X-Y plane, by assuming one axis (either X or Y axis) has gravity acting upon it.
   *
   * Example: Assume there is a X-Y plane with RT(on X axis) and Intensity(on Y axis). Given a point on the plane, you can make it 'drop down',
   *          to a minimum intensity, by applying a gravity on the Y axis. You can also make the point 'fly up'.
   */
  class Gravitator
  {
  public:

    using AreaXYType = PlotCanvas::GenericArea::AreaXYType;

    /**
     * \brief C'tor to apply gravity on any axis
     * \param axis Which axis
     */
    Gravitator(DIM axis)
    {
      setGravityAxis(axis);
    }

    /**
     * \brief Convenience c'tor, which picks the Intensity dimension from a DimMapper as gravity axis
     * \param unit_mapper Pick the Intensity axis from this mapper (or throw exception). See setGravityAxis().
     */
    Gravitator(const DimMapper<2>& unit_mapper)
    {
      setIntensityAsGravity(unit_mapper);
    }

    /// Which axis is pulling a point downwards (e.g. when plotting sticks)
    /// Note that pulling (see gravitateDown()) will only change the value for the gravity axis.
    /// E.g. with gravity on Y, a Point(X=10, Y=10), will be pulled to Point(X=10, Y=min)
    /// @param axis Either X, or Y
    /// @throws Exception::InvalidValue if @p axis is not X or Y
    void setGravityAxis(DIM axis)
    {
      if (axis != DIM::X && axis != DIM::Y)
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Not a valid axis for 1D plotting", String((int)axis));
      }
      gravity_axis_ = axis;
    }

    /**
     * @brief Convenience function, which picks the Intensity dimension from a DimMapper as gravity axis
     * @param unit_mapper
     * @throw Exception::NotImplemented if @p unit_mapper does not have an Intensity dimension
     */
    void setIntensityAsGravity(const DimMapper<2>& unit_mapper)
    {
      if (unit_mapper.getDim(DIM::X).getUnit() == DIM_UNIT::INT)
      {
        setGravityAxis(DIM::X);
      }
      if (unit_mapper.getDim(DIM::Y).getUnit() == DIM_UNIT::INT)
      {
        setGravityAxis(DIM::Y);
      }
      /// if 1D view has no intensity dimension, go think about what dimension should be gravitational...
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);      
    }

    /// Which axis is affected by gravity?
    DIM getGravityAxis() const
    {
      return gravity_axis_;
    }

    /**
     * \brief Swap gravity axis (from X to Y, or vice versa)
     * \return An orthogonal gravity model
     */
    Gravitator swap() const
    {
      auto r = *this;
      r.setGravityAxis( (r.getGravityAxis() == DIM::X) ? DIM::Y : DIM::X);
      return r;
    }

    /// Pull the point @p p to the current gravity axis, i.e. the lowest point on the Area
    ///
    /// @param p A X-Y data point 
    /// @param area An area which contains the min/max range of X and Y axis
    /// @return A X-Y data point identical to @p p, but with its gravity-axis value changed to the minimum given in @p area
    QPoint gravitateMin(QPoint p, const AreaXYType& area) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() = area.minX();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() = area.minY();
      }
      return p;
    }

    /// Add value of @p delta's gravity dimension to the gravity dimension of point @p p. Other dimensions remain untouched.
    ///
    /// @param p A X-Y data point
    /// @param delta A distance, of which we only use the gravity dimension's part.
    /// @return A X-Y data point identical to @p p, but with its gravity-axis value changed by adding delta.
    QPoint gravitateWith(QPoint p, const QPoint& delta) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() += delta.x();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() += delta.y();
      }
      return p;
    }

    /// Same as gravitateWith()
    template<UInt D>
    DPosition<D> gravitateWith(DPosition<D> p, const DPosition<D>& delta) const
    {
      p[(int)gravity_axis_] += delta[(int)gravity_axis_];
      return p;
    }

    /// Change the value of @p p's gravity dimension to the value of @p targets'. Other dimensions remain untouched.
    ///
    /// @param p A X-Y data point
    /// @param target A target value, of which we only use the gravity dimension's part.
    /// @return A X-Y data point identical to @p p, but with its gravity-axis value changed by target's value.
    QPoint gravitateTo(QPoint p, const QPoint& target) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() = target.x();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() = target.y();
      }
      return p;
    }
    
    /// Same as gravitateTo()
    template<UInt D>
    DPosition<D> gravitateTo(DPosition<D> p, const DPosition<D>& target) const
    {
      p[(int)gravity_axis_] = target[(int)gravity_axis_];
      return p;
    }


    /// Opposite of gravitateMin()
    QPoint gravitateMax(QPoint p, const AreaXYType& area) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() = area.maxX();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() = area.maxY();
      }
      return p;
    }

    /// Pull the point @p p to zero (0) on the current gravity axis.
    ///
    /// @param p A X-Y data point
    /// @return A X-Y data point with its gravity axis set to '0'
    QPoint gravitateZero(QPoint p) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() = 0;
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() = 0;
      }
      return p;
    }

    /// Pull the point @p p to zero (0) on the current gravity axis.
    ///
    /// @param p A X-Y data point
    /// @return A X-Y data point with its gravity axis set to '0'
    template<UInt D>
    DPosition<D> gravitateZero(DPosition<D> p) const
    {
      p[(int)gravity_axis_] = 0;
      return p;
    }

    /// Pull the point @p p to NAN on the current gravity axis.
    ///
    /// @param p A X-Y data point
    /// @return A X-Y data point with its gravity axis set to NAN
    template<UInt D>
    DPosition<D> gravitateNAN(DPosition<D> p) const
    {
      p[(int)gravity_axis_] = std::numeric_limits<float>::quiet_NaN();
      return p;
    }

    /// Get the value of the gravity dimension
    ///
    /// @param p A X-Y data point
    /// @return Either the X or Y component, depending on gravity
    int gravityValue(const QPoint& p) const
    {
      if (gravity_axis_ == DIM::X)
      {
        return p.x();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        return p.y();
      }
      // never reached, but make compilers happy
      return 0;
    }
    /// Get the value of the gravity dimension
    ///
    /// @param p A X-Y data point
    /// @return Either the X or Y component, depending on gravity
    template<UInt D>
    int gravityValue(const DPosition<D>& p) const
    {
      return p[(int)gravity_axis_];
    }

    /// Get the difference of values in the gravity dimension
    ///
    /// @param start The start point in XY coordinates
    /// @param end The end point in XY coordinates
    /// @return The difference of (end-start) in the X or Y component, depending on gravity
    template<UInt D>
    auto gravityDiff(const DPosition<D>& start, const DPosition<D>& end) const
    {
      return end[(int)gravity_axis_] - start[(int)gravity_axis_];
    }

  private:
    /// Where are points in the X-Y plane projected onto when drawing lines?
    DIM gravity_axis_;
  };

  /**
      @brief Canvas for visualization of one or several spectra.

      @image html Plot1DCanvas.png

      The example image shows %Plot1DCanvas displaying a raw data layer and a peak data layer.

      @htmlinclude OpenMS_Plot1DCanvas.parameters

      @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot1DCanvas :
    public PlotCanvas
  {
    Q_OBJECT

public:
    /// Label modes (percentage or absolute) of x axis and y axis
    enum LabelMode
    {
      LM_XABSOLUTE_YABSOLUTE,
      LM_XPERCENT_YABSOLUTE,
      LM_XABSOLUTE_YPERCENT,
      LM_XPERCENT_YPERCENT
    };

    /// extra empty margin added on top to ensure annotations and 100% y-axis label are properly drawn
    constexpr static double TOP_MARGIN{1.09};

    /// Default constructor
    Plot1DCanvas(const Param& preferences, const DIM gravity_axis = DIM::Y, QWidget* parent = nullptr);
    /// Destructor
    ~Plot1DCanvas() override;
    
    /// returns the layer data of the layer @p index
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    const LayerData1DBase& getLayer(Size index) const;
    /// returns the layer data of the layer @p index
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    LayerData1DBase& getLayer(Size index);

    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase 
    const LayerData1DBase& getCurrentLayer() const;
    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    LayerData1DBase& getCurrentLayer();

    /// Get the dimension on which gravity is currently acting upon (usually it's the Y axis' unit)
    const DimBase& getGravityDim() const;

    /// Get the dimension on which gravity is currently not acting upon (the orthogonal axis; usually it's the X axis' unit)
    const DimBase& getNonGravityDim() const;

    /// add a chromatogram layer
    /// @param chrom_exp_sptr An MSExperiment with chromatograms
    /// @param ondisc_sptr OnDisk experiment, as fallback to read the chromatogram from, should @p chrom_exp_sptr.getChromatograms(index) be empty
    /// @param chrom_annotation If OSWData was loaded, pass the shared_pointer from the LayerData. Otherwise leave empty.
    /// @param index Index of the chromatogram to show
    /// @param filename For file change watcher (can be empty, if need be)
    /// @param basename Name of layer (usually the basename of the file)
    /// @param basename_extra Optional suffix of the layer name (e.g. a peptide sequence, or an index '[39]).
    /// @return true on success, false if data was missing etc
    /// @note: this does NOT trigger layerActivated signal for efficiency-reasons. Do it manually afterwards!
    bool addChromLayer(ExperimentSharedPtrType chrom_exp_sptr,
                       ODExperimentSharedPtrType ondisc_sptr, 
                       OSWDataSharedPtrType chrom_annotation,
                       const int index,
                       const String& filename, 
                       const String& basename,
                       const String& basename_extra);

    
    ///Enumerate all available paint styles
    enum DrawModes
    {
      DM_PEAKS,         ///< draw data as peak
      DM_CONNECTEDLINES ///< draw as connected lines
    };

    /// Returns the draw mode of the current layer
    DrawModes getDrawMode() const;

    /// Sets draw mode of the current layer
    void setDrawMode(DrawModes mode);

    // Docu in base class
    void showCurrentLayerPreferences() override;

    /// Returns whether flipped layers exist or not
    bool flippedLayersExist();

    /// Flips the layer with @p index up/downwards
    void flipLayer(Size index);

    /// Returns whether this widget is currently in mirror mode
    bool mirrorModeActive() const;

    /// Sets whether this widget is currently in mirror mode
    void setMirrorModeActive(bool b);

    /// For convenience - calls dataToWidget
    void dataToWidget(const DPosition<2>& peak, QPoint& point, bool flipped = false);
    /// For convenience - calls dataToWidget
    void dataToWidget(const DPosition<2>& xy_point, DPosition<2>& point, bool flipped);

    /// Calls PlotCanvas::dataToWidget_(), takes mirror mode into account
    void dataToWidget(double x, double y, QPoint& point, bool flipped = false);

    /// For convenience - calls widgetToData
    PointXYType widgetToData(const QPoint& pos);

    /// Calls PlotCanvas::widgetToData_(), takes mirror mode into account
    PointXYType widgetToData(double x, double y);

    /**
       @brief converts a distance in axis values to pixel values
    */
    inline void dataToWidgetDistance(double x, double y, QPoint& point)
    {
      dataToWidget_(x, y, point);
      // subtract the 'offset'
      QPoint zero;
      dataToWidget_(0, 0, zero);
      point -= zero;
    }

    /**
      @brief compute distance in data coordinates (unit axis as shown) when moving @p x/y pixel in chart/widget coordinates
    */
    inline PointXYType widgetToDataDistance(double x, double y)
    {
      PointXYType point = Plot1DCanvas::widgetToData(x, y); // call the 1D version, otherwise intensity&mirror modes will not be honored
      // subtract the 'offset'
      PointXYType zero = Plot1DCanvas::widgetToData(0, 0); // call the 1D version, otherwise intensity&mirror modes will not be honored
      point -= zero;
      return point;
    }

    /**
     * \brief Pushes a data point back into the valid data range of the current layer area. Useful for annotation items which were mouse-dragged outside the range by the user.
     * \tparam T A data point, e.g. Peak1D, which may be outside the data area
     * \param data_point
     * \param layer_index The layer of the above data_point (to obtain the data range of the layer)
     */
    template <class T>
    void pushIntoDataRange(T& data_point, const int layer_index)
    { // note: if this is needed for anything other than the 1D Canvas, you need to make sure to call the correct widgetToData/ etc functions --- they work a bit different, depending on Canvas
      auto xy_unit = unit_mapper_.map(data_point); // datatype to xy
      pushIntoDataRange(xy_unit, layer_index);
      unit_mapper_.fromXY(xy_unit, data_point);    // xy to datatype
    }

    /**
     * \brief Pushes a data point back into the valid data range of the current layer area. Useful for annotation items which were mouse-dragged outside the range by the user.
     * \param xy_unit A pair (X and Y coordinate) with values in the units currently used on the axis
     * \param layer_index The layer of the above data_point (to obtain the data range of the layer)
     */
    //template<>  // specialization does not compile when declared within the class on GCC -- even though it should; but I'm not moving it outside! :)
    void pushIntoDataRange(PointXYType& xy_unit, const int layer_index)
    { // note: if this is needed for anything other than the 1D Canvas, you need to make sure to call the correct widgetToData/ etc functions --- they work a bit different, depending on Canvas
      auto p_range = unit_mapper_.fromXY(xy_unit);
      const auto all_range = getLayer(layer_index).getRange();
      p_range.pushInto(all_range);
      xy_unit = unit_mapper_.mapRange(p_range).minPosition();
    }
    
    /// Display a static text box on the top right
    void setTextBox(const QString& html);

    /// ----- Annotations

    /// Add an annotation item for the given peak
    Annotation1DItem* addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color);

    /// ----- Alignment
    /// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
    void performAlignment(Size layer_index_1, Size layer_index_2, const Param& param);

    /// Resets alignment_
    void resetAlignment();

    /// Returns the number of aligned pairs of peaks
    Size getAlignmentSize();

    /// Returns the score of the alignment
    double getAlignmentScore() const;

    /// Returns aligned_peaks_indices_
    std::vector<std::pair<Size, Size> > getAlignedPeaksIndices();

    /// Sets current spectrum index of current layer to @p index
    void activateSpectrum(Size index, bool repaint = true);

    /// Set's the Qt PenStyle of the active layer
    void setCurrentLayerPeakPenStyle(Qt::PenStyle ps);

    /// Actual painting takes place here
    void paint(QPainter* paint_device, QPaintEvent* e);

    /// interesting (e.g., high-intensity) get live annotated with m/s's
    void setDrawInterestingMZs(bool enable);

    /// Return true if interesting m/s are annotated
    bool isDrawInterestingMZs() const;

    // Show/hide ion ladder on top right corner (Identification view)
    void setIonLadderVisible(bool show);

    // Returns true if ion ladder is visible
    bool isIonLadderVisible() const;

    /**
     * \brief Get gravity manipulation object to apply gravity to points
     * \return Gravitator
     */
    const Gravitator& getGravitator() const
    {
      return gr_;
    }

signals:
    /// Requests to display all spectra in 2D plot
    void showCurrentPeaksAs2D();

    /// Requests to display all spectra in 3D plot
    void showCurrentPeaksAs3D();

    /// Requests to display this spectrum (=frame) in ion mobility plot
    void showCurrentPeaksAsIonMobility(const MSSpectrum& spec);

    /// Requests to display all spectra as DIA
    void showCurrentPeaksAsDIA(const Precursor& pc, const MSExperiment& exp);

public slots:
    // Docu in base class
    void activateLayer(Size layer_index) override;
    // Docu in base class
    void removeLayer(Size layer_index) override;
    // Docu in base class
    void updateLayer(Size i) override;
    // Docu in base class
    void horizontalScrollBarChange(int value) override;

protected slots:

    /// Reacts on changed layer parameters
    void currentLayerParamtersChanged_();

protected:

    
    /**
      @brief Convert chart to widget coordinates

      Translates chart (unit) coordinates to widget (pixel) coordinates.
      @param x the chart coordinate x
      @param y the chart coordinate y
      @param point returned widget coordinates
    */
    void dataToWidget_(double x, double y, QPoint& point)
    {
      const auto& xy = visible_area_.getAreaXY();
      const auto h_px = height();
      const auto w_px = width();

      point.setX(int((x - xy.minX()) / xy.width() * w_px));

      if (intensity_mode_ != PlotCanvas::IM_LOG)
      {
        point.setY(int((xy.maxY() - y) / xy.height() * h_px));
      }
      else // IM_LOG
      {
        point.setY(h_px - int(std::log10((y - xy.minY()) + 1) / std::log10(xy.height() + 1) * h_px));
      }
    }

    void dataToWidget_(const DPosition<2>& xy, QPoint& point)
    {
      dataToWidget_(xy.getX(), xy.getY(), point);
    }

    QPoint dataToWidget_(const DPosition<2>& xy)
    {
      QPoint point;
      dataToWidget_(xy.getX(), xy.getY(), point);
      return point;
    }

    // Docu in base class
    bool finishAdding_() override;

    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawCoordinates_(QPainter& painter, const PeakIndex& peak);
    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end);
    
    /// Draws the alignment on @p painter
    void drawAlignment_(QPainter& painter);

    /// internal method, called before calling parent function PlotCanvas::changeVisibleArea_
    void changeVisibleAreaCommon_(const UnitRange& new_area, bool repaint, bool add_to_stack);

    // Docu in base class
    void changeVisibleArea_(VisibleArea new_area, bool repaint = true, bool add_to_stack = false) override;

    /**
        @brief Changes visible area interval

        This method is for convenience only. It calls changeVisibleArea_(const VisibleArea&, bool, bool) .
    */
    void changeVisibleArea_(const AreaXYType& new_area, bool repaint = true, bool add_to_stack = false);

    /**
        @brief Changes visible area interval

        This method is for convenience only. It calls changeVisibleArea_(const VisibleArea&, bool, bool) .
    */
    void changeVisibleArea_(const UnitRange& new_area, bool repaint = true, bool add_to_stack = false);


    /// Draws a highlighted peak; if draw_elongation is true, the elongation line is drawn (for measuring)
    void drawHighlightedPeak_(Size layer_index, const PeakIndex& peak, QPainter& painter, bool draw_elongation = false);


    /**
        @brief Zooms fully out and resets the zoom stack

        Sets the visible area to the initial value, such that all data (for the current spec/chrom/...) is shown.

        @param repaint If @em true a repaint is forced. Otherwise only the new area is set.
    */
    void resetZoom(bool repaint = true) override
    {
      zoomClear_();
      PlotCanvas::changeVisibleArea_(visible_area_.cloneWith(overall_data_range_1d_), repaint, true);
    }
        
    /// Recalculates the current scale factor based on the specified layer (= 1.0 if intensity mode != IM_PERCENTAGE)
    void recalculatePercentageFactor_(Size layer_index);

    /**
        @brief Recalculates the overall_data_range_ (by calling PlotCanvas::recalculateRanges_)
               plus the overall_data_range_1d_ (which only takes into account the current spec/chrom/.. of all layers)

        A small margin is added to each side of the range in order to display all data.
    */
    void recalculateRanges_() override
    {
      PlotCanvas::recalculateRanges_(); // for: overall_data_range_
      // the same thing for: overall_data_range_1d_
      RangeType& layer_range_1d = overall_data_range_1d_;
      layer_range_1d.clearRanges();

      for (Size layer_index = 0; layer_index < getLayerCount(); ++layer_index)
      {
        layer_range_1d.extend(getLayer(layer_index).getRange1D());
      }
      // add 4% margin (2% left, 2% right) to all dimensions, except the current gravity axes's minimum (usually intensity)
      layer_range_1d.scaleBy(1.04);

      // set minimum intensity to 0 (avoid negative intensities and show full height of peaks in case their common minimum is large)
      auto& gravity_range = getGravityDim().map(layer_range_1d);
      gravity_range.setMin(0);
      
      // make sure that each dimension is not a single point (axis widget won't like that)
      // (this needs to be the last command to ensure this property holds when leaving the function!)
      layer_range_1d.minSpanIfSingular(1);
    }

    // Docu in base class
    void updateScrollbars_() override;
    // Docu in base class
    void intensityModeChange_() override;
    
    
    /// Adjust the gravity axis (usually y-axis with intensity) according to the given range on the x-axis 
    /// (since the user cannot freely choose the limits of this axis in 1D View)
    RangeAllType correctGravityAxisOfVisibleArea_(UnitRange area);
    
    /** @name Reimplemented QT events */
    //@{
    void paintEvent(QPaintEvent* e) override;
    void mousePressEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;
    //@}

    // docu in base class
    void zoomForward_() override;
    // docu in base class
    void zoom_(int x, int y, bool zoom_in) override;
    // docu in base class
    void translateLeft_(Qt::KeyboardModifiers m) override;
    // docu in base class
    void translateRight_(Qt::KeyboardModifiers m) override;
    // docu in base class
    void translateForward_() override;
    // docu in base class
    void translateBackward_() override;

    // docu in base class
    void paintGridLines_(QPainter& painter) override;

    /// Find peak next to the given position
    PeakIndex findPeakAtPosition_(QPoint);

    /// Shows dialog and calls addLabelAnnotation_
    void addUserLabelAnnotation_(const QPoint& screen_position);
    /// Adds an annotation item at the given screen position
    void addLabelAnnotation_(const QPoint& screen_position, const QString& label_text);
    /// Shows dialog and calls addPeakAnnotation_
    void addUserPeakAnnotation_(PeakIndex near_peak);

    /// Ensure that all annotations are within data range
    void ensureAnnotationsWithinDataRange_();

    friend class Painter1DChrom;
    friend class Painter1DPeak;
    friend class Painter1DIonMobility;

    /////////////////////
    ////// data members
    /////////////////////

    /// The data range (m/z, RT and intensity) of the current(!) spec/chrom for all layers
    RangeType overall_data_range_1d_;

    /// Draw modes (for each layer) - sticks or connected lines
    std::vector<DrawModes> draw_modes_;
    /// Draw style (for each layer)
    std::vector<Qt::PenStyle> peak_penstyle_;

    /// start point of "ruler" in pixel coordinates for measure mode
    QPoint measurement_start_point_px_;
    /// Indicates whether this widget is currently in mirror mode
    bool mirror_mode_ = false;
    /// Indicates whether annotation items are just being moved on the canvas
    bool moving_annotations_ = false;
    /// Indicates whether an alignment is currently visualized
    bool show_alignment_ = false;
    /// Layer index of the first alignment layer
    Size alignment_layer_1_;
    /// Layer index of the second alignment layer
    Size alignment_layer_2_;
    /// Stores the alignment as MZ values of pairs of aligned peaks in both spectra
    std::vector<std::pair<double, double> > aligned_peaks_mz_delta_;
    /// Stores the peak indices of pairs of aligned peaks in both spectra
    std::vector<std::pair<Size, Size> > aligned_peaks_indices_;
    /// Stores the score of the last alignment
    double alignment_score_ = 0.0;
    /// whether the ion ladder is displayed on the top right corner in ID view
    bool ion_ladder_visible_ = true;
    /// annotate interesting peaks with m/z's
    bool draw_interesting_MZs_ = false;
    /// The text box in the upper left corner with the current data coordinates of the cursor
    QTextDocument text_box_content_;
    /// handles pulling/pushing of points to the edges of the widget
    Gravitator gr_;
  };




} // namespace OpenMS

