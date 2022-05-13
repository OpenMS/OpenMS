// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
    /// E.g. with gravity on Y, an Point(X=10, Y=10), will be pulled to Point(X=10, Y=min)
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
    template<int D>
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
    template<int D>
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
    template<int D>
    DPosition<D> gravitateZero(DPosition<D> p) const
    {
      p[(int)gravity_axis_] = 0;
      return p;
    }

    /// Pull the point @p p to NAN on the current gravity axis.
    ///
    /// @param p A X-Y data point
    /// @return A X-Y data point with its gravity axis set to NAN
    QPoint gravitateNAN(QPoint p) const
    {
      if (gravity_axis_ == DIM::X)
      {
        p.rx() = std::numeric_limits<float>::quiet_NaN();
      }
      else if (gravity_axis_ == DIM::Y)
      {
        p.ry() = std::numeric_limits<float>::quiet_NaN();
      }
      return p;
    }

    /// Pull the point @p p to NAN on the current gravity axis.
    ///
    /// @param p A X-Y data point
    /// @return A X-Y data point with its gravity axis set to NAN
    template<int D>
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

    /// Get the difference of values in the gravity dimension
    ///
    /// @param start The start point in XY coordinates
    /// @param end The end point in XY coordinates
    /// @return The difference of (end-start) in the X or Y component, depending on gravity
    template<int D>
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
    
    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    const LayerData1DBase& getLayer(Size index) const override
    {
      return dynamic_cast<const LayerData1DBase&>(layers_.getLayer(index));
    }
    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    LayerData1DBase& getCurrentLayer(Size index)
    {
      return dynamic_cast<LayerData1DBase&>(layers_.getLayer(index));
    }

    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase 
    const LayerData1DBase& getCurrentLayer() const override
    {
      return dynamic_cast<const LayerData1DBase&>(layers_.getCurrentLayer());
    }
    /// returns the layer data of the active layer
    /// @throws std::bad_cast exception if the current layer is not a LayerData1DBase
    LayerData1DBase& getCurrentLayer()
    {
      return dynamic_cast<LayerData1DBase&>(layers_.getCurrentLayer());
    }

    const DimBase& getGravityDim() const
    {
      return unit_mapper_.getDim(getGravitator().getGravityAxis());
    }

    const DimBase& getNonGravityDim() const
    {
      return unit_mapper_.getDim(getGravitator().swap().getGravityAxis());
    }

    /// add a chromatogram layer
    /// @param chrom_exp_sptr An MSExperiment with chromatograms
    /// @param ondisc_sptr OnDisk experiment, as fallback to read the chromatogram from, should @p chrom_exp_sptr.getChromatograms(index) be empty
    /// @param OSWDataSharedPtrType If OSWData was loaded, pass the shared_pointer from the LayerData. Otherwise leave empty.
    /// @param index Index of the chromatogram to show
    /// @param filename For file change watcher (can be empty, if need be)
    /// @param caption Name of layer
    /// @param multiple_select .... not sure ...
    /// @return true on success, false if data was missing etc
    /// @note: this does NOT trigger layerActivated signal for efficiency-reasons. Do it manually afterwards!
    bool addChromLayer(ExperimentSharedPtrType chrom_exp_sptr,
                       ODExperimentSharedPtrType ondisc_sptr, 
                       OSWDataSharedPtrType chrom_annotation,
                       const int index,
                       const String& filename, 
                       const String& caption, 
                       const bool multiple_select);

    
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

    // Docu in base class
    void saveCurrentLayer(bool visible) override;

    /// Returns whether flipped layers exist or not
    bool flippedLayersExist();

    /// Flips the layer with @p index up/downwards
    void flipLayer(Size index);

    /// Returns whether this widget is currently in mirror mode
    bool mirrorModeActive() const;

    /// Sets whether this widget is currently in mirror mode
    void setMirrorModeActive(bool b);

    /// For convenience - calls dataToWidget
    void dataToWidget(const DPosition<2>& peak, QPoint& point, bool flipped = false, bool percentage = true);
    /// For convenience - calls dataToWidget
    void dataToWidget(const DPosition<2>& xy_point, DPosition<2>& point, bool flipped, bool percentage);

    /// Calls PlotCanvas::dataToWidget_(), takes mirror mode into account
    void dataToWidget(double x, double y, QPoint& point, bool flipped = false, bool percentage = false);

    /// For convenience - calls widgetToData
    PointXYType widgetToData(const QPoint& pos, bool percentage = false);

    /// Calls PlotCanvas::widgetToData_(), takes mirror mode into account
    PointXYType widgetToData(double x, double y, bool percentage = false);

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

    template<>
    void pushIntoDataRange<PointXYType>(PointXYType& xy_unit, const int layer_index)
    { // note: if this is needed for anything other than the 1D Canvas, you need to make sure to call the correct widgetToData/ etc functions --- they work a bit different, depending on Canvas
      auto p_range = unit_mapper_.fromXY(xy_unit);
      const auto& all_range = getLayer(layer_index).getRange();
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

    /// Requests to display all spectra in ion mobility plot
    void showCurrentPeaksAsIonMobility(const MSSpectrum& spec);

    /// Requests to display all spectra as DIA
    void showCurrentPeaksAsDIA();

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
    // Docu in base class
    bool finishAdding_() override;

    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawCoordinates_(QPainter& painter, const PeakIndex& peak);
    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end);
    
    /// Draws the alignment on @p painter
    void drawAlignment_(QPainter& painter);

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

    /// Recalculates the current scale factor based on the specified layer (= 1.0 if intensity mode != IM_PERCENTAGE)
    void updatePercentageFactor_(Size layer_index);

    // Docu in base class
    void recalculateSnapFactor_() override;
    // Docu in base class
    void updateScrollbars_() override;
    // Docu in base class
    void intensityModeChange_() override;

    
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

    friend class Painter1DPeak;

    /////////////////////
    ////// data members
    /////////////////////

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

