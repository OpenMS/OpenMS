// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#ifndef OPENMS_VISUAL_SPECTRUM1DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM1DCANVAS_H

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>

#include <QTextDocument>

// STL
#include <vector>
#include <utility>

// QT
class QAction;

namespace OpenMS
{
  /**
      @brief Canvas for visualization of one or several spectra.

      @image html Spectrum1DCanvas.png

      The example image shows %Spectrum1DCanvas displaying a raw data layer and a peak data layer.

      @todo Use spectrum StringDataArray with name 'label' for peak annotations (Hiwi, Johannes)

      @htmlinclude OpenMS_Spectrum1DCanvas.parameters

      @ingroup SpectrumWidgets
  */
  class OPENMS_GUI_DLLAPI Spectrum1DCanvas :
    public SpectrumCanvas
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

    /// Default constructor
    Spectrum1DCanvas(const Param & preferences, QWidget * parent = nullptr);
    /// Destructor
    ~Spectrum1DCanvas() override;

    ///Enumerate all available paint styles
    enum DrawModes
    {
      DM_PEAKS,                                 ///< draw data as peak
      DM_CONNECTEDLINES                 ///< draw as connected lines
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
    bool mirrorModeActive();

    /// Sets whether this widget is currently in mirror mode
    void setMirrorModeActive(bool b);

    /// For convenience - calls dataToWidget
    void dataToWidget(const PeakType & peak, QPoint & point, bool flipped = false, bool percentage = true);

    /// Calls SpectrumCanvas::dataToWidget_(), takes mirror mode into account
    void dataToWidget(double x, double y, QPoint & point, bool flipped = false, bool percentage = false);

    /// For convenience - calls widgetToData
    PointType widgetToData(const QPoint & pos, bool percentage = false);

    /// Calls SpectrumCanvas::widgetToData_(), takes mirror mode into account
    PointType widgetToData(double x, double y, bool percentage = false);

    /// Display a static text box on the top right
    void setTextBox(const QString& html);

    /// ----- Annotations

    /// Add an annotation item for the given peak
    Annotation1DItem * addPeakAnnotation(const PeakIndex& peak_index, const QString& text, const QColor& color);

    /// Draws all annotation items of @p layer_index on @p painter
    void drawAnnotations(Size layer_index, QPainter & painter);

    /// ----- Alignment

    /// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
    void performAlignment(Size layer_index_1, Size layer_index_2, const Param & param);

    /// Resets alignment_
    void resetAlignment();

    /// Draws the alignment on @p painter
    void drawAlignment(QPainter & painter);

    /// Returns the number of aligned pairs of peaks
    Size getAlignmentSize();

    /// Returns the score of the alignment
    double getAlignmentScore();

    /// Returns aligned_peaks_indices_
    std::vector<std::pair<Size, Size> > getAlignedPeaksIndices();

    /// Sets current spectrum index of current layer to @p index
    void activateSpectrum(Size index, bool repaint = true);

    /// is the widget shown vertically? (for projections)
    void setSwappedAxis(bool swapped);

    /// Set's the Qt PenStyle of the active layer
    void setCurrentLayerPeakPenStyle(Qt::PenStyle ps);

    /// Actual painting takes place here
    void paint(QPainter * paint_device, QPaintEvent * e);
signals:
    /// Requests to display all spectra in 2D plot
    void showCurrentPeaksAs2D();
    /// Requests to display all spectra in 3D plot
    void showCurrentPeaksAs3D();

public slots:
    // Docu in base class
    void activateLayer(Size layer_index) override;
    // Docu in base class
    void removeLayer(Size layer_index) override;
    //docu in base class
    void updateLayer(Size i) override;

    /**
        @brief Sets the visible area.

        Sets the visible area to a new value. Note that it does not emit visibleAreaChanged()
        @param range the new visible area
    */
    void setVisibleArea(DRange<2> range);         //Do not change this to AreaType the signal needs QT needs the exact type...
    // Docu in base class
    void horizontalScrollBarChange(int value) override;

protected slots:

    /// Reacts on changed layer parameters
    void currentLayerParamtersChanged_();

protected:
    // Docu in base class
    bool finishAdding_() override;

    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawCoordinates_(QPainter & painter, const PeakIndex & peak);
    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawDeltas_(QPainter & painter, const PeakIndex & start, const PeakIndex & end);

    /**
        @brief Changes visible area interval

        This method is for convenience only. It calls changeVisibleArea_(const AreaType&, bool, bool) .
    */
    void changeVisibleArea_(double lo, double hi, bool repaint = true, bool add_to_stack = false);

    /// Draws a highlighted peak; if draw_elongation is true, the elongation line is drawn (for measuring)
    void drawHighlightedPeak_(Size layer_index, const PeakIndex & peak, QPainter & painter, bool draw_elongation = false);

    /// Draws a dashed line using the highlighted peak color parameter
    void drawDashedLine_(const QPoint & from, const QPoint & to, QPainter & painter);

    /// Recalculates the current scale factor based on the specified layer (= 1.0 if intensity mode != IM_PERCENTAGE)
    void updatePercentageFactor_(Size layer_index);

    /**
        @brief Sets the visible area

        Changes the visible area, adjusts the zoom stack and notifies interested clients about the change.
        If parts of the area are outside of the data area, the new area will be adjusted.

        @param new_area The new visible area.
        @param repaint if repainting of the widget should be triggered
        @param add_to_stack If the new area is to add to the zoom_stack_
    */
    void changeVisibleArea_(const AreaType & new_area, bool repaint = true, bool add_to_stack = false) override;
    // Docu in base class
    void recalculateSnapFactor_() override;
    // Docu in base class
    void updateScrollbars_() override;
    // Docu in base class
    void intensityModeChange_() override;

    /// Draw modes (for each spectrum)
    std::vector<DrawModes> draw_modes_;
    /// Draw style
    std::vector<Qt::PenStyle> peak_penstyle_;

    /// start point of "ruler" for measure mode
    QPoint measurement_start_point_;
    /// Indicates whether this widget is currently in mirror mode
    bool mirror_mode_;

    /// Indicates whether annotation items are just being moved on the canvas
    bool moving_annotations_;

    /// Indicates whether an alignment is currently visualized
    bool show_alignment_;
    /// Layer index of the first alignment layer
    Size alignment_layer_1_;
    /// Layer index of the second alignment layer
    Size alignment_layer_2_;
    /// Stores the alignment as MZ values of pairs of aligned peaks in both spectra
    std::vector<std::pair<double, double> > aligned_peaks_mz_delta_;
    /// Stores the peak indices of pairs of aligned peaks in both spectra
    std::vector<std::pair<Size, Size> > aligned_peaks_indices_;

    /// Stores the score of the last alignment
    double alignment_score_;
    /// is this widget showing data with swapped m/z and RT axis? (for drawCoordinates_ only)
    bool is_swapped_;

    /// Find peak next to the given position
    PeakIndex findPeakAtPosition_(QPoint);

    /// Shows dialog and calls addLabelAnnotation_
    void addUserLabelAnnotation_(const QPoint & screen_position);
    /// Adds an annotation item at the given screen position
    void addLabelAnnotation_(const QPoint & screen_position, QString label_text);
    /// Shows dialog and calls addPeakAnnotation_
    void addUserPeakAnnotation_(PeakIndex near_peak);

    /// Ensure that all annotations are within data range
    void ensureAnnotationsWithinDataRange_();

    QTextDocument text_box_content_;

    /** @name Reimplemented QT events */
    //@{
    void paintEvent(QPaintEvent * e) override;
    void mousePressEvent(QMouseEvent * e) override;
    void mouseReleaseEvent(QMouseEvent * e) override;
    void mouseMoveEvent(QMouseEvent * e) override;
    void keyPressEvent(QKeyEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    //@}

    ///Go forward in zoom history
    void zoomForward_() override;
    /// docu in base class
    void zoom_(int x, int y, bool zoom_in) override;
    //docu in base class
    void translateLeft_(Qt::KeyboardModifiers m) override;
    //docu in base class
    void translateRight_(Qt::KeyboardModifiers m) override;
    //docu in base class
    void paintGridLines_(QPainter & painter) override;
  };
} // namespace OpenMS

#endif
