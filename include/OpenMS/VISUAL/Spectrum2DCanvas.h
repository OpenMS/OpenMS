// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_SPECTRUM2DCANVAS_H
#define OPENMS_VISUAL_SPECTRUM2DCANVAS_H

// OpenMS
#include <OpenMS/VISUAL/SpectrumCanvas.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>
#include <OpenMS/KERNEL/PeakIndex.h>

// QT
class QPainter;
class QMouseEvent;

namespace OpenMS
{
  /**
  @brief Canvas for 2D-visualization of peak map, feature map and consensus map data

  This widget displays a 2D representation of a set of peaks, features or consensus elements.

      @image html Spectrum2DCanvas.png

      The example image shows %Spectrum2DCanvas displaying a peak layer and a feature layer.

      @htmlinclude OpenMS_Spectrum2DCanvas.parameters

      @improvement Add RT interpolation mode for high zoom in 2D View (Hiwi)

      @improvement Snap also to min intensity (Hiwi)

  @ingroup SpectrumWidgets
*/
  class OPENMS_GUI_DLLAPI Spectrum2DCanvas :
    public SpectrumCanvas
  {
    Q_OBJECT

public:
    /// Default constructor
    Spectrum2DCanvas(const Param & preferences, QWidget * parent = 0);

    /// Destructor
    ~Spectrum2DCanvas();

    // Docu in base class
    virtual void showCurrentLayerPreferences();

    // Docu in base class
    virtual void saveCurrentLayer(bool visible);

    /// Merges the features in @p map into the features layer @p i
    void mergeIntoLayer(Size i, FeatureMapSharedPtrType map);

    /// Merges the consensus features in @p map into the features layer @p i
    void mergeIntoLayer(Size i, ConsensusMapSharedPtrType map);

    /// Merges the peptide identifications in @p peptides into the peptide layer @p i
    void mergeIntoLayer(Size i, std::vector<PeptideIdentification> & peptides);

    /// recalculates the dot gradient of the active layer
    void recalculateCurrentLayerDotGradient();

signals:
    /// Sets the data for the horizontal projection
    void showProjectionHorizontal(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes);
    /// Sets the data for the vertical projection
    void showProjectionVertical(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes);
    /// Shows the number of peaks and the intensity sum of the projection
    void showProjectionInfo(int, double, double);
    /// Signal emitted when the projections are to be shown/hidden
    void toggleProjections();
    /// Requests to display the spectrum with index @p index in 1D
    void showSpectrumAs1D(int index);
    void showSpectrumAs1D(std::vector<int, std::allocator<int> > indices);
    /// Requests to display all spectra in 3D plot
    void showCurrentPeaksAs3D();

public slots:
    // Docu in base class
    void activateLayer(Size layer_index);
    // Docu in base class
    void removeLayer(Size layer_index);
    //docu in base class
    virtual void updateLayer(Size i);
    // Docu in base class
    virtual void horizontalScrollBarChange(int value);
    // Docu in base class
    virtual void verticalScrollBarChange(int value);
    /**
    @brief Updates the projection data and emits some related signals.

    Emitted signals are showProjectionHorizontal(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes) and
    showProjectionVertical(ExperimentSharedPtrType, Spectrum1DCanvas::DrawModes).

    @see projection_mz_
    @see projection_rt_
  */
    void updateProjections();

protected slots:

    /// Reacts on changed layer paramters
    void currentLayerParametersChanged_();

protected:
    // Docu in base class
    bool finishAdding_();

    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawCoordinates_(QPainter & painter, const PeakIndex & peak);
    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawDeltas_(QPainter & painter, const PeakIndex & start, const PeakIndex & end);

    /** @name Reimplemented QT events */
    //@{
    void mousePressEvent(QMouseEvent * e);
    void mouseReleaseEvent(QMouseEvent * e);
    void mouseMoveEvent(QMouseEvent * e);
    void paintEvent(QPaintEvent * e);
    void contextMenuEvent(QContextMenuEvent * e);
    void keyPressEvent(QKeyEvent * e);
    void keyReleaseEvent(QKeyEvent * e);
    void mouseDoubleClickEvent(QMouseEvent * e);
    //@}

    // Docu in base class
    virtual void updateScrollbars_();

    /**
    @brief Paints individual peaks.

    Calls different painting methods depending on the layer type and the density of displayed peaks

    @param layer_index The index of the layer.
    @param p The QPainter to paint on.
  */
    void paintDots_(Size layer_index, QPainter & p);

    void paintAllIntensities_(Size layer_index, DoubleReal pen_width, QPainter & painter);

    /**
    @brief Paints maximum intensity of individual peaks.

    Paints the peaks as small ellipses. The peaks are colored according to the
    selected dot gradient.

    @param layer_index The index of the layer.
    @param p The QPainter to paint on.
  */
    void paintMaximumIntensities_(Size layer_index, Size rt_pixel_count, Size mz_pixel_count, QPainter & p);

    /**
    @brief Paints the precursor peaks.

    @param layer_index The index of the layer.
    @param painter The QPainter to paint on.
   */
    void paintPrecursorPeaks_(Size layer_index, QPainter & painter);

    /**
    @brief Paints feature data.

    @param layer_index The index of the layer.
    @param p The QPainter to paint on.
   */
    void paintFeatureData_(Size layer_index, QPainter & p);

    /**
    @brief Paints convex hulls (one for each mass trace) of a features layer.

    @param layer_index Index of the layer.
    @param p The QPainter to paint on.
  */
    void paintTraceConvexHulls_(Size layer_index, QPainter & p);

    /**
    @brief Paints the convex hulls (one for each feature) of a features layer.

    @param layer_index Index of the layer.
    @param p The QPainter to paint on.
  */
    void paintFeatureConvexHulls_(Size layer_index, QPainter & p);

    /**
    @brief Paints peptide identifications (for idXML and unassigned peptides in featureXML).

    @param layer_index Index of the layer.
    @param p The QPainter to paint on.
  */
    void paintIdentifications_(Size layer_index, QPainter & p);

    /**
    @brief Paints the consensus elements of a consensus features layer.

    @param layer_index Index of the layer.
    @param p The QPainter to paint on.
  */
    void paintConsensusElements_(Size layer_index, QPainter & p);

    /**
    @brief Paints one consensus element of a consensus features layer.

    @param layer_index Index of the layer.
    @param cf Reference to the feature to be painted.
    @param p The QPainter to paint on.
    @param use_buffer Flag to switch between painting on the buffer and screen.
  */
    void paintConsensusElement_(Size layer_index, const ConsensusFeature & cf, QPainter & p, bool use_buffer);

    /**
    @brief checks if any element of a consensus feature is currently visible.

    @param layer_index Index of the layer.
    @param ce The ConsensusFeature that needs checking
  */
    bool isConsensusFeatureVisible_(const ConsensusFeature & ce, Size layer_index);

    /**
    @brief Paints convex hulls (one for each mass trace) for a single feature.

    @param hulls Reference to convex hull vector.
    @param p The QPainter to paint on.
  */
    void paintConvexHulls_(const std::vector<ConvexHull2D> & hulls, bool hasIdentifications, QPainter & p);

    // Docu in base class
    virtual void intensityModeChange_();
    // DOcu in base class
    virtual void recalculateSnapFactor_();

    /// m/z projection data
    ExperimentType projection_mz_;
    /// RT projection data
    ExperimentType projection_rt_;

    /**
    @brief Returns the position on color @p gradient associated with given intensity @p.

    Takes intensity modes into account.
  */
    inline Int precalculatedColorIndex_(Real val, const MultiGradient & gradient, DoubleReal snap_factor)
    {
      Real gradientPos;
      switch (intensity_mode_)
      {
      case IM_NONE:
        gradientPos = val;
        break;
      case IM_PERCENTAGE:
        gradientPos = val * percentage_factor_;
        break;
      case IM_SNAP:
        gradientPos = val * snap_factor;
        break;
      case IM_LOG:
        gradientPos = std::log(val + 1);
        break;
      default:
        throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
      }
      return gradient.precalculatedColorIndex( gradientPos );
    }

    /**
    @brief Returns the color associated with @p val for the gradient @p gradient.

    Takes intensity modes into account.
    */
    inline QColor heightColor_(Real val, const MultiGradient & gradient, DoubleReal snap_factor)
    {
      return gradient.precalculatedColorByIndex( precalculatedColorIndex_( val, gradient, snap_factor ) );
    }

    /**
      @brief Convert chart to widget coordinates

      Translates chart coordinates to widget coordinates.
      @param x the chart coordinate x
      @param y the chart coordinate y
      @param point returned widget coordinates
    */
    inline void dataToWidget_(DoubleReal x, DoubleReal y, QPoint & point)
    {
      if (!isMzToXAxis())
      {
        point.setX(int((y - visible_area_.minY()) / visible_area_.height() * width()));
        point.setY(height() - int((x - visible_area_.minX()) / visible_area_.width() * height()));
      }
      else
      {
        point.setX(int((x - visible_area_.minX()) / visible_area_.width() * width()));
        point.setY(height() - int((y - visible_area_.minY()) / visible_area_.height() * height()));
      }
    }

    /// recalculates the dot gradient of a layer
    void recalculateDotGradient_(Size layer);

    /// Highlights a single peak and prints coordinates to screen
    void highlightPeak_(QPainter & p, const PeakIndex & peak);

    /// Returns the nearest peak to position @p pos
    PeakIndex findNearestPeak_(const QPoint & pos);

    /// Paints a peak icon for feature and consensus feature peaks
    void paintIcon_(const QPoint & pos, const QRgb & color, const String & icon, Size s, QPainter & p) const;

    /// the nearest peak/feature to the mouse cursor
    PeakIndex selected_peak_;
    /// start peak/feature of measuring mode
    PeakIndex measurement_start_;

    /// translates the visible area by a given offset specified in fractions of current visible area
    virtual void translateVisibleArea_( DoubleReal mzShiftRel, DoubleReal rtShiftRel );

    //docu in base class
    virtual void translateLeft_();
    //docu in base class
    virtual void translateRight_();
    //docu in base class
    virtual void translateForward_();
    //docu in base class
    virtual void translateBackward_();

    /// Finishes context menu after customization to peaks, features or consensus features
    void finishContextMenu_(QMenu * context_menu, QMenu * settings_menu);
  };
}

#endif
