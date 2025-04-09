// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

// OpenMS
#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>
#include <OpenMS/KERNEL/PeakIndex.h>

// QT
class QPainter;
class QMouseEvent;
class QAction;
class QMenu;

namespace OpenMS
{
  /**
    @brief Canvas for 2D-visualization of peak map, feature map and consensus map data

    This widget displays a 2D representation of a set of peaks, features or consensus elements.

    @image html Plot2DCanvas.png

    The example image shows %Plot2DCanvas displaying a peak layer and a feature layer.

    @htmlinclude OpenMS_Plot2DCanvas.parameters

    @improvement Add RT interpolation mode for high zoom in 2D View (Hiwi)

    @improvement Snap also to min intensity (Hiwi)

    @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot2DCanvas :
    public PlotCanvas
  {
    Q_OBJECT

public:
    /// Default C'tor hidden
    Plot2DCanvas() = delete;

    /// Constructor
    Plot2DCanvas(const Param& preferences, QWidget* parent = nullptr);

    /// Destructor
    ~Plot2DCanvas() override;

    // Docu in base class
    void showCurrentLayerPreferences() override;

    /// Merges the features in @p map into the features layer @p i
    void mergeIntoLayer(Size i, const FeatureMapSharedPtrType& map);

    /// Merges the consensus features in @p map into the features layer @p i
    void mergeIntoLayer(Size i, const ConsensusMapSharedPtrType& map);

    /// Merges the peptide identifications in @p peptides into the peptide layer @p i
    void mergeIntoLayer(Size i, std::vector<PeptideIdentification>& peptides);

    /// recalculates the dot gradient of the active layer
    void recalculateCurrentLayerDotGradient();

signals:
    /// Requests to show projections for the @p source_layer. Emitted after calling pickProjectionLayer().
    void showProjections(const LayerDataBase* source_layer);
    /// Signal emitted when the projections are to be shown/hidden (e.g. via context menu)
    void toggleProjections();
    /// Requests to display the spectrum with index @p index in 1D
    void showSpectrumAsNew1D(int index);
    void showChromatogramsAsNew1D(std::vector<int, std::allocator<int> > indices);
    /// Requests to display all spectra in 3D plot
    void showCurrentPeaksAs3D();
    /// Requests to display this spectrum (=frame) in ion mobility plot
    void showCurrentPeaksAsIonMobility(const MSSpectrum& spec);


public slots:
    // Docu in base class
    void activateLayer(Size layer_index) override;
    // Docu in base class
    void removeLayer(Size layer_index) override;
    //docu in base class
    void updateLayer(Size i) override;
    // Docu in base class
    void horizontalScrollBarChange(int value) override;
    // Docu in base class
    void verticalScrollBarChange(int value) override;

    /// Picks an appropriate layer for projection and emits the signal showProjections(LayerDataBase*).
    void pickProjectionLayer();

protected slots:

    /// Reacts on changed layer parameters
    void currentLayerParametersChanged_();

protected:
    // Docu in base class
    bool finishAdding_() override;

    /// Collects fragment ion scans in the indicated RT/mz area and adds them to the indicated action
    bool collectFragmentScansInArea_(const RangeType& range, QAction* a, QMenu* msn_scans, QMenu* msn_meta);

    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawCoordinates_(QPainter& painter, const PeakIndex& peak);
    /// Draws the coordinates (or coordinate deltas) to the widget's upper left corner
    void drawDeltas_(QPainter& painter, const PeakIndex& start, const PeakIndex& end);

    /** @name Reimplemented QT events */
    //@{
    void mousePressEvent(QMouseEvent* e) override;
    void mouseReleaseEvent(QMouseEvent* e) override;
    void mouseMoveEvent(QMouseEvent* e) override;
    void paintEvent(QPaintEvent* e) override;
    void contextMenuEvent(QContextMenuEvent* e) override;
    void keyPressEvent(QKeyEvent* e) override;
    void keyReleaseEvent(QKeyEvent* e) override;
    void mouseDoubleClickEvent(QMouseEvent* e) override;
    //@}

    // Docu in base class
    void updateScrollbars_() override;

    // Docu in base class
    void intensityModeChange_() override;
    // Docu in base class
    void recalculateSnapFactor_() override;

    /**
      @brief Returns the position on color @p gradient associated with given intensity.

      Takes intensity modes into account.
    */
    Int precalculatedColorIndex_(float val, const MultiGradient& gradient, double snap_factor)
    {
      float gradientPos = val;
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
      }
      return gradient.precalculatedColorIndex(gradientPos);
    }

    /**
      @brief Returns the color associated with @p val for the gradient @p gradient.

      Takes intensity modes into account.
    */
    QColor heightColor_(float val, const MultiGradient& gradient, double snap_factor)
    {
      return gradient.precalculatedColorByIndex(precalculatedColorIndex_(val, gradient, snap_factor));
    }

    /** 
      @brief For a certain dimension: computes the size a data point would need, such that the image
             reaches a certain coverage
      
      Internally, this function makes use of the members 'canvas_coverage_min_' (giving the fraction (e.g. 20%) of area which should be covered by data)
      and 'pen_size_max_' (maximum allowed number of pixels per data point).

      @param ratio_data2pixel The current ratio of # data points vs. # pixels of image
      @param pen_size In/Out param: gives the initial pen size, and is increased (up to @p MAX_PEN_SIZE) to reach desired coverage given by 'canvas_coverage_min_'
      @return The factor by which @p pen_size increased (gives a hint of how many data points should be merged to avoid overplotting)
    */
    double adaptPenScaling_(double ratio_data2pixel, double& pen_size) const;
    
    /// Recalculates the dot gradient of a layer
    void recalculateDotGradient_(Size layer);

    /// Highlights a single peak and prints coordinates to screen
    void highlightPeak_(QPainter& p, const PeakIndex& peak);

    /// Returns the nearest peak to position @p pos
    PeakIndex findNearestPeak_(const QPoint& pos);

    /// Paints a peak icon for feature and consensus feature peaks
    void paintIcon_(const QPoint& pos, const QRgb& color, const String& icon, Size s, QPainter& p) const;

    /// Translates the visible area by a given offset specified in fractions of current visible area
    void translateVisibleArea_(double x_axis_rel, double y_axis_rel);

    /**
      @brief Convert chart to widget coordinates

      Translates chart coordinates to widget coordinates.
      @param x the chart coordinate x
      @param y the chart coordinate y
      @return A point in widget coordinates
    */
    QPoint dataToWidget_(double x, double y) const
    {
      QPoint point;
      const auto& xy = visible_area_.getAreaXY();
      point.setX(int((x - xy.minX()) / xy.width() * width()));
      point.setY(int((xy.maxY() - y) / xy.height() * height()));
      return point;
    }

    QPoint dataToWidget_(const DPosition<2>& xy)
    {
      return dataToWidget_(xy.getX(), xy.getY());
    }

    //docu in base class
    void translateLeft_(Qt::KeyboardModifiers m) override;
    //docu in base class
    void translateRight_(Qt::KeyboardModifiers m) override;
    //docu in base class
    void translateForward_() override;
    //docu in base class
    void translateBackward_() override;

    /// Finishes context menu after customization to peaks, features or consensus features
    void finishContextMenu_(QMenu* context_menu, QMenu* settings_menu);

    friend class Painter2DBase;
    friend class Painter2DChrom;
    friend class Painter2DConsensus;
    friend class Painter2DIdent;
    friend class Painter2DFeature;
    friend class Painter2DIonMobility;
    friend class Painter2DPeak;

    /// the nearest peak/feature to the mouse cursor
    PeakIndex selected_peak_;
    /// start peak/feature of measuring mode
    PeakIndex measurement_start_;

    /// stores the linear color gradient for non-log modes
    MultiGradient linear_gradient_;
    
    double pen_size_min_; ///< minimum number of pixels for one data point
    double pen_size_max_; ///< maximum number of pixels for one data point
    double canvas_coverage_min_; ///< minimum coverage of the canvas required; if lower, points are upscaled in size
  };
}

