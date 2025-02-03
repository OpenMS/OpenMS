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

// STL
#include <vector>

// OpenMS
#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/Plot1DCanvas.h>

class QAction;
class QSpacerItem;

namespace OpenMS
{
  class Plot1DCanvas;

  /**
      @brief Widget for visualization of several spectra

      The widget is composed of a scroll bar, an AxisWidget and a Plot1DCanvas as central widget.

      @image html Plot1DWidget.png

      The example image shows %Plot1DWidget displaying a raw data layer and a peak data layer.

      @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot1DWidget :
    public PlotWidget
  {
    Q_OBJECT

public:    
    /// Default constructor
    Plot1DWidget(const Param& preferences, const DIM gravity_axis = DIM::Y, QWidget* parent = nullptr);
    /// Destructor
    ~Plot1DWidget() override;

    // docu in base class
    Plot1DCanvas* canvas() const override
    {
      return static_cast<Plot1DCanvas*>(canvas_);
    }

    // Docu in base class
    void setMapper(const DimMapper<2>& mapper) override
    {
      canvas_->setMapper(mapper);
    }

    // Docu in base class
    void hideAxes() override;

    // Docu in base class
    void showLegend(bool show) override;

    /// Switches to mirror view, displays another y-axis for the second spectrum
    void toggleMirrorView(bool mirror);
    
    /// Performs an alignment of the layers with @p layer_index_1 and @p layer_index_2
    void performAlignment(Size layer_index_1, Size layer_index_2, const Param & param);

    /// Resets the alignment
    void resetAlignment();

    // Docu in base class
    void saveAsImage() override;

    // Docu in base class
    virtual void renderForImage(QPainter& painter);

signals:
    /// Requests to display the whole spectrum in 2D view
    void showCurrentPeaksAs2D();

    /// Requests to display the whole spectrum in 3D view
    void showCurrentPeaksAs3D();

    /// Requests to display the whole spectrum in ion mobility view
    void showCurrentPeaksAsIonMobility(const MSSpectrum& spec);

    /// Requests to display a full DIA window
    void showCurrentPeaksAsDIA(const Precursor& pc, const MSExperiment& exp);

public slots:
    // Docu in base class
    void showGoToDialog() override;

protected:
    // Docu in base class
    void recalculateAxes_() override;

    /// The second y-axis for the mirror view
    AxisWidget * flipped_y_axis_;

    /// Spacer between the two y-axes in mirror mode (needed when visualizing an alignment)
    QSpacerItem * spacer_;

  };
} // namespace OpenMS

