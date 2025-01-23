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
#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/Plot2DCanvas.h>

class QGroupBox;
class QLabel;
class QCheckBox;

namespace OpenMS
{
  class Plot1DWidget;

  /**
      @brief Widget for 2D-visualization of peak map and feature map data

      The widget is composed of two scroll bars, two AxisWidget and a Plot2DCanvas as central widget.

      @image html Plot2DWidget.png

      The example image shows %Plot2DWidget displaying a peak layer and a feature layer.

      @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot2DWidget :
    public PlotWidget
  {
    Q_OBJECT
public:
    /// Main managed data type (experiment)
    typedef LayerDataBase::ExperimentSharedPtrType ExperimentSharedPtrType;

    /// Default constructor
    Plot2DWidget(const Param & preferences, QWidget * parent = nullptr);
    /// Destructor
    ~Plot2DWidget() override = default;

    // docu in base class
    Plot2DCanvas* canvas() const override
    {
      return static_cast<Plot2DCanvas*>(canvas_);
    }
            
    /// const reference to the horizontal projection
    const Plot1DWidget* getProjectionOntoX() const;
    /// const reference to the vertical projection
    const Plot1DWidget* getProjectionOntoY() const;

    /// Returns if one of the projections is visible (or both are visible)
    bool projectionsVisible() const;

    /// set the mapper for the canvas and the projections (the non-marginal projection
    /// axis will be mapped to Intensity).
    /// Also tries to guess drawing (sticks vs line) and intensity mode
    void setMapper(const DimMapper<2>& mapper) override
    {
      canvas_->setMapper(mapper); // update canvas
      // ... and projections: the projected Dim becomes intensity
      projection_onto_X_->setMapper(DimMapper<2>({mapper.getDim(DIM::X).getUnit(), DIM_UNIT::INT}));
      projection_onto_Y_->setMapper(DimMapper<2>({DIM_UNIT::INT, mapper.getDim(DIM::Y).getUnit()}));

          // decide on default draw mode, depending on main axis unit (e.g. m/z or RT)
      auto set_style = [&](const DIM_UNIT main_unit_1d, Plot1DCanvas* canvas) {
        switch (main_unit_1d)
        { // this may not be optimal for every unit. Feel free to change behavior.
          case DIM_UNIT::MZ:
            // to show isotope distributions as sticks
            canvas->setDrawMode(Plot1DCanvas::DM_PEAKS);
            canvas->setIntensityMode(PlotCanvas::IM_PERCENTAGE);
            break;
          // all other units
          default:
            canvas->setDrawMode(Plot1DCanvas::DM_CONNECTEDLINES);
            canvas->setIntensityMode(PlotCanvas::IM_SNAP);
            break;
        }
      };
      set_style(mapper.getDim(DIM::X).getUnit(), projection_onto_Y_->canvas());
      set_style(mapper.getDim(DIM::Y).getUnit(), projection_onto_X_->canvas());
    }

public slots:
    // Docu in base class
    void recalculateAxes_() override;
    /// Shows/hides the projections
    void toggleProjections();
    // Docu in base class
    void showGoToDialog() override;

signals:
    /**
        @brief Signal emitted whenever the visible area changes.

        @param area The new visible area.
    */
    void visibleAreaChanged(DRange<2> area);
    /// Requests to display the spectrum with index @p index in 1D
    void showSpectrumAsNew1D(int index);
    void showChromatogramsAsNew1D(std::vector<int, std::allocator<int> > indices);
    /// Requests to display all spectra as 1D
    void showCurrentPeaksAs3D();
    /// Requests to display this spectrum (=frame) in ion mobility plot
    void showCurrentPeaksAsIonMobility(const MSSpectrum& spec);


protected:
    /// shows projections information
    void projectionInfo_(int peaks, double intensity, double max);

    /// Vertical projection widget
    Plot1DWidget * projection_onto_X_;
    /// Horizontal projection widget
    Plot1DWidget * projection_onto_Y_;
    /// Group box that shows information about the projections
    QGroupBox * projection_box_;
    /// Number of peaks of the projection
    QLabel * projection_peaks_;
    /// Intensity sum of the projection
    QLabel * projection_sum_;
    /// Intensity maximum of the projection
    QLabel * projection_max_;
    /// Checkbox that indicates that projections should be automatically updated (with a slight delay)
    QCheckBox * projections_auto_;
    /// Timer that triggers auto-update of projections
    QTimer * projections_timer_;

private slots:
    /// extracts the projections from the @p source_layer and displays them
    void showProjections_(const LayerDataBase* source_layer);
    /// slot that monitors the visible area changes and triggers the update of projections
    void autoUpdateProjections_();
  };
}

