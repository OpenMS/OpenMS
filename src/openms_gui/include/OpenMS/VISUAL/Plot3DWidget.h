// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Cornelia Friedle $
// --------------------------------------------------------------------------

#pragma once

// OpenMS_GUI config
#include <OpenMS/VISUAL/OpenMS_GUIConfig.h>

#include <OpenMS/VISUAL/PlotWidget.h>
#include <OpenMS/VISUAL/Plot3DCanvas.h>

namespace OpenMS
{
  class Plot3DCanvas;
  /**
      @brief Widget for 3D-visualization of map data

      @image html Plot3DWidget.png

      @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot3DWidget :
    public PlotWidget
  {
    Q_OBJECT

public:
    /// Constructor
    Plot3DWidget(const Param & preferences, QWidget * parent = nullptr);

    /// Destructor
    ~Plot3DWidget() override;

    // docu in base class
    Plot3DCanvas* canvas() const override
    {
      return static_cast<Plot3DCanvas*>(canvas_);
    }

    // Docu in base class
    void setMapper(const DimMapper<2>& /*mapper*/) override
    { // 3D widget currently only handles MSExperiment. That's it.
      throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
    }

    // Docu in base class
    void recalculateAxes_() override;

    //docu in base class
    bool isLegendShown() const override;
    //docu in base class
    void showLegend(bool show) override;

signals:
    /// Requests to display all spectra in 2D plot
    void showCurrentPeaksAs2D();

public slots:
    // Docu in base class
    void showGoToDialog() override;
  };

} //namespace

