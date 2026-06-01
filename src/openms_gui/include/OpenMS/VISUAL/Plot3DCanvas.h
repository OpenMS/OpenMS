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

// OpenMS
#include <OpenMS/VISUAL/PlotCanvas.h>
#include <OpenMS/VISUAL/MultiGradient.h>


class QPainter;
class QOpenGLWidget;
class QResizeEvent;

namespace OpenMS
{
  class Plot3DOpenGLCanvas;

  /**
    @brief Canvas for 3D-visualization of peak map data

    The Plot3DCanvas uses the helper class Plot3DOpenGLCanvas for the
    actual 3D rendering.  Deriving Plot3DCanvas directly from QGLWidget is
    not possible due to the "Deadly Diamond" shape of inheritance.

    @image html Plot3DWidget.png

    @htmlinclude OpenMS_Plot3DCanvas.parameters

    @ingroup PlotWidgets
  */
  class OPENMS_GUI_DLLAPI Plot3DCanvas :
    public PlotCanvas
  {
    Q_OBJECT

    friend class Plot3DOpenGLCanvas;

public:

    /// Constructor
    Plot3DCanvas(const Param & preferences, QWidget * parent = nullptr);
    /// Destructor
     ~Plot3DCanvas() override;

    ///Different shade modes
    enum ShadeModes
    {
      SHADE_FLAT = 0,
      SHADE_SMOOTH = 1
    };

    ///returns the Plot3DOpenGLcanvas
    Plot3DOpenGLCanvas * openglwidget() const;

    ///@name Reimplemented Qt events
    //@{
    void resizeEvent(QResizeEvent * e) override;
    void contextMenuEvent(QContextMenuEvent * e) override;
    //@}
    /// Returns if the legend is shown
    bool isLegendShown() const;
    ///Shows/hides the legend
    void showLegend(bool);
    ///pointer to the SpectrumOpenGLCanvas implementation
    Plot3DOpenGLCanvas * openglcanvas_;

    // docu in base class
    void showCurrentLayerPreferences() override;

signals:

    /// Requests to display all spectra in 2D plot
    void showCurrentPeaksAs2D();

public slots:

    // Docu in base class
    void activateLayer(Size layer_index) override;
    // Docu in base class
    void removeLayer(Size layer_index) override;
    // Docu in base class
    void updateLayer(Size i) override;
    // Docu in base class
    void intensityModeChange_() override;

protected slots:

    /// Reacts on changed layer parameters
    void currentLayerParamtersChanged_();

protected:
    // Docu in base class
    bool finishAdding_() override;

    // Reimplementation in order to update the OpenGL widget
    void update_(const char * caller_name = nullptr) override;

    ///whether the legend is shown or not
    bool legend_shown_;

    ///stores the linear color gradient for non-log modes
    MultiGradient linear_gradient_;
  };

} //namespace
