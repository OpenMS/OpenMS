// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream>

//OpenMS
#include <OpenMS/VISUAL/Plot3DWidget.h>
#include <OpenMS/VISUAL/Plot3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Plot2DGoToDialog.h>


using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  Plot3DWidget::Plot3DWidget(const Param & preferences, QWidget * parent) :
    PlotWidget(preferences, parent)
  {
    setCanvas_(new Plot3DCanvas(preferences, this));

    x_axis_->hide();
    y_axis_->hide();

    // delegate signals from canvas
    connect(canvas(), SIGNAL(showCurrentPeaksAs2D()), this, SIGNAL(showCurrentPeaksAs2D()));
  }

  Plot3DWidget::~Plot3DWidget() = default;

  void Plot3DWidget::recalculateAxes_()
  {
  }

  void Plot3DWidget::showLegend(bool show)
  {
    canvas()->showLegend(show);
  }

  bool Plot3DWidget::isLegendShown() const
  {
    return static_cast<const Plot3DCanvas *>(canvas_)->isLegendShown();
  }

  void Plot3DWidget::showGoToDialog()
  {
    Plot2DGoToDialog goto_dialog(this, canvas_->getMapper().getDim(DIM::X).getDimNameShort(), canvas_->getMapper().getDim(DIM::Y).getDimNameShort());
    auto va = canvas()->getVisibleArea().getAreaXY();
    goto_dialog.setRange(va);

    auto all_area_xy = canvas_->getMapper().mapRange(canvas_->getDataRange());
    goto_dialog.setMinMaxOfRange(all_area_xy);
    
    goto_dialog.enableFeatureNumber(false);
    if (goto_dialog.exec())
    {
      canvas()->setVisibleArea(goto_dialog.getRange());
    }
  }

} //namespace
