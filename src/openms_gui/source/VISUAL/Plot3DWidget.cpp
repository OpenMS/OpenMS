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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <iostream>

//OpenMS
#include <OpenMS/VISUAL/Plot3DWidget.h>
#include <OpenMS/VISUAL/Plot3DOpenGLCanvas.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Plot2DGoToDialog.h>

//QT
#include <QtWidgets/QGridLayout>

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

  Plot3DWidget::~Plot3DWidget()
  {

  }

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
    Plot2DGoToDialog goto_dialog(this);
    const DRange<3> & area = canvas()->getDataRange();
    goto_dialog.setRange(area.minY(), area.maxY(), area.minX(), area.maxX());
    goto_dialog.setMinMaxOfRange(canvas()->getDataRange().minY(), canvas()->getDataRange().maxY(), canvas()->getDataRange().minX(), canvas()->getDataRange().maxX());
    goto_dialog.enableFeatureNumber(false);
    if (goto_dialog.exec())
    {
      goto_dialog.fixRange(); // in case user did something invalid
      PlotCanvas::AreaType area (goto_dialog.getMinMZ(), goto_dialog.getMinRT(), goto_dialog.getMaxMZ(), goto_dialog.getMaxRT());
      correctAreaToObeyMinMaxRanges_(area);
      canvas()->setVisibleArea(area);
    }
  }

} //namespace
