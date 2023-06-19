// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/PlotWidget.h>

#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/HistogramDialog.h>
#include <OpenMS/VISUAL/DIALOGS/LayerStatisticsDialog.h>
#include <OpenMS/VISUAL/VISITORS/LayerStatistics.h>

#include <QCloseEvent>
#include <QtCore/QMimeData>
#include <QtWidgets/QFileDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/QScrollBar>
#include <QScrollBar>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  const char PlotWidget::RT_AXIS_TITLE[] = "Time [s]";
  const char PlotWidget::MZ_AXIS_TITLE[] = "m/z";
  const char PlotWidget::INTENSITY_AXIS_TITLE[] = "Intensity";
  const char PlotWidget::IM_MS_AXIS_TITLE[] = "Ion Mobility [ms]";
  const char PlotWidget::IM_ONEKZERO_AXIS_TITLE[] = "Ion Mobility [1/K0]";

  PlotWidget::PlotWidget(const Param& /*preferences*/, QWidget* parent) :
    QWidget(parent),
    canvas_(nullptr)
  {
    setAttribute(Qt::WA_DeleteOnClose);
    grid_ = new QGridLayout(this);
    grid_->setSpacing(0);
    grid_->setMargin(1);
    // axes
    y_axis_ = new AxisWidget(AxisPainter::LEFT, "", this);
    x_axis_ = new AxisWidget(AxisPainter::BOTTOM, "", this);
    // scrollbars
    x_scrollbar_ = new QScrollBar(Qt::Horizontal, this); // left is small value, right is high value
    y_scrollbar_ = new QScrollBar(Qt::Vertical, this);   // top is low value, bottom is high value (however, our coordinate system is inverse for the y-Axis!)
    // We achieve the desired behavior by setting negative min/max ranges within the scrollbar when updateVScrollbar() is called.
    // Thus, y_scrollbar_->valueChanged() will report negative values (which you need to multiply by -1 to get the correct value).
    // Remember this when implementing verticalScrollBarChange() in your canvas class (currently only used in Plot2DCanvas)
    // Do NOT use the build-in functions to invert a scrollbar, since implementation can be incomplete depending on style and platform
    // y_scrollbar_->setInvertedAppearance(true);
    // y_scrollbar_->setInvertedControls(true);

    setMinimumSize(250, 250); //Canvas (200) + AxisWidget (30) + ScrollBar (20)
    setSizePolicy(QSizePolicy::MinimumExpanding, QSizePolicy::MinimumExpanding);

    setAcceptDrops(true);
  }

  void PlotWidget::setCanvas_(PlotCanvas* canvas, UInt row, UInt col)
  {
    canvas_ = canvas;
    setFocusProxy(canvas_);
    grid_->addWidget(canvas_, row, col);
    grid_->addWidget(y_axis_, row, col - 1);
    grid_->addWidget(x_axis_, row + 1, col);
    connect(canvas_, &PlotCanvas::visibleAreaChanged, this, &PlotWidget::updateAxes);
    connect(canvas_, &PlotCanvas::recalculateAxes, this, &PlotWidget::updateAxes);
    connect(canvas_, &PlotCanvas::changeLegendVisibility, this, &PlotWidget::changeLegendVisibility);

    grid_->addWidget(y_scrollbar_, row, col - 2);
    grid_->addWidget(x_scrollbar_, row + 2, col);
    x_scrollbar_->hide();
    y_scrollbar_->hide();
    connect(canvas_, &PlotCanvas::updateHScrollbar, this, &PlotWidget::updateHScrollbar);
    connect(canvas_, &PlotCanvas::updateVScrollbar, this, &PlotWidget::updateVScrollbar);
    connect(x_scrollbar_, &QScrollBar::valueChanged, canvas_, &PlotCanvas::horizontalScrollBarChange);
    connect(y_scrollbar_, &QScrollBar::valueChanged, canvas_, &PlotCanvas::verticalScrollBarChange);
    connect(canvas_, &PlotCanvas::sendStatusMessage, this, &PlotWidget::sendStatusMessage);
    connect(canvas_, &PlotCanvas::sendCursorStatus, this, &PlotWidget::sendCursorStatus);

    canvas_->setPlotWidget(this);
  }

  PlotWidget::~PlotWidget() = default;

  Int PlotWidget::getActionMode() const
  {
    return canvas_->getActionMode();
  }

  void PlotWidget::setIntensityMode(PlotCanvas::IntensityModes mode)
  {
    if (canvas_->getIntensityMode() != mode)
    {
      canvas_->setIntensityMode(mode);
    }
  }

  void PlotWidget::showStatistics()
  {
    LayerStatisticsDialog lsd(this, canvas_->getCurrentLayer().getStats());
    lsd.exec();
  }

  void PlotWidget::showIntensityDistribution(const Histogram<>& dist)
  {
    HistogramDialog dw(dist);
    dw.setLegend(PlotWidget::INTENSITY_AXIS_TITLE);
    dw.setLogMode(true);
    if (dw.exec() == QDialog::Accepted)
    {
      DataFilters filters;

      if (dw.getLeftSplitter() > dist.minBound())
      {
        DataFilters::DataFilter filter;
        filter.value = dw.getLeftSplitter();
        filter.field = DataFilters::INTENSITY;
        filter.op = DataFilters::GREATER_EQUAL;
        filters.add(filter);
      }

      if (dw.getRightSplitter() < dist.maxBound())
      {
        DataFilters::DataFilter filter;
        filter.value = dw.getRightSplitter();
        filter.field = DataFilters::INTENSITY;
        filter.op = DataFilters::LESS_EQUAL;
        filters.add(filter);
      }

      canvas_->setFilters(filters);
    }
  }

  void PlotWidget::showMetaDistribution(const String& name, const Histogram<>& dist)
  {
    HistogramDialog dw(dist);
    dw.setLegend(name);

    if (dw.exec() == QDialog::Accepted)
    {
      DataFilters filters;

      if (dw.getLeftSplitter() > dist.minBound())
      {
        DataFilters::DataFilter filter;
        filter.value = dw.getLeftSplitter();
        filter.field = DataFilters::META_DATA;
        filter.meta_name = name;
        filter.op = DataFilters::GREATER_EQUAL;
        filter.value_is_numerical = true;
        filters.add(filter);
      }

      if (dw.getRightSplitter() < dist.maxBound())
      {
        DataFilters::DataFilter filter;
        filter.value = dw.getRightSplitter();
        filter.field = DataFilters::META_DATA;
        filter.meta_name = name;
        filter.op = DataFilters::LESS_EQUAL;
        filter.value_is_numerical = true;
        filters.add(filter);
      }

      canvas_->setFilters(filters);
    }
  }

  void PlotWidget::showLegend(bool show)
  {
    y_axis_->showLegend(show);
    x_axis_->showLegend(show);
    update();
  }

  void PlotWidget::saveAsImage()
  {
    QString file_name = QFileDialog::getSaveFileName(this, "Save File", "", "Images (*.bmp *.png *.jpg *.gif)");
    QString old_stylesheet = this->styleSheet();
    // Make the whole widget (including the usually natively styled AxisWidgets) white
    this->setStyleSheet("background: white");
    bool x_visible = x_scrollbar_->isVisible();
    bool y_visible = y_scrollbar_->isVisible();
    x_scrollbar_->hide();
    y_scrollbar_->hide();
    QPixmap pixmap = this->grab();
    x_scrollbar_->setVisible(x_visible);
    y_scrollbar_->setVisible(y_visible);
    pixmap.save(file_name);
    this->setStyleSheet(old_stylesheet);
  }

  void PlotWidget::updateAxes()
  {
    recalculateAxes_();
  }

  void PlotWidget::intensityModeChange_()
  {

  }

  bool PlotWidget::isLegendShown() const
  {
    // Both are shown or hidden, so we simply return the state of the x-axis
    return x_axis_->isLegendShown();
  }

  void PlotWidget::hideAxes()
  {
    y_axis_->hide();
    x_axis_->hide();
  }

  void updateScrollbar(QScrollBar* scroll, float f_min, float disp_min, float disp_max, float f_max)
  {
    if ((disp_min == f_min && disp_max == f_max) || (disp_min < f_min && disp_max > f_max))
    {
      scroll->hide();
    }
    else
    {
      // block signals as this causes repainting due to rounding (QScrollBar works with int ...)
      auto local_min = min(f_min, disp_min);
      auto local_max = max(f_max, disp_max);
      auto vis_span = disp_max - disp_min;
      scroll->blockSignals(true);
      //scroll->setMinimum(static_cast<int>(local_min));
      //scroll->setMaximum(static_cast<int>(std::ceil(local_max - disp_max + disp_min)));
      scroll->setRange(int(local_min), int(std::ceil(local_max - vis_span)));
      scroll->setValue(int(disp_min)); // emits valueChanged, which will call this function here unless signal is blocked
      scroll->setPageStep(vis_span);
      scroll->blockSignals(false);
      scroll->show();
    }
  }

  void PlotWidget::updateHScrollbar(float f_min, float disp_min, float disp_max, float f_max)
  {
    updateScrollbar(x_scrollbar_, f_min, disp_min, disp_max, f_max);
  }

  void PlotWidget::updateVScrollbar(float f_min, float disp_min, float disp_max, float f_max)
  {
    updateScrollbar(y_scrollbar_, f_min, disp_min, disp_max, f_max);
  }

  void PlotWidget::changeLegendVisibility()
  {
    showLegend(!isLegendShown());
  }

  void PlotWidget::closeEvent(QCloseEvent* e)
  {
    for (UInt l = 0; l < canvas()->getLayerCount(); ++l)
    {
      //modified => ask if it should be saved
      const LayerDataBase& layer = canvas()->getLayer(l);
      if (layer.modified)
      {
        QMessageBox::StandardButton result = QMessageBox::question(this, "Save?", (String("Do you want to save your changes to layer '") + layer.getName() +  "'?").toQString(), QMessageBox::Ok | QMessageBox::Discard);
        if (result == QMessageBox::Ok)
        {
          canvas()->activateLayer(l);
          canvas()->saveCurrentLayer(false);
        }
      }
    }
    e->accept();
  }

  void PlotWidget::dragEnterEvent(QDragEnterEvent* event)
  {
    if (event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
    }
  }

  void PlotWidget::dragMoveEvent(QDragMoveEvent* event)
  {
    if (event->mimeData()->hasUrls())
    {
      event->acceptProposedAction();
    }
  }

  void PlotWidget::dropEvent(QDropEvent* event)
  {
    emit dropReceived(event->mimeData(), dynamic_cast<QWidget*>(event->source()), this->getWindowId());
    event->acceptProposedAction();
  }

  void PlotWidget::paintEvent(QPaintEvent * /*event*/)
  {
    QStyleOption opt;
    opt.init(this);
    QPainter p(this);
    // apply style options and draw the widget using current stylesheets
    style()->drawPrimitive(QStyle::PE_Widget, &opt, &p, this);
  }

} //namespace OpenMS
