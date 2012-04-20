// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Timo Sachsenberg $
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/Spectrum1DWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>
#include <OpenMS/VISUAL/DIALOGS/Spectrum1DGoToDialog.h>
#include <QtGui/QSpacerItem>
#include <QtGui/QScrollBar>
#include <QtGui/QFileDialog>
#include <QtGui/QPainter>
#include <QtGui/QPaintEvent>
#include <QtSvg/QtSvg>
#include <QtSvg/QSvgGenerator>

using namespace std;

namespace OpenMS
{
  using namespace Internal;
  using namespace Math;

  Spectrum1DWidget::Spectrum1DWidget(const Param & preferences, QWidget * parent) :
    SpectrumWidget(preferences, parent)
  {
    //set the label mode for the axes  - side effect
    setCanvas_(new Spectrum1DCanvas(preferences, this));

    x_axis_->setLegend(String(Peak2D::shortDimensionName(Peak2D::MZ)) + " [" + String(Peak2D::shortDimensionUnit(Peak2D::MZ)) + "]");
    x_axis_->setAllowShortNumbers(false);
    y_axis_->setLegend("Intensity");
    y_axis_->setAllowShortNumbers(true);
    y_axis_->setMinimumWidth(50);

    flipped_y_axis_ = new AxisWidget(AxisPainter::LEFT, "Intensity", this);
    flipped_y_axis_->setInverseOrientation(true);
    flipped_y_axis_->setAllowShortNumbers(true);
    flipped_y_axis_->setMinimumWidth(50);
    flipped_y_axis_->hide();

    spacer_ = new QSpacerItem(0, 0);

    //Delegate signals
    connect(canvas(), SIGNAL(showCurrentPeaksAs2D()), this, SIGNAL(showCurrentPeaksAs2D()));
    connect(canvas(), SIGNAL(showCurrentPeaksAs3D()), this, SIGNAL(showCurrentPeaksAs3D()));
  }

  void Spectrum1DWidget::recalculateAxes_()
  {
    AxisWidget * mz_axis;
    AxisWidget * it_axis;

    //determine axes
    if (canvas()->isMzToXAxis())
    {
      mz_axis = x_axis_;
      it_axis = y_axis_;
    }
    else
    {
      mz_axis = y_axis_;
      it_axis = x_axis_;
    }

    // recalculate gridlines
    mz_axis->setAxisBounds(canvas()->getVisibleArea().minX(), canvas()->getVisibleArea().maxX());
    switch (canvas()->getIntensityMode())
    {
    case SpectrumCanvas::IM_NONE:
      if (it_axis->isLogScale())
      {
        it_axis->setLogScale(false);
        flipped_y_axis_->setLogScale(false);
      }

      it_axis->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
      flipped_y_axis_->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
      break;

    case SpectrumCanvas::IM_PERCENTAGE:
      if (it_axis->isLogScale())
      {
        it_axis->setLogScale(false);
        flipped_y_axis_->setLogScale(false);
      }

      it_axis->setAxisBounds(canvas()->getVisibleArea().minY() / canvas()->getDataRange().maxY() * 100.0, canvas()->getVisibleArea().maxY() / canvas()->getDataRange().maxY() * 100.0);
      flipped_y_axis_->setAxisBounds(canvas()->getVisibleArea().minY() / canvas()->getDataRange().maxY() * 100.0, canvas()->getVisibleArea().maxY() / canvas()->getDataRange().maxY() * 100.0);
      break;

    case SpectrumCanvas::IM_SNAP:
      if (it_axis->isLogScale())
      {
        it_axis->setLogScale(false);
        flipped_y_axis_->setLogScale(false);
      }

      it_axis->setAxisBounds(canvas()->getVisibleArea().minY() / canvas()->getSnapFactor(), canvas()->getVisibleArea().maxY() / canvas()->getSnapFactor());
      flipped_y_axis_->setAxisBounds(canvas()->getVisibleArea().minY() / canvas()->getSnapFactor(), canvas()->getVisibleArea().maxY() / canvas()->getSnapFactor());
      break;

    case SpectrumCanvas::IM_LOG:
      if (!it_axis->isLogScale())
      {
        it_axis->setLogScale(true);
        flipped_y_axis_->setLogScale(true);
      }

      it_axis->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
      flipped_y_axis_->setAxisBounds(canvas()->getVisibleArea().minY(), canvas()->getVisibleArea().maxY());
      break;

    default:
      throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
    }
  }

  Histogram<> Spectrum1DWidget::createIntensityDistribution_() const
  {
    //initialize histogram
    DoubleReal min = canvas_->getCurrentMinIntensity();
    DoubleReal max = canvas_->getCurrentMaxIntensity();
    if (min == max)
    {
      min -= 0.01;
      max += 0.01;
    }
    Histogram<> tmp(min, max, (max - min) / 500.0);

    for (ExperimentType::SpectrumType::ConstIterator it = (*canvas_->getCurrentLayer().getPeakData())[0].begin(); it != (*canvas_->getCurrentLayer().getPeakData())[0].end(); ++it)
    {
      tmp.inc(it->getIntensity());
    }
    return tmp;
  }

  Histogram<> Spectrum1DWidget::createMetaDistribution_(const String & name) const
  {
    Histogram<> tmp;
    //float arrays
    const ExperimentType::SpectrumType::FloatDataArrays & f_arrays = (*canvas_->getCurrentLayer().getPeakData())[0].getFloatDataArrays();
    for (ExperimentType::SpectrumType::FloatDataArrays::const_iterator it = f_arrays.begin(); it != f_arrays.end(); ++it)
    {
      if (it->getName() == name)
      {
        //determine min and max of the data
        Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
        for (Size i = 0; i < it->size(); ++i)
        {
          if ((*it)[i] < min)
            min = (*it)[i];
          if ((*it)[i] > max)
            max = (*it)[i];
        }
        if (min >= max)
          return tmp;

        //create histogram
        tmp.reset(min, max, (max - min) / 500.0);
        for (Size i = 0; i < it->size(); ++i)
        {
          tmp.inc((*it)[i]);
        }
      }
    }
    //integer arrays
    const ExperimentType::SpectrumType::IntegerDataArrays & i_arrays = (*canvas_->getCurrentLayer().getPeakData())[0].getIntegerDataArrays();
    for (ExperimentType::SpectrumType::IntegerDataArrays::const_iterator it = i_arrays.begin(); it != i_arrays.end(); ++it)
    {
      if (it->getName() == name)
      {
        //determine min and max of the data
        Real min = numeric_limits<Real>::max(), max = -numeric_limits<Real>::max();
        for (Size i = 0; i < it->size(); ++i)
        {
          if ((*it)[i] < min)
            min = (*it)[i];
          if ((*it)[i] > max)
            max = (*it)[i];
        }
        if (min >= max)
          return tmp;

        //create histogram
        tmp.reset(min, max, (max - min) / 500.0);
        for (Size i = 0; i < it->size(); ++i)
        {
          tmp.inc((*it)[i]);
        }
      }
    }
    //fallback if no array with that name exists
    return tmp;
  }

  Spectrum1DWidget::~Spectrum1DWidget()
  {

  }

  void Spectrum1DWidget::showGoToDialog()
  {
    Spectrum1DGoToDialog goto_dialog(this);
    goto_dialog.setRange(canvas()->getVisibleArea().minX(), canvas()->getVisibleArea().maxX());
    if (goto_dialog.exec())
    {
      canvas()->setVisibleArea(SpectrumCanvas::AreaType(goto_dialog.getMin(), 0, goto_dialog.getMax(), 0));
    }
  }

  void Spectrum1DWidget::showLegend(bool show)
  {
    y_axis_->showLegend(show);
    flipped_y_axis_->showLegend(show);
    x_axis_->showLegend(show);
    update();
  }

  void Spectrum1DWidget::hideAxes()
  {
    y_axis_->hide();
    flipped_y_axis_->hide();
    x_axis_->hide();
  }

  void Spectrum1DWidget::toggleMirrorView(bool mirror)
  {
    if (mirror)
    {
      grid_->addItem(spacer_, 1, 1);
      grid_->addWidget(flipped_y_axis_, 2, 1);
      grid_->removeWidget(canvas());
      grid_->removeWidget(x_axis_);
      grid_->removeWidget(x_scrollbar_);
      grid_->addWidget(canvas(), 0, 2, 3, 1);       // rowspan = 3
      grid_->addWidget(x_axis_, 3, 2);
      grid_->addWidget(x_scrollbar_, 4, 2);
      flipped_y_axis_->show();
    }
    else
    {
      grid_->removeWidget(canvas());
      grid_->removeWidget(flipped_y_axis_);
      flipped_y_axis_->hide();
      grid_->removeItem(spacer_);
      grid_->removeWidget(x_axis_);
      grid_->removeWidget(x_scrollbar_);
      grid_->addWidget(canvas(), 0, 2);
      grid_->addWidget(x_axis_, 1, 2);
      grid_->addWidget(x_scrollbar_, 2, 2);
    }
  }

  void Spectrum1DWidget::performAlignment(Size layer_index_1, Size layer_index_2, const Param & param)
  {
    spacer_->changeSize(0, 10);
    grid_->removeWidget(y_axis_);
    grid_->removeWidget(flipped_y_axis_);
    grid_->addWidget(y_axis_, 0, 1);
    grid_->addWidget(flipped_y_axis_, 2, 1);

    canvas()->performAlignment(layer_index_1, layer_index_2, param);
  }

  void Spectrum1DWidget::resetAlignment()
  {
    spacer_->changeSize(0, 0);
    grid_->removeWidget(y_axis_);
    grid_->removeWidget(flipped_y_axis_);
    grid_->addWidget(y_axis_, 0, 1);
    grid_->addWidget(flipped_y_axis_, 2, 1);
  }

  void Spectrum1DWidget::saveAsImage()
  {
    QString file_name = QFileDialog::getSaveFileName(this, "Save File", "", "Raster images *.bmp *.png *.jpg *.gif (*.bmp *.png *.jpg *.gif);;Vector images *.svg (*.svg)");
    bool x_visible = x_scrollbar_->isVisible();
    bool y_visible = y_scrollbar_->isVisible();
    x_scrollbar_->hide();
    y_scrollbar_->hide();

    if (file_name.contains(".svg", Qt::CaseInsensitive)) // svg vector format
    {
      QSvgGenerator generator;
      generator.setFileName(file_name);
      generator.setSize(QSize(this->width(), this->height()));
      generator.setViewBox(QRect(0, 0, this->width() - 1, this->height() - 1));
      generator.setTitle(file_name);
      generator.setDescription("TOPPView generated SVG");
      QPainter painter;
      painter.begin(&generator);

      painter.save();
      painter.translate(QPoint(y_axis_->pos()));
      y_axis_->paint(&painter, new QPaintEvent(y_axis_->contentsRect()));
      painter.restore();

      painter.save();
      painter.translate(QPoint(canvas_->pos()));
      dynamic_cast<Spectrum1DCanvas *>(canvas_)->paint(&painter, new QPaintEvent(canvas_->contentsRect()));
      painter.restore();

      painter.save();
      painter.translate(QPoint(x_axis_->pos()));
      x_axis_->paint(&painter, new QPaintEvent(x_axis_->contentsRect()));
      painter.restore();

      painter.end();
      x_scrollbar_->setVisible(x_visible);
      y_scrollbar_->setVisible(y_visible);
    }
    else   // raster graphics formats
    {
      QPixmap pixmap = QPixmap::grabWidget(this);
      x_scrollbar_->setVisible(x_visible);
      y_scrollbar_->setVisible(y_visible);
      pixmap.save(file_name);
    }
  }

} //namespace
