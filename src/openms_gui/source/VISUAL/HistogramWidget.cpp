// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


// Qt
#include <QResizeEvent>
#include <QMouseEvent>
#include <QPaintEvent>
#include <QPainter>
#include <QtWidgets/QMenu>

// STL
#include <iostream>

// OpenMS
#include <OpenMS/VISUAL/HistogramWidget.h>
#include <OpenMS/VISUAL/AxisWidget.h>

using namespace std;

namespace OpenMS
{
  using namespace Math;

  HistogramWidget::HistogramWidget(const Histogram<> & distribution, QWidget * parent) :
    QWidget(parent),
    dist_(distribution),
    show_splitters_(false),
    moving_splitter_(0),
    margin_(30),
    buffer_(),
    log_mode_(false)
  {
    left_splitter_ =  dist_.minBound();
    right_splitter_ = dist_.maxBound();
    setSizePolicy(QSizePolicy::Minimum, QSizePolicy::Minimum);
    setMinimumSize(600, 450);
    bottom_axis_ = new AxisWidget(AxisPainter::BOTTOM, "", this);
    bottom_axis_->setMargin(margin_);
    bottom_axis_->setTickLevel(2);
    bottom_axis_->setAxisBounds(dist_.minBound(), dist_.maxBound());

    //signals and slots
    setContextMenuPolicy(Qt::CustomContextMenu);
    connect(this, SIGNAL(customContextMenuRequested(const QPoint &)), this, SLOT(showContextMenu(const QPoint &)));
  }

  HistogramWidget::~HistogramWidget()
  {
    delete(bottom_axis_);
  }

  double HistogramWidget::getLeftSplitter() const
  {
    return left_splitter_;
  }

  double HistogramWidget::getRightSplitter() const
  {
    return right_splitter_;
  }

  void HistogramWidget::showSplitters(bool on)
  {
    show_splitters_ = on;
  }

  void HistogramWidget::setRightSplitter(double pos)
  {
    right_splitter_ = min(dist_.maxBound(), pos);
  }

  void HistogramWidget::setLeftSplitter(double pos)
  {
    left_splitter_ = max(dist_.minBound(), pos);
  }

  void HistogramWidget::setLegend(const String & legend)
  {
    bottom_axis_->setLegend(legend);
  }

  void HistogramWidget::mousePressEvent(QMouseEvent * e)
  {
    if (show_splitters_ && e->button() == Qt::LeftButton)
    {
      //left
      Int p = margin_ + UInt(((left_splitter_ - dist_.minBound()) / (dist_.maxBound() - dist_.minBound())) * (width() - 2 * margin_));
      //cout << "Mouse: " << e->x() << " p: " << p << " splitter: " << left_splitter_ << endl;
      if (e->x() >= p && e->x() <= p + 5)
      {
        moving_splitter_ = 1;
      }

      //right
      p = margin_ + UInt(((right_splitter_ - dist_.minBound()) / (dist_.maxBound() - dist_.minBound())) * (width() - 2 * margin_));
      if (e->x() <= p && e->x() >= p - 5)
      {
        moving_splitter_ = 2;
      }
    }
    else
    {
      e->ignore();
    }
  }

  void HistogramWidget::mouseMoveEvent(QMouseEvent * e)
  {
    if (show_splitters_ && (e->buttons() & Qt::LeftButton))
    {
      //left
      if (moving_splitter_ == 1)
      {
        left_splitter_ = double(Int(e->x()) - Int(margin_)) / (width() - 2 * margin_) * (dist_.maxBound() - dist_.minBound()) + dist_.minBound();
        //upper bound
        if (left_splitter_ > right_splitter_ - (dist_.maxBound() - dist_.minBound()) / 50.0)
        {
          left_splitter_ = right_splitter_ - (dist_.maxBound() - dist_.minBound()) / 50.0;
        }
        //lower bound
        if (left_splitter_ < dist_.minBound())
        {
          left_splitter_ = dist_.minBound();
        }
        update();
      }

      //right
      if (moving_splitter_ == 2)
      {

        right_splitter_ = double(Int(e->x()) - Int(margin_)) / (width() - 2 * margin_ + 2) * (dist_.maxBound() - dist_.minBound()) + dist_.minBound();
        //upper bound
        if (right_splitter_ < left_splitter_ + (dist_.maxBound() - dist_.minBound()) / 50.0)
        {
          right_splitter_ = left_splitter_ + (dist_.maxBound() - dist_.minBound()) / 50.0;
        }
        //lower bound
        if (right_splitter_ > dist_.maxBound())
        {
          right_splitter_ = dist_.maxBound();
        }
        update();
      }
    }
    else
    {
      e->ignore();
    }
  }

  void HistogramWidget::mouseReleaseEvent(QMouseEvent * e)
  {
    if (show_splitters_)
    {
      moving_splitter_ = 0;
    }
    else
    {
      e->ignore();
    }
  }

  void HistogramWidget::paintEvent(QPaintEvent * /*e*/)
  {
    //histogram from buffer
    QPainter painter2(this);
    painter2.drawPixmap(margin_, 0, buffer_);

    //y-axis label
    painter2.rotate(270);
    painter2.setPen(Qt::black);
    QString label = "count";
    if (log_mode_)
    {
      label = "log ( count )";
    }
    painter2.drawText(0, 0, -height(), margin_, Qt::AlignHCenter | Qt::AlignVCenter, label);
    painter2.end();

    //draw splitters
    if (show_splitters_)
    {
      QPainter painter(this);
      painter.setPen(Qt::black);
      QFont label_font;
      label_font.setPointSize(8);

      //cout << "Left splitter: " << left_splitter_<< " dist: " << dist_.minBound() << endl;
      //cout << "Right splitter: " << right_splitter_<< " dist: " << dist_.maxBound() << endl;

      //left
      UInt p =  UInt(((left_splitter_ - dist_.minBound()) / (dist_.maxBound() - dist_.minBound())) * (width() - 2 * margin_)) + margin_;
      //cout << "Left splitter position: " << p << endl;
      painter.drawLine(p, margin_ - 8, p, height() - bottom_axis_->height());
      painter.drawLine(p, margin_ - 8, p + 5, margin_ - 8);
      painter.drawLine(p + 5, margin_ - 8, p, margin_ - 3);
      painter.setFont(label_font);
      painter.drawText(p, margin_ - 8, "lower boundary");
      painter.setFont(QFont());

      //right
      p = UInt(((right_splitter_ - dist_.minBound()) / (dist_.maxBound() - dist_.minBound())) * (width() - 2 * margin_)) + margin_;
      painter.drawLine(p, margin_ - 8, p, height() - bottom_axis_->height());
      painter.drawLine(p, margin_ - 8, p - 5, margin_ - 8);
      painter.drawLine(p - 5, margin_ - 8, p, margin_ - 3);
      painter.setFont(label_font);
      painter.drawText(p, margin_ - 8, "upper boundary");
      painter.setFont(QFont());
    }
  }

  void HistogramWidget::resizeEvent(QResizeEvent * /*e*/)
  {
    buffer_ = QPixmap(width() - margin_, height() - bottom_axis_->height());
    bottom_axis_->setGeometry(margin_, height() - bottom_axis_->height(), width() - margin_, bottom_axis_->height());
    invalidate_();
  }

  void HistogramWidget::invalidate_()
  {
    //apply log trafo if needed
    Math::Histogram<> dist(dist_);
    if (log_mode_)
    {
      dist.applyLogTransformation(100.0);
    }

    QPainter painter(&buffer_);
    buffer_.fill(palette().window().color());
    UInt w = buffer_.width();
    UInt h = buffer_.height();
    UInt pen_width = std::min(margin_, UInt(0.5 * w / dist.size()));

    //draw distribution
    QPen pen;
    pen.setWidth(pen_width);
    pen.setColor(QColor(100, 125, 175));
    painter.setPen(pen);

    for (Size i = 0; i < dist.size(); ++i)
    {
      if (dist[i] != 0)
      {
        UInt bin_pos = UInt((double(i) / (dist.size() - 1)) * (w - margin_));
        UInt bin_height = UInt(((double)dist[i] / dist.maxValue()) * (h - margin_));
        painter.drawLine(bin_pos + 1, h, bin_pos + 1, h - bin_height);
      }
    }

    //calculate total intensity
    double total_sum = 0;
    for (Size i = 0; i < dist.size(); ++i)
    {
      total_sum += dist[i];
    }

    // draw part of total intensity
    painter.setPen(Qt::red);
    QPoint last_point(1, h);
    QPoint point;
    double int_sum = 0;
    for (Size i = 0; i < dist.size(); ++i)
    {
      int_sum += dist[i];
      point.setX(UInt((double(i) / (dist.size() - 1)) * (w - margin_)));
      point.setY(UInt((1 - (int_sum / total_sum)) * (h - margin_) + margin_));
      painter.drawLine(last_point, point);
      last_point = point;
    }
    // draw coord system (on top distribution)
    painter.setPen(Qt::black);
    painter.drawLine(0, h - 1, w - margin_ + Int(0.5 * pen_width), h - 1);

    update();
  }

  void HistogramWidget::showContextMenu(const QPoint & pos)
  {
    //create menu
    QMenu menu(this);
    QAction * action = menu.addAction("Normal mode");
    if (!log_mode_)
    {
      action->setEnabled(false);
    }
    action = menu.addAction("Log mode");
    if (log_mode_)
    {
      action->setEnabled(false);
    }
    //execute
    QAction * result = menu.exec(mapToGlobal(pos));
    //change according to selected value
    if (result != nullptr)
    {
      if (result->text() == "Normal mode")
      {
        setLogMode(false);
      }
      else if (result->text() == "Log mode")
      {
        setLogMode(true);
      }
    }
  }

  void HistogramWidget::setLogMode(bool log_mode)
  {
    log_mode_ = log_mode;
    if (!buffer_.isNull())
    {
      invalidate_();
    }
  }

} //namespace OpenMS
