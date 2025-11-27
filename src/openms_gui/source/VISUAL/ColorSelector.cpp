// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/CONCEPT/Types.h>

//qt includes
#include <QPainter>
#include <QtWidgets/QColorDialog>
#include <QPaintEvent>
#include <QMouseEvent>

using namespace std;

namespace OpenMS
{

  ColorSelector::ColorSelector(QWidget * parent) :
    QWidget(parent),
    color_(255, 255, 255)
  {
    setSizePolicy(QSizePolicy::Fixed, QSizePolicy::Fixed);
  }

  QSize ColorSelector::sizeHint() const
  {
    return QSize(15, 15);
  }

  ColorSelector::~ColorSelector() = default;

  void ColorSelector::paintEvent(QPaintEvent * /*e*/)
  {
    Int size = std::min(width(), height());
    QPainter painter(this);
    painter.setPen(QColor(0, 0, 0));
    painter.drawRect(0, 0, size - 1, size - 1);
    painter.setPen(QColor(255, 255, 255));
    painter.drawRect(1, 1, size - 3, size - 3);

    painter.fillRect(2, 2, size - 4, size - 4, color_);
  }

  void ColorSelector::mousePressEvent(QMouseEvent * e)
  {
    if (e->button() != Qt::LeftButton)
    {
      e->ignore();
      return;
    }
    QColor tmp = QColorDialog::getColor(color_, this);
    if (tmp.isValid())
    {
      color_ = tmp;
      repaint();
    }
  }

  const QColor & ColorSelector::getColor()
  {
    return color_;
  }

  void ColorSelector::setColor(const QColor & col)
  {
    color_ = col;
    repaint();
  }

} //namespace OpenMS
