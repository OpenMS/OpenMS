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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

// OpenMS includes
#include <OpenMS/VISUAL/ColorSelector.h>
#include <OpenMS/CONCEPT/Types.h>

//qt includes
#include <QtGui/QPainter>
#include <QtGui/QColorDialog>
#include <QtGui/QPaintEvent>
#include <QtGui/QMouseEvent>

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

  ColorSelector::~ColorSelector()
  {

  }

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
