// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <OpenMS/VISUAL/ScaleWidget.h>
#include <iostream>

#include <qtooltip.h>

namespace OpenMS
{

  ScaleWidget::ScaleWidget(QWidget* parent, const char* name)
    : QWidget(parent,name),painter_(),Marker_(0.5),good_(Qt::green),bad_(Qt::red),markerC_(Qt::white),width_(0),height_(0)
  {
  }

  void ScaleWidget::setMarker(double x)
  {
    if ( x <= 1 && x >= 0 )
    {
      QToolTip::remove(this);
      Marker_ = x;
      repaint();
      QToolTip::add(this,QString("%1").arg(x));
    }
  }

  void ScaleWidget::reSize(int width)
  {
    width_ = width;
    height_ = width/3;
    update();
    updateGeometry();
  }

  void ScaleWidget::paintEvent(QPaintEvent* )
  {
    painter_.begin(this);
    painter_.setPen(Qt::NoPen);
    int steps = width_/4;
    painter_.setBrush(QBrush(good_));
    painter_.drawRect(0,0,width_,height_);
    for ( int i = 0; i < steps-1; ++i )
    {
      double fraction = (double)i/steps;
      painter_.setBrush(QBrush(QColor(
            (int) ( (1-fraction) * bad_.red() + fraction * good_.red() )  ,
            (int) ( (1-fraction) * bad_.green() +  fraction * good_.green() ) ,
            (int) ( (1-fraction) * bad_.blue() + fraction * good_.blue() ) ) ) );
      painter_.drawRect(width_/steps*i,0,width_/steps,height_);
    }
    painter_.setBrush(QBrush(QColor(markerC_)));
    QPointArray arrows = QPointArray(6);
    arrows.setPoint(0,(int)(Marker_*width_+width_/50),0);
    arrows.setPoint(1,(int)(Marker_*width_-width_/50),0);
    arrows.setPoint(2,(int)(Marker_*width_),(int)(height_/5));
    painter_.drawPolygon(arrows,1,0,3);
    arrows.setPoint(3,(int)(Marker_*width_+width_/50),height_);
    arrows.setPoint(4,(int)(Marker_*width_-width_/50),height_);
    arrows.setPoint(5,(int)(Marker_*width_),(int)(4*height_/5));
    painter_.drawPolygon(arrows,0,3,3);
    painter_.end();
  }

  void ScaleWidget::resizeEvent(QResizeEvent* e)
  {
    reSize(e->size().width());
  }

  QSize ScaleWidget::sizeHint() const
  {
    return QSize(width_,height_);
  }

}
