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
// $Id: SimMatrixWidget.C,v 1.4 2006/03/29 12:30:29 andreas_bertsch Exp $
// $Author: andreas_bertsch $
// $Maintainer:  $
// --------------------------------------------------------------------------
#include <qpainter.h>
#include <OpenMS/VISUAL/SimMatrixWidget.h>
#include <iostream>

using namespace std;

namespace OpenMS
{

  SimMatrixWidget::SimMatrixWidget(QWidget *parent, const char* name)
    : QWidget(parent, name, WStaticContents),painter(),hlColor(Qt::yellow),hlColor2(Qt::white),gColor(Qt::green),bColor(Qt::red),zoomstack_()
  {
    zoom = 1;
    matrix = 0;
    setMouseTracking (1);
    setFocusPolicy(QWidget::StrongFocus);
    lastx = 0;
    lasty = 0;
    setSizePolicy(QSizePolicy(QSizePolicy::Preferred,QSizePolicy::Preferred));
  }

  void SimMatrixWidget::setMatrix(std::vector<std::vector < double > > * matr)
  {
    matrix = matr;
    if ( (int) matrix->size() < this->size().width()) setZoomFactor(this->size().width()/matrix->size());
    lasty = lastx = 0;
    update();
  }

  QSize SimMatrixWidget::sizeHint() const
  {
    QSize size ;
    if ( matrix != 0 ) size = zoom * QSize(matrix->size(),matrix->size());
    else size = QSize(100,100);
    return size;
  }

  void SimMatrixWidget::setZoomFactor(int newZoom)
  {
    if (newZoom < 1)
      newZoom = 1;

    if (newZoom != zoom) 
    {
      zoom = newZoom;
      update();
      updateGeometry();
    }
  }

  void SimMatrixWidget::paintEvent(QPaintEvent *)
  {
    painter.begin(this);
    painter.fillRect(0,0,size().width(),size().height(),QBrush(QColor(Qt::white)));
    if (matrix == 0 ) 
    {
      painter.drawLine(25,25,75,75);
      painter.drawLine(25,75,75,25);
      painter.end();
      return;
    }
    //int newzoom = minimum(size().width() / (ymax_-ymin_),size().width() / (xmax_-xmin_));
    //setZoomFactor(newzoom);
    for (uint i = 0; i < matrix->size(); ++i) 
    {
      for (uint j = 0; j < matrix->size(); ++j) 
        drawMatrixEntry(&painter, i, j);
    }
    painter.end();
  }

  void SimMatrixWidget::highlight(int row)
  {
    repaint();
    painter.begin(this); 
    QColor color;
    double value;
    for ( uint i = 0; i < matrix->size();++i)
    {
      value = (*matrix).at(row).at(i);
      color = QColor(
        (int) ( ((1-value) * bColor.red() + value * gColor.red())/2 + hlColor2.red()/2 )  ,
        (int) ( ((1-value) * bColor.green() +  value * gColor.green())/2 + hlColor2.green()/2 ) , 
        (int) ( ((1-value) * bColor.blue() + value*gColor.blue())/2 + hlColor2.blue()/2 ) 
        );
      drawMatrixEntry(&painter,row,i,&color);
    }
    painter.end();
  }

  void SimMatrixWidget::drawMatrixEntry(QPainter *painter, int i, int j, QColor* color)
  {
    double value = (*matrix).at(i).at(j);
    if ( !color ) color = new QColor( 
        (int) ( (1-value) * bColor.red() + value * gColor.red() )  ,
        (int) ( (1-value) * bColor.green() +  value * gColor.green() ) , 
        (int) ( (1-value) * bColor.blue() + value * gColor.blue() ) 
        );
    painter->fillRect(zoom * (i), zoom * (j),zoom, zoom, *color); 
  }

  void SimMatrixWidget::mouseMoveEvent(QMouseEvent *event)
  {
    if (!matrix) return;
    int x = event->pos().x()/zoom;
    int y = event->pos().y()/zoom;
    if ( x < 0 || (uint)x >= matrix->size() || y < 0 || (uint)y >= matrix->size()) return;
    else emit(info(x,y));
  }

  void SimMatrixWidget::mouseReleaseEvent(QMouseEvent* event)
  {
    if (!matrix) return;
    if (event->state() & Qt::ControlButton ) 
    {
      zoomstack_.push_back(std::vector<int>(4));
  //    (*zoomstack_.rbegin())[0] = xmin_;
  //    (*zoomstack_.rbegin())[1] = xmax_;
  //    (*zoomstack_.rbegin())[2] = ymin_;
  //    (*zoomstack_.rbegin())[3] = ymax_;
  //    displayMatrix(lastx,event->pos().x()/zoom+xmin_,lasty,event->pos().y()/zoom+ymin_);
    }
    else if (event->state() & Qt::ShiftButton && !zoomstack_.empty())
    {
      zoomstack_.pop_back();
    }
  }

  void SimMatrixWidget::mousePressEvent(QMouseEvent *event)
  {
    if (!matrix) return;
    uint x = event->pos().x()/zoom;
    uint y = event->pos().y()/zoom;
    if (x >= matrix->size() || y >= matrix->size()) return;
    activateCell(x,y);
  }

  QSize SimMatrixWidget::minimumSize() const 
  {
    QSize size ;
    if ( matrix != 0 ) size = QSize(matrix->size(),matrix->size() );
    else size = QSize(100,100);
    return size;
  }

  QSize SimMatrixWidget::sizeIncrement() const 
  {
    return QSize(matrix->size(), matrix->size());
  }

  void SimMatrixWidget::resizeEvent(QResizeEvent* event)
  {
    if (!event || !matrix ) return;
    uint maxresize = event->size().width();
    if ( maxresize  > matrix->size() )
    {
      setZoomFactor(maxresize/matrix->size());
    }
    else setZoomFactor(1);
  }

  void SimMatrixWidget::keyPressEvent ( QKeyEvent * key )
  {
    if ( !matrix ) return;
    int x = lastx;
    int y = lasty;
    if ( key->key() == Qt::Key_Up )
    {
      if ( lasty > 0 ) y--;
    }
    else if ( key->key() == Qt::Key_Down )
    {
      if ( lasty < (int)matrix->size() -1 ) y++;
    }
    else if ( key->key() == Qt::Key_Left )
    {
      if ( lastx > 0 ) x--;
    }
    else if ( key->key() == Qt::Key_Right )
    {
      if ( lastx <  (int)matrix->size() -1 ) x++;
    }
    else 
    {
      key->ignore();
      return;
    }
    key->accept();
    activateCell(x,y);
  }

  void SimMatrixWidget::activateCell(int x, int y )
  {
    highlight(x,y);
    emit score((*matrix)[x][y]);
    emit matrixEntry(x,y);//translate
    emit info(x,y);//translate
  }

  void SimMatrixWidget::highlight(int x, int y)
  {
    painter.begin(this);
    drawMatrixEntry(&painter,lastx,lasty);
    drawMatrixEntry(&painter,x,y,&hlColor);
    lastx = x; 
    lasty = y;
    painter.end();
    painter.flush();
  }

}


