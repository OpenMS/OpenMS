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
// $Id: SimMatrixWidget.h,v 1.4 2006/02/21 13:46:35 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_SIMMATRIXWIDGET_H
#define OPENMS_VISUAL_SIMMATRIXWIDGET_H
#include <qwidget.h>
#include <vector>
#include <map>
#include <qpainter.h>

namespace OpenMS
{
  /**
  shows a similarity matrix <br>
  standardcolors: green=similar, red=dissimilar<br>
  white=highlightcolumn, yellow=highlightentry<br>
   */
  class SimMatrixWidget : public QWidget
  {
    Q_OBJECT
  public:
    SimMatrixWidget( QWidget *parent = 0, const char* name = 0 );
    /**
    set the similarity matrix <br>
    values have to be between 0 and 1
    @param matr pointer to the similarity matrix
    */
    void setMatrix ( std::vector<std::vector<double > >* matr);
    ///set zoom factor
    void setZoomFactor ( int );

    /** @brief read accessor for zoomfactor<br> */
    int zoomFactor () const { return zoom; }
    QSize sizeHint () const;
    QSize minimumSize () const;
    QSize sizeIncrement () const;
  public slots:
    void resizeEvent ( QResizeEvent* );
    void highlight(int);
  protected:
    void mousePressEvent ( QMouseEvent * );
    void mouseMoveEvent ( QMouseEvent * );
    void mouseReleaseEvent( QMouseEvent * );
    void paintEvent ( QPaintEvent * );
    void keyPressEvent ( QKeyEvent * ); 
  private:
    void drawMatrixEntry ( QPainter* , int , int , QColor* = 0); 
    std::vector<std::vector<double> > * matrix;
    int zoom,lastx,lasty;
    QPainter painter;
    void activateCell ( int , int );
    void highlight ( int , int );
    QColor hlColor,hlColor2,gColor,bColor;
    std::vector<std::vector<int> > zoomstack_; 
  signals:
    void info ( int , int );
    void matrixEntry ( int , int );
    void score ( double );
  };
}
#endif //OPENMS_VISUAL_SIMMATRIXWIDGET_H
