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
// $Id: ScaleWidget.h,v 1.3 2006/03/28 08:03:26 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer:  $
// --------------------------------------------------------------------------
#ifndef OPENMS_VISUAL_SCALEWIDGET_H
#define OPENMS_VISUAL_SCALEWIDGET_H

#include <qwidget.h>
#include <qpainter.h>

namespace OpenMS
{

  /**
   a scale between 0 and 1 <br>
   0 is red, 1 is green <br>
  */
  class ScaleWidget : public QWidget
  {
    Q_OBJECT
  public:
    ScaleWidget ( QWidget* = 0,  const char* = 0);
    QSize sizeHint() const;
  public slots:
    void setMarker ( double );
    void reSize ( int );
  protected:
    void resizeEvent( QResizeEvent* );
    void paintEvent ( QPaintEvent * );
  private:
    QPainter painter_;
    double Marker_;
    QColor good_;
    QColor bad_;
    QColor markerC_;
    int width_,height_;
  };

}
#endif //OPENMS_VISUAL_SCALEWIDGET_H
