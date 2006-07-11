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

#ifndef OPENMS_VISUAL_BINNEDREPWIDGET_H
#define OPENMS_VISUAL_BINNEDREPWIDGET_H

#include <qwidget.h>
#include <qpainter.h>
#include <qpixmap.h>

#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include<vector>

namespace OpenMS
{

  class BinnedRepWidget : public QWidget
  {
    Q_OBJECT
  public:
    BinnedRepWidget(QWidget* parent = 0, const char* name = 0);
    void loadBinnedRep(const BinnedRep& );
    void loadSpectrum(const MSSpectrum< DPeak<1> >&);
    void draw_upside_down(bool b = true);
    void setVisibleArea(double,double);
  protected:
    void paintEvent( QPaintEvent* );
    void resizeEvent( QResizeEvent* );
    void mouseMoveEvent ( QMouseEvent * );
  private:
    void updateBuffer_(int,int);

    std::vector<QRect> ToolTips; // QToolTips are not easy to get rid off
    BinnedRep* binrep_;
    QPainter painter_;
    unsigned int width_;
    unsigned int height_;
    QPixmap* buffer_;
    bool draw_upside_down_;
    int left_padding_;

  signals:
    void boxInfo(const QString&);
    void posInfo(const QString&);
    void binInfo(const QString&);
  };

}
#endif //OPENMS_VISUAL_BINNEDREPWIDGET_H
