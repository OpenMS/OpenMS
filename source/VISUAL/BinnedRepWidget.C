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

#include <OpenMS/VISUAL/BinnedRepWidget.h>
#include <OpenMS/COMPARISON/CLUSTERING/BinnedRep.h>

#include <qpixmap.h>

#include <iostream>

namespace OpenMS
{

  BinnedRepWidget::BinnedRepWidget(QWidget* parent , const char* name) 
  : QWidget(parent,name),binrep_(0),width_(1),height_(1),buffer_(),draw_upside_down_(0),left_padding_(0)
  {
    resize(parent->width(), parent->height());
    updateBuffer_(parent->width(),parent->height());
    setMouseTracking(1);
  }

  void BinnedRepWidget::paintEvent(QPaintEvent* )
  {
    painter_.begin(buffer_);
    painter_.setWindow(QRect(0,0,width_,height_));
    painter_.fillRect(0, 0, width_, height_, QColor(255, 255, 255));
    painter_.setBrush(QBrush(Qt::blue));
    painter_.setPen(QPen(Qt::darkBlue,1,SolidLine));
    ToolTips.clear();
    if ( binrep_ )
    {
      for ( uint i = 0; i < binrep_->size(); ++i)
      {
        if ( !draw_upside_down_)
        {
          QRect rect = painter_.xForm(QRect(i*10+left_padding_,110,10,(int)(-(*binrep_)[i]*100)).normalize());
          painter_.drawRect(i*10+2+left_padding_,110,6,(int)(-(*binrep_)[i]*100));
          ToolTips.push_back(rect);
        }
        else 
        {
          QRect rect = painter_.xForm(QRect(i*10+left_padding_,10,10,(int)((*binrep_)[i]*100)).normalize());
          painter_.drawRect(i*10+2+left_padding_,10,6,(int)((*binrep_)[i]*100)); 
          ToolTips.push_back(rect);
        }
      }
    }
    painter_.end();
    bitBlt(this, 0 , 0, buffer_);
  }


  void BinnedRepWidget::loadBinnedRep(const BinnedRep& data)
  {
    width_ = data.size()*10;
    height_ = 120;
    binrep_ = new BinnedRep(data);
    repaint();
  }

  void BinnedRepWidget::loadSpectrum(const MSSpectrum< DPeak<1> >& spectrum)
  {
    binrep_ = new BinnedRep(10,3);
    (*binrep_) << spectrum;
    width_ = binrep_->size()*10;
    height_ = 120;
    repaint();
  }

  void BinnedRepWidget::resizeEvent(QResizeEvent* e)
  {
    updateBuffer_(e->size().width(),e->size().height());
  }

  void BinnedRepWidget::updateBuffer_(int width, int height)
  {
    if (buffer_ == 0)
    {
      buffer_ = new QPixmap(width, height);
      buffer_->setOptimization(QPixmap::BestOptim);
      buffer_->fill();
    } 
    else
    {
      buffer_->resize(width, height);
      buffer_->fill();
    }
  }

  void BinnedRepWidget::draw_upside_down(bool b )
  {
    draw_upside_down_ = b;
  }

  //only sets the Area larger than the binrep!
  //only usefull for compatible BinnedReps (i.e. same binsize)
  void BinnedRepWidget::setVisibleArea(double a, double b)
  {
    if ( !binrep_ ) return;
    if ( binrep_->min() > a )
    {
      left_padding_ = ((int)( binrep_->min()/binrep_->getBinSize() - a/binrep_->getBinSize() +0.9 ) )*10; // +1 means round up , todo better
      width_+= left_padding_;
    }
    else left_padding_ = 0;
    if ( binrep_->max() < b )
    {
      width_+=  ((int) (b/binrep_->getBinSize() - binrep_->max()/binrep_->getBinSize() ))*10;
    }
    repaint();
  }
  void BinnedRepWidget::mouseMoveEvent ( QMouseEvent * e )
  {
    for ( uint i = 0; i < ToolTips.size() ; ++i )
    {
      if ( ToolTips[i].contains(e->pos()) )
      {
        emit(posInfo(QString("x: %1 ").arg(binrep_->min()+i*binrep_->getBinSize())));
        emit(boxInfo(QString("normalized height: %1").arg((*binrep_)[i])));
      }
    }
    if (binrep_ ) emit(binInfo(QString("Bin: spectrum.id: %1 Binsize: %2 BinSpread: %3").arg(binrep_->id()).arg(binrep_->getBinSize()).arg(binrep_->getBinSpread())));
  }

}

