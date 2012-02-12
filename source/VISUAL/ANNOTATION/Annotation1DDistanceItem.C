// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/ANNOTATION/Annotation1DDistanceItem.h>
#include <OpenMS/VISUAL/Spectrum1DCanvas.h>

#include <QtCore/QPoint>
#include <QtGui/QPainter>

using namespace std;

namespace OpenMS
{	

	Annotation1DDistanceItem::Annotation1DDistanceItem(const QString& text, const PointType& start_point, const PointType& end_point)
		: Annotation1DItem(text),
			start_point_(start_point),
			end_point_(end_point)
	{
	}
	
	Annotation1DDistanceItem::Annotation1DDistanceItem(const Annotation1DDistanceItem& rhs)
		: Annotation1DItem(rhs)
	{
		start_point_ = rhs.getStartPoint();
		end_point_ = rhs.getEndPoint();
	}
	
	Annotation1DDistanceItem::~Annotation1DDistanceItem()
	{
	}
	
	void Annotation1DDistanceItem::draw(Spectrum1DCanvas* const canvas, QPainter& painter, bool flipped)
	{
		//translate mz/intensity to pixel coordinates
		QPoint start_p, end_p;
		canvas->dataToWidget(start_point_.getX(), start_point_.getY(), start_p, flipped, true);
		canvas->dataToWidget(end_point_.getX(), end_point_.getY(), end_p, flipped, true);
		
		// compute bounding box on the specified painter
		if (canvas->isMzToXAxis())
		{
			bounding_box_ = QRectF(QPointF(start_p.x(), start_p.y()), QPointF(end_p.x(), end_p.y()+4)); // +4 for lower half of arrow heads
		}
		else
		{
			bounding_box_ = QRectF(QPointF(start_p.x()-4, start_p.y()), QPointF(end_p.x(), end_p.y()));
		}
		
		// find out how much additional space is needed for the text:
		QRectF text_boundings = painter.boundingRect(QRectF(), Qt::AlignCenter, text_);
		if (canvas->isMzToXAxis())
		{
			bounding_box_.setTop(bounding_box_.top() - text_boundings.height());
		}
		else
		{
			bounding_box_.setRight(bounding_box_.right() + text_boundings.width());
		}
		// if text doesn't fit between peaks, enlarge bounding box:
		if (canvas->isMzToXAxis())
		{
			if (text_boundings.width() > bounding_box_.width())
			{
				float additional_space = (text_boundings.width() - bounding_box_.width()) / 2;
				bounding_box_.setLeft(bounding_box_.left() - additional_space);
				bounding_box_.setRight(bounding_box_.right() + additional_space);
			}
		}
		else
		{
			if (text_boundings.height() > bounding_box_.height())
			{
				float additional_space = (text_boundings.height() - bounding_box_.height()) / 2;
				bounding_box_.setTop(bounding_box_.top() - additional_space);
				bounding_box_.setBottom(bounding_box_.bottom() + additional_space);
			}
		}
		
		// draw line
		painter.drawLine(start_p, end_p);

    // draw ticks
    if ( !ticks_.empty() )
    {
      for(vector<DoubleReal>::iterator it = ticks_.begin(); it != ticks_.end(); ++it)
      {
        QPoint tick;
        canvas->dataToWidget(*it, start_point_.getY(), tick, flipped, true);
        painter.drawLine(tick.x(), tick.y()-4, tick.x(), tick.y()+4);
      }
    }
		
    // draw arrow heads and the ends if they won't overlap
    if ((start_p - end_p).manhattanLength() > 10)
    {
      if (canvas->isMzToXAxis())
      {
        painter.drawLine(start_p, QPoint(start_p.x()+5, start_p.y()-4));
        painter.drawLine(start_p, QPoint(start_p.x()+5, start_p.y()+4));
        painter.drawLine(end_p, QPoint(end_p.x()-5, end_p.y()-4));
        painter.drawLine(end_p, QPoint(end_p.x()-5, end_p.y()+4));
      }
      else
      {
        painter.drawLine(start_p, QPoint(start_p.x()+4, start_p.y()-5));
        painter.drawLine(start_p, QPoint(start_p.x()-4, start_p.y()-5));
        painter.drawLine(end_p, QPoint(end_p.x()+4, end_p.y()+5));
        painter.drawLine(end_p, QPoint(end_p.x()-4, end_p.y()+5));
      }
    }

		if (!canvas->isMzToXAxis())
		{
			bounding_box_.setWidth(bounding_box_.width()+10.0);
		}
		
		painter.drawText(bounding_box_, Qt::AlignHCenter, text_);
		
		if (selected_)
		{
			drawBoundingBox_(painter);
		}
	}
	
	void Annotation1DDistanceItem::move(const PointType& delta)
	{
		// shift vertical position according to y-component of delta
		start_point_.setY(start_point_.getY()+delta.getY());
		end_point_.setY(end_point_.getY()+delta.getY());
	}

	
	void Annotation1DDistanceItem::setStartPoint(const PointType& p)
	{
		start_point_ = p;
	}

	void Annotation1DDistanceItem::setEndPoint(const PointType& p)
	{
		end_point_ = p;
	}
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getStartPoint() const
	{
		return start_point_;
	}
	
	const Annotation1DDistanceItem::PointType& Annotation1DDistanceItem::getEndPoint() const
	{
		return end_point_;
	}

  void Annotation1DDistanceItem::setTicks(const std::vector<DoubleReal> &ticks)
  {
    ticks_ = ticks;
  }

	void Annotation1DDistanceItem::ensureWithinDataRange(Spectrum1DCanvas* const canvas)
	{
		// can only be moved vertically, so check only y-position
		DRange<3> data_range = canvas->getDataRange();
		CoordinateType y_pos = start_point_.getY() * canvas->getPercentageFactor();
		
		if (y_pos < data_range.minPosition()[1])
		{
			start_point_.setY(data_range.minPosition()[1] / canvas->getPercentageFactor());
			end_point_.setY(start_point_.getY());
		}
		if (y_pos > data_range.maxPosition()[1])
		{
			start_point_.setY(data_range.maxPosition()[1] / canvas->getPercentageFactor());
			end_point_.setY(start_point_.getY());
		}
	}
	
}//Namespace




