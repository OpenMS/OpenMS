// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/PeakIcon.h>

#include <QtGui/QPainter>

namespace OpenMS
{

	void PeakIcon::drawIcon(Icon icon, QPainter& painter, const QRect& r)
	{
		switch(icon)
			{
				case IT_ELLIPSE:
					drawEllipse(painter, r);
					break;
				case IT_TRIANGLE:
					drawTriangle(painter, r);
					break;
				case IT_ASTERIX:
					drawAsterix(painter, r);
					break;
				case IT_SQUARE:
					drawRectangle(painter, r);
					break;
				default:
					break;
			}
	}
	
	
	void PeakIcon::drawEllipse(QPainter& painter, const QRect& r)
	{
		painter.save();
	  painter.drawEllipse(r);
		painter.restore();
	}
	
	void PeakIcon::drawTriangle(QPainter& painter, const QRect& r)
	{
		painter.save();
		painter.drawLine(r.x(), r.y() + r.height(), r.x() + int(r.width()/2), r.y());
	  painter.drawLine(r.x() + int(r.width()/2), r.y(), r.x() + r.width(), r.y() + r.height());
		painter.drawLine(r.x(), r.y() + r.height(), r.x() + r.width(), r.y() + r.height());
		painter.restore();
	}
	
	void PeakIcon::drawAsterix(QPainter& painter, const QRect& r)
	{
		painter.save();
		painter.drawLine(r.x(), r.y(), r.x() + r.width(), r.y() + r.height());
		painter.drawLine(r.x(), r.y() + r.height(), r.x() + r.width(), r.y());
		painter.drawLine(r.x(), r.y() + int(r.height()/2), r.x() + r.width(), int(r.height()/2));
		painter.drawLine(r.x() + int(r.width()/2), r.y(), r.x() + int(r.width()/2), r.y() + r.height());
		painter.restore();
	}
	
	void PeakIcon::drawRectangle(QPainter& painter, const QRect& r)
	{
		painter.save();
	  painter.drawRect(r);
	  painter.restore();
	}

} //namespace
