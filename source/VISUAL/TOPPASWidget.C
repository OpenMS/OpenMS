// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

// OpenMS
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>


// Qt
#include <QtGui/QDragEnterEvent>
#include <QtGui/QDragMoveEvent>
#include <QtGui/QDropEvent>
#include <QtCore/QMimeData>

using namespace std;

namespace OpenMS
{
	TOPPASWidget::TOPPASWidget(const Param& /*preferences*/, QWidget* parent)
		:	QGraphicsView(parent),
			scene_(new TOPPASScene())
	{
		setAttribute(Qt::WA_DeleteOnClose);
		setRenderHint(QPainter::Antialiasing);
		setScene(scene_);
		setAcceptDrops(true);
		setDragMode(QGraphicsView::RubberBandDrag);
	}
	
	TOPPASWidget::~TOPPASWidget()
	{
		emit aboutToBeDestroyed(window_id);
	}
	
	TOPPASScene* TOPPASWidget::getScene()
	{
		return scene_;
	}
	
	void TOPPASWidget::dragEnterEvent(QDragEnterEvent* event)
	{
		// TODO: test mime type
		event->accept();
	}
	
	void TOPPASWidget::dragMoveEvent(QDragMoveEvent* event)
	{
		// TODO: test mime type
		event->accept();
	}
	
	void TOPPASWidget::dropEvent(QDropEvent* event)
	{
		QPointF scene_pos = mapToScene(event->pos());
		emit toolDroppedOnWidget(scene_pos.x(), scene_pos.y());
	}
	
} //Namespace

