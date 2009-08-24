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
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASTreeView.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <QtGui/QDrag>
#include <QtGui/QApplication>
#include <QtCore/QMimeData>

#include <iostream>

using namespace std;

namespace OpenMS
{

	TOPPASTreeView::TOPPASTreeView(QWidget* parent)
		: QTreeWidget(parent)
	{
		// we drag by ourselves:
		setDragEnabled(false);
	}
	
	TOPPASTreeView::~TOPPASTreeView()
	{
		
	}
	
	void TOPPASTreeView::mousePressEvent(QMouseEvent* event)
	{
		QTreeWidget::mousePressEvent(event);
		
		if (event->button() == Qt::LeftButton)
		{
			drag_start_pos_ = event->pos();
		}
	}
	
	void TOPPASTreeView::mouseMoveEvent(QMouseEvent* event)
	{
		QTreeWidget::mouseMoveEvent(event);
		
		if (!(event->buttons() & Qt::LeftButton))
		{
			return;
		}
		if ((event->pos() - drag_start_pos_).manhattanLength() < QApplication::startDragDistance())
		{
			return;
		}
		if (currentItem() && currentItem()->childCount() > 0)
		{
			// drag item is a tool with types - one of the types must be selected
			return;
		}
	
		QDrag* drag = new QDrag(this);
		QMimeData* mime_data = new QMimeData;
		
		mime_data->setText("currently_unused_mime_data");
		drag->setMimeData(mime_data);
		
		// start drag
		drag->exec(Qt::CopyAction);
	}

} //namespace OpenMS	

