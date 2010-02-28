// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

// OpenMS
#include <OpenMS/CONCEPT/Types.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>


// Qt
#include <QtGui/QDragEnterEvent>
#include <QtGui/QDragMoveEvent>
#include <QtGui/QDropEvent>
#include <QtCore/QMimeData>

using namespace std;

namespace OpenMS
{
	TOPPASWidget::TOPPASWidget(const Param& /*preferences*/, QWidget* parent, const String& tmp_path)
		:	QGraphicsView(parent),
			scene_(new TOPPASScene(this, tmp_path))
	{
		setAttribute(Qt::WA_DeleteOnClose);
		setRenderHint(QPainter::Antialiasing);
		setScene(scene_);
		setAcceptDrops(true);
		setDragMode(QGraphicsView::RubberBandDrag);
		setFocusPolicy(Qt::StrongFocus);
	}
	
	TOPPASWidget::~TOPPASWidget()
	{
		emit aboutToBeDestroyed(window_id);
	}
	
	TOPPASScene* TOPPASWidget::getScene()
	{
		return scene_;
	}
	
	void TOPPASWidget::zoom(bool zoom_in)
	{
		qreal factor = 1.2;
		if (zoom_in)
		{
			factor = 1.0 / factor;
		}
		scale(factor, factor);
		
		QRectF items_rect = scene_->itemsBoundingRect();
		scene_->setSceneRect(items_rect.united(mapToScene(rect()).boundingRect()));
	}
	
	void TOPPASWidget::wheelEvent(QWheelEvent* event)
	{
		zoom(event->delta() < 0);
	}

	
	void TOPPASWidget::dragEnterEvent(QDragEnterEvent* event)
	{
		// TODO: test mime type/source? where?
		event->acceptProposedAction();
	}
	
	void TOPPASWidget::dragMoveEvent(QDragMoveEvent* event)
	{
		// TODO: test mime type/source? where?
		event->acceptProposedAction();
	}
	
	void TOPPASWidget::dropEvent(QDropEvent* event)
	{
		// TODO: test mime type/source? where?
		QPointF scene_pos = mapToScene(event->pos());
		emit toolDroppedOnWidget(scene_pos.x(), scene_pos.y());
		event->acceptProposedAction();
	}
	
	void TOPPASWidget::keyPressEvent(QKeyEvent* e)
	{
		if (e->key() == Qt::Key_C && e->modifiers() == Qt::ControlModifier)
		{
			scene_->copySelected();
			e->accept();
		}
		else if (e->key() == Qt::Key_X && e->modifiers() == Qt::ControlModifier)
		{
			scene_->copySelected();
			scene_->removeSelected();
			e->accept();
		}
		else if (e->key() == Qt::Key_V && e->modifiers() == Qt::ControlModifier)
		{
			scene_->paste();
			e->accept();
		}
		else if (e->key() == Qt::Key_Control)
		{
			setDragMode(QGraphicsView::ScrollHandDrag);
			e->accept();
		}
		else if (e->key() == Qt::Key_Delete || e->key() == Qt::Key_Backspace)
		{
			scene_->removeSelected();
			e->accept();
		}
		else if (e->key() == Qt::Key_F5)
		{
			scene_->runPipeline();
			e->accept();
		}
		else if (e->key() == Qt::Key_Plus)
		{
			zoom(false);
			e->accept();
		}
		else if (e->key() == Qt::Key_Minus)
		{
			zoom(true);
			e->accept();
		}
	}
	
	void TOPPASWidget::keyReleaseEvent(QKeyEvent* e)
	{
		if (e->key() == Qt::Key_Control)
		{
			setDragMode(QGraphicsView::RubberBandDrag);
			e->accept();
		}
	}
	
	void TOPPASWidget::leaveEvent(QEvent* /*e*/)
	{
		
	}

	void TOPPASWidget::enterEvent(QEvent* /*e*/)
	{
		setFocus();
	}
	
	void TOPPASWidget::resizeEvent(QResizeEvent* event)
	{
		QGraphicsView::resizeEvent(event);
		if (scene_)
		{
			QRectF items_rect = scene_->itemsBoundingRect();
			scene_->setSceneRect(items_rect.united(mapToScene(viewport()->rect()).boundingRect()));
		}
	}
	
	void TOPPASWidget::closeEvent(QCloseEvent* e)
	{
		bool close = scene_->saveIfChanged();
		if (close)
		{
			e->accept();
		}
		else
		{
			e->ignore();
		}
	}

} //Namespace

