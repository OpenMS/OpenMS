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

#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFileDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QFile>

namespace OpenMS
{
	TOPPASOutputFileVertex::TOPPASOutputFileVertex()
		:	TOPPASVertex(),
			file_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::TOPPASOutputFileVertex(const QString& file)
		:	TOPPASVertex(),
			file_(file)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::TOPPASOutputFileVertex(const TOPPASOutputFileVertex& rhs)
		:	TOPPASVertex(rhs),
			file_(rhs.file_)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::~TOPPASOutputFileVertex()
	{
	
	}
	
	TOPPASOutputFileVertex& TOPPASOutputFileVertex::operator= (const TOPPASOutputFileVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);		
		
		file_ = rhs.file_;
		
		return *this;
	}
	
	void TOPPASOutputFileVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		showFileDialog();
	}
	
	void TOPPASOutputFileVertex::showFileDialog()
	{
		TOPPASOutputFileDialog tofd(file_);
		if (tofd.exec())
		{
			file_ = tofd.getFilename();
		}
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
	}
	
	void TOPPASOutputFileVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
			pen.setColor(Qt::darkBlue);
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
	
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
 		painter->drawPath(path);
 		
 		pen.setColor(pen_color_);
 		painter->setPen(pen);
		QString text = "Output file";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
	}
	
	QRectF TOPPASOutputFileVertex::boundingRect() const
	{
		return QRectF(-71,-41,142,82);
	}
	
	QPainterPath TOPPASOutputFileVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-71.0, -41.0, 142.0, 81.0, 20, 20);
		return shape;
	}
	
	void TOPPASOutputFileVertex::startComputation()
	{
		finished_ = false;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (ttv)
			{
				ttv->runRecursively();
			}
		}
	}
	
	const QString& TOPPASOutputFileVertex::getFilename()
	{
		return file_;
	}
	
	void TOPPASOutputFileVertex::finished()
	{
		// rename tmp out file if a file name was specified
		if (file_ != "")
		{
			TOPPASEdge* e = *inEdgesBegin();
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
			const QVector<QStringList>& output_files = tv->getOutputFileNames();
			int param_index = e->getSourceOutParam();
			QString tmp_file_name = output_files[param_index].first();
			
			if (QFile::exists(file_))
			{
				QFile::remove(file_);
			}
			QFile::rename(tmp_file_name, file_);
			
			emit outputFileWritten(String(file_));
		}
		finished_ = true;
		emit iAmDone();
	}
	
	void TOPPASOutputFileVertex::inEdgeHasChanged()
	{
		// we do not need to forward the change (we have no childs)
	}
	
	void TOPPASOutputFileVertex::contextMenuEvent(QGraphicsSceneContextMenuEvent* event)
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		ts->unselectAll();
		setSelected(true);
		
		QMenu menu;
		menu.addAction("Change file");
		menu.addAction("Remove");
		
		QAction* selected_action = menu.exec(event->screenPos());
		if (selected_action)
		{
			QString text = selected_action->text();
			if (text == "Change file")
			{
				showFileDialog();
			}
			else if (text == "Remove")
			{
				ts->removeSelected();
			}
			event->accept();
		}
		else
		{
			event->ignore();	
		}
	}
	
	bool TOPPASOutputFileVertex::isFinished()
	{
		return finished_;
	}
}

