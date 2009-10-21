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

#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
	TOPPASInputFileListVertex::TOPPASInputFileListVertex()
		:	TOPPASVertex(),
			files_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const QStringList& files)
		:	TOPPASVertex(),
			files_(files)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs)
		:	TOPPASVertex(rhs),
			files_(rhs.files_)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::~TOPPASInputFileListVertex()
	{
	
	}
	
	TOPPASInputFileListVertex& TOPPASInputFileListVertex::operator= (const TOPPASInputFileListVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		
		files_ = rhs.files_;
		
		return *this;
	}
	
	void TOPPASInputFileListVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		showFilesDialog();
	}
	
	void TOPPASInputFileListVertex::showFilesDialog()
	{
		TOPPASInputFilesDialog tifd(files_);
		if (tifd.exec())
		{
			tifd.getFilenames(files_);
			qobject_cast<TOPPASScene*>(scene())->setChanged(true);
		}
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
		
		emit somethingHasChanged();
	}
	
	const QStringList& TOPPASInputFileListVertex::getFilenames()
	{
		return files_;
	}
	
	void TOPPASInputFileListVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
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
		QString text = "Input files";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
	}
	
	QRectF TOPPASInputFileListVertex::boundingRect() const
	{
		return QRectF(-71,-41,142,82);
	}
	
	QPainterPath TOPPASInputFileListVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-71.0, -41.0, 142.0, 81.0, 20, 20);
		return shape;
	}
	
	bool TOPPASInputFileListVertex::fileNamesValid(const QStringList& /*files*/)
	{
		// some more checks TODO..
		return true;
	}
	
	void TOPPASInputFileListVertex::openInTOPPView()
	{
		QProcess* p = new QProcess();
		p->setProcessChannelMode(QProcess::ForwardedChannels);
		p->start("TOPPView", files_);
	}
	
}
