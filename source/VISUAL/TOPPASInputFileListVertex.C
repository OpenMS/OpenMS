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

#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>

namespace OpenMS
{
	TOPPASInputFileListVertex::TOPPASInputFileListVertex()
		:	TOPPASVertex(),
			files_(),
			key_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const QStringList& files)
		:	TOPPASVertex(),
			files_(files),
			key_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs)
		:	TOPPASVertex(rhs),
			files_(rhs.files_),
			key_()
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
		key_ = rhs.key_;
		
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
			qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
			emit somethingHasChanged();
		}
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
		QString text = QString::number(files_.size())+" input file"
										+(files_.size() == 1 ? "" : "s");
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);

		//topo sort number
    qreal x_pos = -63.0;
    qreal y_pos = -19.0;
    painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
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

    QString toppview_executable;
    toppview_executable = "TOPPView";

    p->start(toppview_executable, files_);
    if(!p->waitForStarted())
    {
      // execution failed
      std::cerr << p->errorString().toStdString() << std::endl;
#if defined(Q_WS_MAC)
      std::cerr << "Please check if TOPPAS and TOPPView are located in the same directory" << std::endl;
#endif

    }
	}
	
	void TOPPASInputFileListVertex::startPipeline()
	{
		if (files_.empty())
		{
			return;
		}
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
			if (ttv)
			{
				if (!ttv->isAlreadyStarted())
				{
					ttv->runToolIfInputReady();
					ttv->setAlreadyStarted(true);
				}
				continue;
			}
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(tv);
			if (mv)
			{
				if (!mv->isAlreadyStarted())
				{
					mv->forwardPipelineExecution();
					mv->setAlreadyStarted(true);
				}
				continue;
			}
		}
	}
	
	void TOPPASInputFileListVertex::checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run)
	{
		__DEBUG_BEGIN_METHOD__
		
		sc_files_per_round_ = files_.size();
		sc_files_total_ = sc_files_per_round_;
		sc_list_length_checked_ = true;
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			tv->checkListLengths(unequal_per_round, unequal_over_entire_run);
		}
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASInputFileListVertex::setKey(const QString& key)
	{
		key_ = key;
		setToolTip(key);
	}
	
	const QString& TOPPASInputFileListVertex::getKey()
	{
		return key_;
	}
	
	void TOPPASInputFileListVertex::setFilenames(const QStringList& files)
	{
		files_ = files;
	}
	
}
