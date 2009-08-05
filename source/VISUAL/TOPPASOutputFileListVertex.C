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

#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QMessageBox>
#include <QtCore/QFile>


namespace OpenMS
{
	TOPPASOutputFileListVertex::TOPPASOutputFileListVertex()
		:	TOPPASVertex(),
			files_(),
			ready_(false)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileListVertex::TOPPASOutputFileListVertex(const QStringList& files)
		:	TOPPASVertex(),
			files_(files),
			ready_(false)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileListVertex::TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex& rhs)
		:	TOPPASVertex(rhs),
			files_(rhs.files_),
			ready_(rhs.ready_)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileListVertex::~TOPPASOutputFileListVertex()
	{
	
	}
	
	TOPPASOutputFileListVertex& TOPPASOutputFileListVertex::operator= (const TOPPASOutputFileListVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);		
		
		files_ = rhs.files_;
		ready_ = rhs.ready_;
		
		return *this;
	}
	
	void TOPPASOutputFileListVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		bool go = true;
		if (inEdgesBegin() == inEdgesEnd())
		{
			go = false;
		}
		else
		{
			TOPPASToolVertex* parent_tv = qobject_cast<TOPPASToolVertex*>((*inEdgesBegin())->getSourceVertex());
			
			//first, update file names recursively
			parent_tv->updateOutputFileNames();
			
			const QVector<QStringList>& output_files = parent_tv->getOutputFileNames();
			int param_index = (*inEdgesBegin())->getSourceOutParam();
			if (param_index == -1 || output_files[param_index].empty())
			{
				go = false;
			}
		}
		
		if (!go)
		{
			QMessageBox::information(0,"No input","The output file names cannot be specified until the number of files is known. You have to set up the rest of the pipeline, first (including input files).");
			return;
		}
		
		TOPPASOutputFilesDialog tofd(this);
		if (tofd.exec())
		{
			tofd.getFilenames(files_);
		}
		updateStatus();
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
	}
	
	void TOPPASOutputFileListVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
	
		QPainterPath path;
		path.addRoundRect(-70.0, -40.0, 140.0, 80.0, 20, 20);		
 		painter->drawPath(path);
 		
		QString text = "Output file list";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
	}
	
	QRectF TOPPASOutputFileListVertex::boundingRect() const
	{
		return QRectF(-71,-41,142,82);
	}
	
	QPainterPath TOPPASOutputFileListVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-71.0, -41.0, 142.0, 81.0, 20, 20);
		return shape;
	}
	
	void TOPPASOutputFileListVertex::startComputation()
	{
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
			if (ttv)
			{
				ttv->runRecursively();
			}
		}
	}
	
	const QStringList& TOPPASOutputFileListVertex::getFilenames()
	{
		return files_;
	}
	
	void TOPPASOutputFileListVertex::finished()
	{
		// rename tmp out file if file names specified
		TOPPASEdge* e = *inEdgesBegin();
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
		const QVector<QStringList>& output_files = tv->getOutputFileNames();
		int param_index = e->getSourceOutParam();
		const QStringList& tmp_file_names = output_files[param_index];
		
		int specified_names_count = files_.size();
		int counter = 0;
		foreach (const QString& file, tmp_file_names)
		{
			if (counter < specified_names_count)
			{
				const QString& save_name = files_[counter++];	
				if (QFile::exists(save_name))
				{
					QFile::remove(save_name);
				}
				QFile::rename(file, save_name);
				emit outputFileWritten(String(save_name));
			}
		}
	}
	
	void TOPPASOutputFileListVertex::inEdgeHasChanged()
	{
		updateStatus();
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
		
		// we do not need to forward the change (we have no childs)
	}
	
	void TOPPASOutputFileListVertex::updateStatus()
	{
		if (inEdgesBegin() == inEdgesEnd())
		{
			ready_ = false;
			return;
		}
		TOPPASEdge* e = *inEdgesBegin();
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
		tv->updateOutputFileNames();
		const QVector<QStringList>& output_files = tv->getOutputFileNames();
		int param_index = e->getSourceOutParam();
		if (param_index == -1)
		{
			ready_ = false;
			return;
		}
		const QStringList& tmp_file_names = output_files[param_index];
		int tmp_files_count = tmp_file_names.size();
		if (tmp_files_count == 0)
		{
			ready_ = false;
			return;
		}
		int specified_names_count = files_.size();
		if (specified_names_count > tmp_files_count)
		{
			// truncate superfluous file names
			int diff = specified_names_count - tmp_files_count;
			for (int i = 0; i < diff; ++i)
			{
				files_.removeLast();
			}
			ready_ = true;
		}
		else
		{
			ready_ = (tmp_files_count == specified_names_count);
		}
	}
	
	bool TOPPASOutputFileListVertex::fileNamesValid(const QStringList& /*files*/)
	{
		// some more checks TODO...
		return true;
	}
	
	bool TOPPASOutputFileListVertex::isReady()
	{
		return ready_;
	}
}

