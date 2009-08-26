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
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>

namespace OpenMS
{
	TOPPASOutputFileVertex::TOPPASOutputFileVertex()
		:	TOPPASVertex()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileVertex::TOPPASOutputFileVertex(const TOPPASOutputFileVertex& rhs)
		:	TOPPASVertex(rhs)
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
		
		return *this;
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
		
		//topo sort number
		qreal x_pos = -62.0;
		qreal y_pos = 28.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
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
	
	void TOPPASOutputFileVertex::finished()
	{
		// copy tmp files to output dir
		TOPPASEdge* e = *inEdgesBegin();
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
		const QVector<QStringList>& output_files = tv->getOutputFileNames();
		int param_index = e->getSourceOutParam();
		QStringList tmp_file_names = output_files[param_index];
		QString parent_dir = qobject_cast<TOPPASScene*>(scene())->getOutDir();
		
		file_ = "";
		if (!tmp_file_names.isEmpty())
		{
			QString old_file = tmp_file_names.first();
			if (!File::exists(old_file))
			{
				std::cerr << "The file '" << String(old_file) << "' does not exist!" << std::endl;
				return;
			}
			QString new_file = parent_dir + QDir::separator() + getOutputDir().toQString()+QDir::separator()+File::basename(old_file).toQString();
			if (new_file.endsWith("_tmp"))
			{
				new_file.truncate(new_file.size() - 4);
			}
			new_file += "_processed";
			// get file type and rename
			FileTypes::Type ft = FileHandler::getTypeByContent(old_file);
			if (ft != FileTypes::UNKNOWN)
			{
				new_file += "."+(FileHandler::typeToName(ft)).toQString();
			}
			if (File::exists(new_file))
			{
				QFile::remove(new_file);
			}
			if (!QFile::copy(old_file, new_file))
			{
				std::cerr << "Could not copy tmp output file " << String(old_file) << " to " << String(new_file) << std::endl;
			}
			else
			{
				file_ = new_file;
				emit outputFileWritten(new_file);
			}
		}
		
		finished_ = true;
		emit iAmDone();
	}
	
	void TOPPASOutputFileVertex::inEdgeHasChanged()
	{
		file_ = "";
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
		TOPPASVertex::inEdgeHasChanged();
	}
	
	void TOPPASOutputFileVertex::contextMenuEvent(QGraphicsSceneContextMenuEvent* event)
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		ts->unselectAll();
		setSelected(true);
		
		QMenu menu;
		
		QAction* open_action = menu.addAction("Open file in TOPPView");
		open_action->setEnabled(false);
		if (file_ != "")
		{
			open_action->setEnabled(true);
		}
		
		menu.addAction("Remove");
		
		QAction* selected_action = menu.exec(event->screenPos());
		if (selected_action)
		{
			QString text = selected_action->text();
			if (text == "Remove")
			{
				ts->removeSelected();
			}
			else if (text == "Open file in TOPPView")
			{
				QProcess* p = new QProcess();
				p->setProcessChannelMode(QProcess::ForwardedChannels);
				p->start("TOPPView", QStringList(file_));
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
	
	String TOPPASOutputFileVertex::getOutputDir()
	{
		String dir = String("TOPPAS_out")+String(QDir::separator())+get3CharsNumber_(topo_nr_);
		
		return dir;
	}
	
	void TOPPASOutputFileVertex::createDirs(const QString& out_dir)
	{
		QDir current_dir(out_dir);
		String new_dir = getOutputDir();
		if (!File::exists(new_dir))
		{
			if (!current_dir.mkpath(new_dir.toQString()))
			{
				std::cerr << "Could not create path " << new_dir << std::endl;
			}
		}
	}
	
	void TOPPASOutputFileVertex::setTopoNr(UInt nr)
	{
		if (topo_nr_ != nr)
		{
			topo_nr_ = nr;
		}
	}
}

