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
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>

#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>

namespace OpenMS
{
	TOPPASOutputFileListVertex::TOPPASOutputFileListVertex()
		:	TOPPASVertex()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASOutputFileListVertex::TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex& rhs)
		:	TOPPASVertex(rhs)
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
		
		return *this;
	}
	
	void TOPPASOutputFileListVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
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
		QString text = "Output files";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
		
		//topo sort number
		qreal x_pos = -62.0;
		qreal y_pos = 28.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
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
	
	void TOPPASOutputFileListVertex::finish()
	{
		// copy tmp files to output dir
		
		TOPPASEdge* e = *inEdgesBegin();
		TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
		const QVector<QStringList>& output_files = tv->getCurrentOutputFileNames();
		int param_index = e->getSourceOutParam();
		const QStringList& tmp_file_names = output_files[param_index];
		QString parent_dir = qobject_cast<TOPPASScene*>(scene())->getOutDir();
		
		createDirs();

		if (!tmp_file_names.isEmpty())
		{
			foreach (const QString& f, tmp_file_names)
			{
				if (!File::exists(f))
				{
					std::cerr << "The file '" << String(f) << "' does not exist!" << std::endl;
					continue;
				}
				QString new_file = parent_dir + QDir::separator() + getOutputDir().toQString()+QDir::separator()+File::basename(f).toQString();
				if (new_file.endsWith("_tmp"))
				{
					new_file.truncate(new_file.size() - 4);
				}
				new_file += "_processed";
				// get file type and rename
				FileTypes::Type ft = FileHandler::getTypeByContent(f);
				if (ft != FileTypes::UNKNOWN)
				{
					new_file += "."+(FileHandler::typeToName(ft)).toQString();
				}
				if (File::exists(new_file))
				{
					QFile::remove(new_file);
				}
				if (!QFile::copy(f, new_file))
				{
					std::cerr << "Could not copy tmp output file " << String(f) << " to " << String(new_file) << std::endl;
				}
				else
				{
					files_.push_back(new_file);
					emit outputFileWritten(new_file);
				}
			}
		}
		
		finished_ = true;
		checkIfSubtreeFinished();
		emit iAmDone();
	}
	
	void TOPPASOutputFileListVertex::inEdgeHasChanged()
	{
		reset(true);
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
		TOPPASVertex::inEdgeHasChanged();
	}
	
	void TOPPASOutputFileListVertex::openInTOPPView()
	{
		QProcess* p = new QProcess();
		p->setProcessChannelMode(QProcess::ForwardedChannels);
		p->start("TOPPView", files_);
	}
	
	bool TOPPASOutputFileListVertex::isFinished()
	{
		return finished_;
	}
	
	String TOPPASOutputFileListVertex::getOutputDir()
	{
		String dir = String("TOPPAS_out")+String(QDir::separator())+get3CharsNumber_(topo_nr_);
		
		return dir;
	}
	
	void TOPPASOutputFileListVertex::createDirs()
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		QString out_dir = ts->getOutDir();
		QDir current_dir(out_dir);
		String new_dir = getOutputDir();
		if (!File::exists(String(out_dir)+String(QDir::separator())+new_dir))
		{
			if (!current_dir.mkpath(new_dir.toQString()))
			{
				std::cerr << "Could not create path " << new_dir << std::endl;
			}
		}
	}
	
	void TOPPASOutputFileListVertex::setTopoNr(UInt nr)
	{
		if (topo_nr_ != nr)
		{
			// topological number changes --> output dir changes --> reset
			reset(true);
			topo_nr_ = nr;
		}
	}
	
	void TOPPASOutputFileListVertex::reset(bool reset_all_files)
	{
		TOPPASVertex::reset();
		finished_ = false;
		
		if (reset_all_files)
		{
			files_.clear();
			// do not actually delete the output files here
		}
	}
}

