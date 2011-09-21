// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QDesktopServices>
#include <QUrl>
#include <QMessageBox>
#include <QSvgRenderer>

namespace OpenMS
{
	TOPPASInputFileListVertex::TOPPASInputFileListVertex()
		:	TOPPASVertex(),
			key_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const QStringList& files)
		:	TOPPASVertex(),
			key_()
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
    setFilenames(files);
	}
	
	TOPPASInputFileListVertex::TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs)
		:	TOPPASVertex(rhs),
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
		
		key_ = rhs.key_;
		
		return *this;
	}
	
	void TOPPASInputFileListVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		showFilesDialog();
	}
	
	void TOPPASInputFileListVertex::showFilesDialog()
	{
    TOPPASInputFilesDialog tifd(this->getFileNames());
		if (tifd.exec())
		{
      QStringList updated_filelist;
			tifd.getFilenames(updated_filelist);
      this->setFilenames(updated_filelist); // to correct filenames (separators etc)
			qobject_cast<TOPPASScene*>(scene())->setChanged(true);
			qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
			emit somethingHasChanged();
		}
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
    QString text = QString::number(getFileNames().size())
                   + " input file"
									 + (getFileNames().size() == 1 ? "" : "s");
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);

		//topo sort number
    qreal x_pos = -63.0;
    qreal y_pos = -19.0;
    painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
	
		// recycling status
    if (this->allow_output_recycling_)
    {
      painter->setPen(Qt::green);
      QSvgRenderer* svg_renderer = new QSvgRenderer(QString(":/Recycling_symbol.svg"), 0);
      svg_renderer->render(painter, QRectF(-7, -32, 14, 14));
    }
		
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
	
	bool TOPPASInputFileListVertex::fileNamesValid()
	{
    QStringList fl = getFileNames();
    foreach (const QString& file, fl)
		{
      if (!File::exists(file))
      {
				return false;
      }
		}
    return true;
	}
	
	void TOPPASInputFileListVertex::openContainingFolder()
	{
    std::set<String> directories;
    QStringList fl = getFileNames();
    for (int i = 0; i < fl.size(); ++i)
    { // collect unique directories
      QFileInfo fi(fl[i]);
      directories.insert(String(QFileInfo(fi.canonicalFilePath()).path()));
    }

    // open them
    for (std::set<String>::const_iterator it=directories.begin();it!=directories.end();++it)
    {
      QString path = QDir::toNativeSeparators(it->toQString());
      if (!QDir(path).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + path, QUrl::TolerantMode))))
      {
        QMessageBox::warning(0, "Open Folder Error", String("The folder " + path + " could not be opened!").toQString());
      }
    }
	}
	
	void TOPPASInputFileListVertex::run()
	{
    round_total_   = (int) output_files_.size(); // for now each file is one round; for the future we might allow to create blocks of files (e.g. for replicate measurements)
    round_counter_ = (int) round_total_;

    this->finished_ = true; // input node is ready to go (file check was already done)

    //std::cerr << "#" << this->getTopoNr() << " set #rounds: " << round_total_ << "\n";

		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			if (tv && !tv->isFinished()) // this tool might have already been called by another path, so do not call it again (as this will throw an error)
			{
				tv->run();
			}
		}
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
    output_files_.clear();
    output_files_.resize(files.size()); // for now, assume one file per round (we could later extend that)
    for (int f=0; f<files.size(); ++f)
    {
      output_files_[f][-1].filenames << QDir::toNativeSeparators(files[f]);
    }
	}
	
}
