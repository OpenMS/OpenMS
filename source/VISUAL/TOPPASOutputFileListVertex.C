// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
#include <QDesktopServices>
#include <QUrl>
#include <QMessageBox>

#include <QCoreApplication>

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
    QString text =  QString::number(files_written_) + "/" 
                  + QString::number(files_total_) + " output file" + (files_total_ == 1 ? "" : "s");
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
		
		//topo sort number
		qreal x_pos = -63.0;
		qreal y_pos = -19.0; 
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
	
	void TOPPASOutputFileListVertex::run()
	{
		__DEBUG_BEGIN_METHOD__
		
		// copy tmp files to output dir
		
    // get incoming edge and preceding vertex
		TOPPASEdge* e = *inEdgesBegin();
		TOPPASVertex* tv = e->getSourceVertex();
    RoundPackages pkg = tv->getOutputFiles();
    if (pkg.empty())
    {
			std::cerr << "A problem occured while grabbing files from the parent tool. This is a bug, please report it!" << std::endl;
			__DEBUG_END_METHOD__
			return;
		}

		String full_dir = createOutputDir(); // create output dir

    round_total_ = (int) pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0;        // once round_counter_ reaches round_total_, we are done

    // clear output file list
		output_files_.clear();
    output_files_.resize(pkg.size()); // #rounds

    files_total_ = 0;
    files_written_ = 0;

    bool dry_run = qobject_cast<TOPPASScene*>(scene())->isDryRun();

    int param_index_src = e->getSourceOutParam();
    int param_index_me = e->getTargetInParam();
    for (Size round=0; round < pkg.size(); ++round)
    {
      foreach (const QString& f, pkg[round][param_index_src].filenames)
		  {
			  if (!dry_run && !File::exists(f))
			  {
				  std::cerr << "The file '" << String(f) << "' does not exist!" << std::endl;
				  continue;
			  }
        QString new_file = full_dir.toQString()
                           + QDir::separator()
                           + File::basename(f).toQString().left(190); // ensure 190 char filename (we might append more and reach ~NTFS limit)
			
        // remove "_tmp<number>" if its a suffix
        QRegExp rx("_tmp\\d+$");
			  int tmp_index = rx.indexIn(new_file);
        if (tmp_index != -1) new_file = new_file.left(tmp_index);
			
        // get file type and rename
			  FileTypes::Type ft = (dry_run ? FileTypes::UNKNOWN : FileHandler::getTypeByContent(f)); // this will access the file physically
			  if (ft != FileTypes::UNKNOWN)
			  {
          QString suffix = QString(".") + FileHandler::typeToName(ft).toQString();
          if (!new_file.endsWith(suffix)) new_file += suffix;
			  }

        // only scheduled for writing
        output_files_[round][param_index_me].filenames.push_back(QDir::toNativeSeparators(new_file));
        ++files_total_; 
		  }
    }

    // do the actual copying
    if (dry_run) // assume the copy worked
    {
      files_written_ = files_total_;
      update(boundingRect()); // repaint
    }
    else
    {
      for (Size round=0; round < pkg.size(); ++round)
      {
        round_counter_ = (int)round; // for global update, in case someone asks
        for (int i = 0; i < pkg[round][param_index_src].filenames.size(); ++i)
		    {
          QString file_from = pkg[round][param_index_src].filenames[i];
          QString file_to   = output_files_[round][param_index_me].filenames[i];
			    if (File::exists(file_to))
			    {
				    if (!QFile::remove(file_to))
            {
              std::cerr << "Could not delete old output file " << String(file_to) << " before overwriting with new one."<< std::endl;
            }
			    }
			    if (!QFile::copy(file_from, file_to))
			    {
				    std::cerr << "Could not copy tmp output file " << String(file_from) << " to " << String(file_to) << std::endl;
			    }
			    else
			    {
            ++files_written_;
				    emit outputFileWritten(file_to);
			    }
        }
        update(boundingRect()); // repaint
      }
    }

		finished_ = true;
		emit iAmDone();
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASOutputFileListVertex::inEdgeHasChanged()
	{
		reset(true);
		qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
		TOPPASVertex::inEdgeHasChanged();
	}

	void TOPPASOutputFileListVertex::openContainingFolder()
	{
    QString path = getFullOutputDirectory().toQString();
    if (!QDir(path).exists() || (!QDesktopServices::openUrl(QUrl("file:///" + path, QUrl::TolerantMode))))
    {
      QMessageBox::warning(0, "Open Folder Error", "The folder " + path + " could not be opened!");
    }
	}

  String TOPPASOutputFileListVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    String dir = String(ts->getOutDir()).substitute("\\", "/");
    return QDir::cleanPath((dir.ensureLastChar('/') + getOutputDir()).toQString() );
  }

	String TOPPASOutputFileListVertex::getOutputDir() const
	{
    String dir = String("TOPPAS_out") + String(QDir::separator()) + get3CharsNumber_(topo_nr_);
    return dir;
	}
	
	String TOPPASOutputFileListVertex::createOutputDir()
	{
    String full_dir = getFullOutputDirectory();
		if (!File::exists(full_dir))
		{
			QDir dir;
      if (!dir.mkpath(full_dir.toQString()))
			{
				std::cerr << "Could not create path " << full_dir << std::endl;
			}
		}
    return full_dir;
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
		__DEBUG_BEGIN_METHOD__
		
		files_total_ = 0;
    files_written_ = 0;

    TOPPASVertex::reset();
		
		if (reset_all_files)
		{
			// do not actually delete the output files here
		}
		
		__DEBUG_END_METHOD__
	}
}

