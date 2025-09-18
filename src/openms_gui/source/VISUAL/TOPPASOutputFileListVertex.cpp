// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtCore>
#include <QtCore/QFile>
#include <QtCore/QDir>
#include <QtWidgets/QMessageBox>

#include <QCoreApplication>

#include <future>

namespace OpenMS
{
  std::unique_ptr<TOPPASVertex> TOPPASOutputFileListVertex::clone() const
  {
    return std::make_unique<TOPPASOutputFileListVertex>(*this);
  }

  String TOPPASOutputFileListVertex::getName() const
    {
      return "OutputFileVertex";
    }

  void TOPPASOutputFileListVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget);

    QString text = QString::number(files_written_) + "/"
                   + QString::number(files_total_) + " output file" + (files_total_ == 1 ? "" : "s");
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    // display file type(s)
    QStringList text_l = TOPPASVertex::TOPPASFilenames(getFileNames()).getSuffixCounts();
    text = text_l.join(" | ");
    // might get very long, especially if node was not reached yet, so trim
    text = text.left(15) + " ...";
    text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), 35 - (int)(text_boundings.height() / 4.0), text);

    // output folder name
    painter->drawText(painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, output_folder_name_).width()/-2, -25, output_folder_name_);
  }
  
  QRectF TOPPASOutputFileListVertex::boundingRect() const
  {
    return QRectF(-71, -41, 142, 82);
  }

  void TOPPASOutputFileListVertex::run()
  {
    __DEBUG_BEGIN_METHOD__

    // copy tmp files to output dir

    // get incoming edge and preceding vertex
    TOPPASEdge* e = *inEdgesBegin();
    RoundPackages pkg = e->getSourceVertex()->getOutputFiles();
    if (pkg.empty())
    {
      std::cerr << "A problem occurred while grabbing files from the parent tool. This is a bug, please report it!" << std::endl;
      __DEBUG_END_METHOD__
      return;
    }

    String full_dir = createOutputDir(); // create output dir

    round_total_ = (int)pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0;             // once round_counter_ reaches round_total_, we are done

    // clear output file list
    output_files_.clear();
    output_files_.resize(pkg.size()); // #rounds

    files_total_ = 0;
    files_written_ = 0;

    const TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    bool dry_run = ts->isDryRun();

    int param_index_src = e->getSourceOutParam();
    int param_index_me = e->getTargetInParam();
    for (Size round = 0; round < pkg.size(); ++round)
    {
      foreach (const QString& f, pkg[round][param_index_src].filenames.get())
      {
        if (! dry_run && ! File::exists(f))
        {
          OPENMS_LOG_ERROR << "The file '" << String(f) << "' does not exist!" << std::endl;
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, f.toStdString());
        }
        QString new_file = full_dir.toQString() + QDir::separator() + File::basename(f).toQString();

        // remove "_tmp<number>" if its a suffix
        QRegularExpression rx("_tmp\\d+$");
        int tmp_index = new_file.indexOf(rx);
        if (tmp_index != -1)
        {
          new_file = new_file.left(tmp_index);
        }

        // get file type and rename
        FileTypes::Type ft = FileTypes::UNKNOWN;
        if (! dry_run)
        {
          TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
          if (ttv)
          {
            QVector<TOPPASToolVertex::IOInfo> source_output_files = ttv->getOutputParameters();
            if (e->getSourceOutParam() < source_output_files.size())
            {
              StringList types = source_output_files[e->getSourceOutParam()].valid_types;
              if (types.size() == 1) // if suffix is unambiguous
              {
                ft = FileTypes::nameToType(types[0]);
              }
              else
              {
                ft = FileHandler::getTypeByContent(f); // this will access the file physically
              }
              // do we know the extension already?
              if (ft == FileTypes::UNKNOWN)
              { // if not, try param value of '<name>_type' (more generic, supporting more than just 'out_type')
                const Param& p = ttv->getParam();
                String out_type = source_output_files[e->getSourceOutParam()].param_name + "_type";
                if (p.exists(out_type)) { ft = FileTypes::nameToType(p.getValue(out_type).toString()); }
              }
            }
          }
        }

        // replace old suffix by new suffix
        FileHandler::swapExtension(new_file, ft);

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
      for (Size round = 0; round < pkg.size(); ++round)
      {
        round_counter_ = (int)round; // for global update, in case someone asks
        for (int i = 0; i < pkg[round][param_index_src].filenames.size(); ++i)
        {
          String file_from = pkg[round][param_index_src].filenames[i];
          String file_to = output_files_[round][param_index_me].filenames[i];
          if (File::exists(file_to))
          {
            if (! QFile::remove(file_to.toQString())) // todo: this goes wrong on first run .... why???
            {
              String msg = "Error: Could not remove old output file '" + file_to + "' for node '"
                           + pkg[round][param_index_src].edge->getTargetVertex()->getName()
                           + "' in preparation to write the new one. Please make sure the file is not open in other applications and try again.";
              OPENMS_LOG_ERROR << msg << std::endl;
              if (ts->isGUIMode()) { QMessageBox::warning(nullptr, tr("File removing failed"), tr(msg.c_str()), QMessageBox::Ok); }
              else
              {
                throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
              }
            }
          }

          // running the copy() in an extra thread, such that the GUI stays responsive

          auto copyFuture = std::async(std::launch::async, &File::copy, file_from, file_to);
          while (copyFuture.wait_for(std::chrono::milliseconds(25)) != std::future_status::ready)
          {
            qApp->processEvents(); // GUI responsiveness
          }

          if (copyFuture.get())
          {
            ++files_written_;
            emit outputFileWritten(file_to);
          }
          else
          {
            String msg = "Error: Could not copy temporary output file '" + file_from + "' for node '"
                         + pkg[round][param_index_src].edge->getTargetVertex()->getName() + "' to " + file_to
                         + "'. Probably the old file still exists (see earlier errors).";
            OPENMS_LOG_ERROR << msg << std::endl;
            if (ts->isGUIMode()) { QMessageBox::warning(nullptr, tr("File copy failed"), tr(msg.c_str()), QMessageBox::Ok); }
            else
            {
              throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
            }
          }
        }
        update(boundingRect()); // repaint
      }
    }

    finished_ = true;

    __DEBUG_END_METHOD__
  }

  
}
