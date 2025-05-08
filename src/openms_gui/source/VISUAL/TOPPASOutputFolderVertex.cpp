// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>

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
  std::unique_ptr<TOPPASVertex> TOPPASOutputFolderVertex::clone() const
  {
    return std::make_unique<TOPPASOutputFolderVertex>(*this);
  }

  String TOPPASOutputFolderVertex::getName() const
    {
      return "OutputFolderVertex";
    }

  void TOPPASOutputFolderVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget);

    QString text = "output folder";
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    // display file type(s) -- not implemented for output folder (we would need to inspect the directory)
    // ...

    // output folder name
    painter->drawText(painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, output_folder_name_).width()/-2, -25, output_folder_name_);
  }

  QRectF TOPPASOutputFolderVertex::boundingRect() const
  {
    return QRectF(-71, -41, 142, 82);
  }

  void TOPPASOutputFolderVertex::run()
  {
    // copy tmp files to output dir

    // get incoming edge and preceding vertex
    TOPPASEdge* e = *inEdgesBegin();
    RoundPackages pkg = e->getSourceVertex()->getOutputFiles(); // this will be output folders
    if (pkg.empty())
    {
      std::cerr << "A problem occurred while grabbing files from the parent tool. This is a bug, please report it!" << std::endl;
      return;
    }

    String full_dir = createOutputDir(); // create output dir

    round_total_ = (int) pkg.size(); // take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0; // once round_counter_ reaches round_total_, we are done

    // clear output file list
    output_files_.clear();
    output_files_.resize(pkg.size()); // #rounds

    files_total_ = 0;
    files_written_ = 0;

    const TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    bool dry_run = ts->isDryRun();

    int param_index_src = e->getSourceOutParam();
    int param_index_me = e->getTargetInParam();
    for (Size round = 0; round < round_total_; ++round)
    {
      for (const QString &src_dir : pkg[round][param_index_src].filenames.get())
      {
        if (!dry_run && !File::isDirectory(src_dir))
        {
          OPENMS_LOG_ERROR << "The directory '" << String(src_dir) << "' does not exist!" << std::endl;
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, src_dir.toStdString());
        }
        String new_dir = full_dir;
        // if its only one round, do not create a subfolder
        if (round_total_ > 1)
        { // the previous TOPP tool node should have placed all files in a subfolder for each round. Use that name
          String last_subfolder = QDir(src_dir).dirName();
          new_dir += last_subfolder.ensureLastChar('/');
        }
        output_files_[round][param_index_me].filenames.push_back(new_dir.toQString());
        // find number of files in 'src_dir'
        QDir dir(src_dir);
        auto nr_of_files = dir.entryInfoList(QDir::Files).size();
        if (!dry_run && nr_of_files == 0)
        { 
          String msg = "Output directory '" + e->getSourceOutParamName() + "' did not yield any files!";
          if (ts->isGUIMode()) 
          {
            QMessageBox::warning(nullptr, tr("No files found"), msg.toQString(), QMessageBox::Ok);
          }
          else
          {
            throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
          }
        }
        files_total_ += nr_of_files; // count number of files in 'src_dir'
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
          String dir_from = pkg[round][param_index_src].filenames[i];
          String dir_to   = output_files_[round][param_index_me].filenames[i];
          // create the output directory (if it does not exist)
          if (!File::makeDir(dir_to))
          {
              String msg = "Could not create output directory '" + dir_to + "' for node '"
                           + pkg[round][param_index_src].edge->getTargetVertex()->getName() + "' . Please make sure the path is writable.";
              OPENMS_LOG_ERROR << msg << std::endl;
              if (ts->isGUIMode())
              {
                QMessageBox::warning(nullptr, tr("Directory creation failed"), tr(msg.c_str()), QMessageBox::Ok);
              }
              else
              {
                throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
              }
          }
          const auto& src_files_to_copy = QDir(dir_from.toQString()).entryInfoList(QDir::Files);
          for (const auto& src_file : src_files_to_copy)
          {
            String file_from = src_file.absoluteFilePath();
            String file_to = dir_to + '/' + src_file.fileName();
            if (File::exists(file_to) // someone may have deleted the file in the meantime, which is fine
              && !QFile::remove(file_to.toQString())) // remove old file (would fail if file does not exist)
            {
              String msg = "Error: Could not remove old output file '" + file_to + "' for node '"
                           + pkg[round][param_index_src].edge->getTargetVertex()->getName()
                           + "' in preparation to write the new one. Please make sure the file is not open in other applications and try again.";
              OPENMS_LOG_ERROR << msg << std::endl;
              if (ts->isGUIMode()) 
              {
                QMessageBox::warning(nullptr, tr("File removing failed"), tr(msg.c_str()), QMessageBox::Ok);
              }
              else
              {
                throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
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
        }
        update(boundingRect()); // repaint after each round
      } // end of round loop
    } // end of dry_run check

    finished_ = true;
  }

} //namespace OpenMS
