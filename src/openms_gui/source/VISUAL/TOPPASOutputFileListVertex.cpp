// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
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
#include <QtConcurrent/QtConcurrent>
#include <QtWidgets/QMessageBox>

#include <QCoreApplication>

namespace OpenMS
{
  TOPPASOutputFileListVertex::TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex& rhs) :
    TOPPASVertex(rhs),
    output_folder_name_() // leave empty! otherwise we have conflicting output folder names
  {
  }

  TOPPASOutputFileListVertex& TOPPASOutputFileListVertex::operator=(const TOPPASOutputFileListVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);
    output_folder_name_ = ""; // leave empty! otherwise we have conflicting output folder names

    return *this;
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

  void TOPPASOutputFileListVertex::setOutputFolderName(const QString& name)
  {
    if (output_folder_name_ != name)
    {
      output_folder_name_ = name;
      emit outputFolderNameChanged(); // enable storing of modified pipeline
    }
  }

  const QString& TOPPASOutputFileListVertex::getOutputFolderName() const
  {
    return output_folder_name_;
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
    for (Size round = 0; round < pkg.size(); ++round)
    {
      foreach(const QString &f, pkg[round][param_index_src].filenames.get())
      {
        if (!dry_run && !File::exists(f))
        {
          OPENMS_LOG_ERROR << "The file '" << String(f) << "' does not exist!" << std::endl;
          throw Exception::FileNotFound(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, f.toStdString());
        }
        QString new_file = full_dir.toQString()
            + QDir::separator()
            + File::basename(f).toQString();

        // remove "_tmp<number>" if its a suffix
        QRegExp rx("_tmp\\d+$");
        int tmp_index = rx.indexIn(new_file);
        if (tmp_index != -1)
        {
          new_file = new_file.left(tmp_index);
        }

        // get file type and rename
        FileTypes::Type ft = FileTypes::UNKNOWN;
        if (!dry_run)
        {
          TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(e->getSourceVertex());
          if (ttv)
          {
            QVector<TOPPASToolVertex::IOInfo> source_output_files;
            ttv->getOutputParameters(source_output_files);
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
                // look for <name>_type (more generic, supporting more than just 'out')
                if (p.exists(out_type))
                {
                  ft = FileTypes::nameToType(p.getValue(out_type));
                }
              }

            }
          }
        }

        if (ft != FileTypes::UNKNOWN)
        { // replace old suffix by new suffix
          String new_suffix = String(".") + FileTypes::typeToName(ft);
          if (!new_file.endsWith(new_suffix.toQString())) new_file = (File::removeExtension(new_file) + new_suffix).toQString();
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
      for (Size round = 0; round < pkg.size(); ++round)
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
              String msg = "Error: Could not remove old output file '" + String(file_to) + "' for node '" + pkg[round][param_index_src].edge->getTargetVertex()->getName() + "' in preparation to write the new one. Please make sure the file is not open in other applications and try again.";
              OPENMS_LOG_ERROR << msg << std::endl;
              if (ts->isGUIMode()) QMessageBox::warning(nullptr, tr("File removing failed"), tr(msg.c_str()), QMessageBox::Ok);
              else throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
            }
          }

          // running the copy() in an extra thread, such that the GUI stays responsive
          QFuture<bool> future = QtConcurrent::run(copy_, file_from, file_to);
          QMutex mutex; mutex.lock();
          QWaitCondition qwait;
          while (!future.isFinished())
          {
            qApp->processEvents(); // GUI responsiveness
            qwait.wait(&mutex, 25); // block for 25ms (enough for GUI responsiveness), so CPU usage remains low
          }
          mutex.unlock();

          if (bool(future.result()))
          {
            ++files_written_;
            emit outputFileWritten(file_to);
          }
          else
          {
            String msg = "Error: Could not copy temporary output file '" + String(file_to) + "' for node '" + pkg[round][param_index_src].edge->getTargetVertex()->getName() + "' to " + String(file_to) + "'. Probably the old file still exists (see earlier errors).";
            OPENMS_LOG_ERROR << msg << std::endl;
            if (ts->isGUIMode()) QMessageBox::warning(nullptr, tr("File copy failed"), tr(msg.c_str()), QMessageBox::Ok);
            else throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, msg); // fail hard for ExecutePipeline
          }
        }
        update(boundingRect()); // repaint
      }
    }

    finished_ = true;

    __DEBUG_END_METHOD__
  }

  bool TOPPASOutputFileListVertex::copy_(const QString& from, const QString& to)
  {
    return QFile::copy(from, to);
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
    GUIHelpers::openFolder(path);
  }

  String TOPPASOutputFileListVertex::getFullOutputDirectory() const
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
    String dir = String(ts->getOutDir()).substitute("\\", "/");
    return QDir::cleanPath((dir.ensureLastChar('/') + getOutputDir()).toQString());
  }

  String TOPPASOutputFileListVertex::getName() const
  {
    return "OutputVertex";
  }

  String TOPPASOutputFileListVertex::getOutputDir() const
  {
    TOPPASEdge* e = *inEdgesBegin();
    TOPPASVertex* tv = e->getSourceVertex();
    String dir;
    if (output_folder_name_.isEmpty()) {
      // create meaningful output name using vertex + TOPP name + output parameter, e.g. "010-FileConverter-out"
      dir = String("TOPPAS_out") + String(QDir::separator()) + get3CharsNumber_(topo_nr_) + "-"
                                                             + tv->getName() + "-" 
                                                             + e->getSourceOutParamName().remove(':');
    }
    else
    {
      dir = String("TOPPAS_out") + String(QDir::separator()) + output_folder_name_;
    }
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
