// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtCore>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>

#include <QCoreApplication>

namespace OpenMS
{
  TOPPASOutputFileListVertex::TOPPASOutputFileListVertex() :
    TOPPASVertex(),
    output_folder_name_()
  {
    pen_color_ = Qt::black;
    brush_color_ = Qt::lightGray;
  }

  TOPPASOutputFileListVertex::TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex& rhs) :
    TOPPASVertex(rhs),
    output_folder_name_() // leave empty! otherwise we have conflicting output folder names
  {
    pen_color_ = Qt::black;
    brush_color_ = Qt::lightGray;
  }

  TOPPASOutputFileListVertex::~TOPPASOutputFileListVertex()
  {
  }

  TOPPASOutputFileListVertex& TOPPASOutputFileListVertex::operator=(const TOPPASOutputFileListVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);
    output_folder_name_ = ""; // leave empty! otherwise we have conflicting output folder names

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
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    // display file type(s)
    Map<QString, Size> suffices;
    foreach(QString fn, getFileNames())
    {
      QStringList l = QFileInfo(fn).completeSuffix().split('.');
      QString suf = ((l.size() > 1 && l[l.size() - 2].size() <= 4) ? l[l.size() - 2] + "." : QString()) + l.back(); // take up to two dots as suffix (the first only if its <=4 chars, e.g. we want ".prot.xml" or ".tar.gz", but not "stupid.filename.with.longdots.mzML")
      ++suffices[suf];
    }
    StringList text_l;
    for (Map<QString, Size>::const_iterator sit = suffices.begin(); sit != suffices.end(); ++sit)
    {
      if (suffices.size() > 1)
        text_l.push_back(String(".") + sit->first + "(" + String(sit->second) + ")");
      else
        text_l.push_back(String(".") + sit->first);
    }
    text = ListUtils::concatenate(text_l, " | ").toQString();
    // might get very long, especially if node was not reached yet, so trim
    text = text.left(15) + " ...";
    text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), 35 - (int)(text_boundings.height() / 4.0), text);

    // topo sort number
    painter->drawText(-63.0, -19.0, QString::number(topo_nr_));
    
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

  QPainterPath TOPPASOutputFileListVertex::shape() const
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

    bool dry_run = qobject_cast<TOPPASScene*>(scene())->isDryRun();

    int param_index_src = e->getSourceOutParam();
    int param_index_me = e->getTargetInParam();
    for (Size round = 0; round < pkg.size(); ++round)
    {
      foreach(const QString &f, pkg[round][param_index_src].filenames)
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
              std::cerr << "Could not delete old output file " << String(file_to) << " before overwriting with new one." << std::endl;
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
            LOG_ERROR << "Could not copy tmp output file " << String(file_from) << " to " << String(file_to) << std::endl;
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
