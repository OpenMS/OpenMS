// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASInputFilesDialog.h>
#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QSvgRenderer>

namespace OpenMS
{
  TOPPASInputFileListVertex::TOPPASInputFileListVertex() :
    TOPPASVertex(),
    key_()
  {
    pen_color_ = Qt::black;
    brush_color_ = Qt::lightGray;
  }

  TOPPASInputFileListVertex::TOPPASInputFileListVertex(const QStringList& files) :
    TOPPASVertex(),
    key_()
  {
    pen_color_ = Qt::black;
    brush_color_ = Qt::lightGray;
    setFilenames(files);
  }

  TOPPASInputFileListVertex::TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs) :
    TOPPASVertex(rhs),
    key_()
  {
    pen_color_ = Qt::black;
    brush_color_ = Qt::lightGray;
    output_files_ = rhs.output_files_; // copy input file paths, too
  }

  TOPPASInputFileListVertex::~TOPPASInputFileListVertex()
  {

  }

  TOPPASInputFileListVertex& TOPPASInputFileListVertex::operator=(const TOPPASInputFileListVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);

    key_ = rhs.key_;
    output_files_ = rhs.output_files_; // copy input file paths, too

    return *this;
  }

  String TOPPASInputFileListVertex::getName() const
  {
    return "InputVertex";
  }

  void TOPPASInputFileListVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent * /*e*/)
  {
    showFilesDialog();
  }

  void TOPPASInputFileListVertex::showFilesDialog()
  {
    TOPPASInputFilesDialog tifd(getFileNames(), cwd_);
    if (tifd.exec())
    {
      QStringList updated_filelist;
      tifd.getFilenames(updated_filelist);
      if (getFileNames() != updated_filelist)
      { // files were changed
        setFilenames(updated_filelist); // to correct filenames (separators etc)
        qobject_cast<TOPPASScene *>(scene())->updateEdgeColors();

        // update cwd
        cwd_ = tifd.getCWD();

        emit parameterChanged(true); // aborts the pipeline (if running) and resets downstream nodes
      }
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

    // display number of input files
    QString text = QString::number(getFileNames().size())
                   + " input file"
                   + (getFileNames().size() == 1 ? "" : "s");
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
    text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), 35 - (int)(text_boundings.height() / 4.0), text);

    // topological sort number
    qreal x_pos = -63.0;
    qreal y_pos = -19.0;
    painter->drawText(x_pos, y_pos, QString::number(topo_nr_));

    // recycling status
    if (this->allow_output_recycling_)
    {
      painter->setPen(Qt::green);
      QSvgRenderer* svg_renderer = new QSvgRenderer(QString(":/Recycling_symbol.svg"), nullptr);
      svg_renderer->render(painter, QRectF(-7, -32, 14, 14));
    }

  }

  QRectF TOPPASInputFileListVertex::boundingRect() const
  {
    return QRectF(-71, -41, 142, 82);
  }

  QPainterPath TOPPASInputFileListVertex::shape() const
  {
    QPainterPath shape;
    shape.addRoundRect(-71.0, -41.0, 142.0, 81.0, 20, 20);
    return shape;
  }

  bool TOPPASInputFileListVertex::fileNamesValid()
  {
    QStringList fl = getFileNames();
    foreach(const QString& file, fl)
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
    for (int i = 0; i < fl.size(); ++i) // collect unique directories
    {
      QFileInfo fi(fl[i]);
      directories.insert(String(QFileInfo(fi.canonicalFilePath()).path()));
    }

    // open them
    for (std::set<String>::const_iterator it = directories.begin(); it != directories.end(); ++it)
    {
      QString path = QDir::toNativeSeparators(it->toQString());
      GUIHelpers::openFolder(path);
    }
  }

  void TOPPASInputFileListVertex::run()
  {
    round_total_   = (int) output_files_.size(); // for now each file is one round; for the future we might allow to create blocks of files (e.g. for replicate measurements)
    round_counter_ = (int) round_total_;

    this->finished_ = true; // input node is ready to go (file check was already done)

    //std::cerr << "#" << this->getTopoNr() << " set #rounds: " << round_total_ << "\n";

    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
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
  }

  const QString& TOPPASInputFileListVertex::getKey()
  {
    return key_;
  }

  void TOPPASInputFileListVertex::setFilenames(const QStringList& files)
  {
    output_files_.clear();

    if (files.empty()) return;

    output_files_.resize(files.size()); // for now, assume one file per round (we could later extend that)
    for (int f = 0; f < files.size(); ++f)
    {
      output_files_[f][-1].filenames.push_back(QDir::toNativeSeparators(files[f]));
    }

    setToolTip(files.join("\n"));

    // set current working dir when opening files to the last file
    cwd_ = File::path(files.back()).toQString();
  }

  void TOPPASInputFileListVertex::outEdgeHasChanged()
  {
    reset();
    qobject_cast<TOPPASScene*>(scene())->updateEdgeColors();
    TOPPASVertex::outEdgeHasChanged();
  }

}
