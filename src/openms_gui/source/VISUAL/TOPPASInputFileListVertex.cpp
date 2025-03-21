// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

namespace OpenMS
{
  TOPPASInputFileListVertex::TOPPASInputFileListVertex(const QStringList& files)
  {
    setFilenames(files);
  }

  std::unique_ptr<TOPPASVertex> TOPPASInputFileListVertex::clone() const
  {
    return std::make_unique<TOPPASInputFileListVertex>(*this);
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
    TOPPASInputFilesDialog tifd(getFileNames(), cwd_, nullptr);
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

  void TOPPASInputFileListVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget);

    // display number of input files
    QString text = QString::number(getFileNames().size())
                   + " input file"
                   + (getFileNames().size() == 1 ? "" : "s");
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    // display file type(s)
    QStringList text_l = TOPPASVertex::TOPPASFilenames(getFileNames()).getSuffixCounts();
    text = text_l.join(" | ");
    // might get very long, especially if node was not reached yet, so trim
    text = text.left(15) + " ...";
    text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), 35 - (int)(text_boundings.height() / 4.0), text);
  }

  QRectF TOPPASInputFileListVertex::boundingRect() const
  {
    return QRectF(-71, -41, 142, 82);
  }

  bool TOPPASInputFileListVertex::fileNamesValid()
  {
    QStringList fl = getFileNames();
    std::set<std::string> unique_names;
    foreach(const QString& file, fl)
    {
      if (!File::exists(file))
      {
        return false;
      }
      QFileInfo fi(file);
      const auto& [it_unique, was_inserted] = unique_names.insert(fi.canonicalFilePath().toStdString());
      if (!was_inserted) // duplicate
      {
        const auto path = *it_unique;  // working around 'error: reference to local binding 'it_unique' declared in enclosing function' on Clang (capture of structured binding problem)
        OPENMS_LOG_ERROR << "File '" << file.toStdString() << "' (resolved to '" << path << "') appears twice in the input list!" << std::endl;
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

    if (files.empty())
    {
      return;
    }
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
