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
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

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

  std::string TOPPASOutputFileListVertex::getName() const
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
    if (auto* pipeline = qobject_cast<TOPPASScene*>(scene())) { pipeline->runPipeline(); }
  }

  
}
