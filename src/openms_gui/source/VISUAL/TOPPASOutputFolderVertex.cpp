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
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

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

  std::string TOPPASOutputFolderVertex::getName() const
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
    if (auto* pipeline = qobject_cast<TOPPASScene*>(scene())) { pipeline->runPipeline(); }
  }

} //namespace OpenMS
