// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Johannes Junker, Chris Bielow, Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>

#include <iostream>

namespace OpenMS
{

  TOPPASSplitterVertex::TOPPASSplitterVertex(const TOPPASSplitterVertex& rhs) = default;

  TOPPASSplitterVertex& TOPPASSplitterVertex::operator=(const TOPPASSplitterVertex& rhs) = default;

  std::unique_ptr<TOPPASVertex> TOPPASSplitterVertex::clone() const
  {
    return std::make_unique<TOPPASSplitterVertex>(*this);
  }

  std::string TOPPASSplitterVertex::getName() const
  {
    return "SplitterVertex";
  }

  void TOPPASSplitterVertex::run()
  {
    if (auto* pipeline = qobject_cast<TOPPASScene*>(scene())) { pipeline->runPipeline(); }
  }

  void TOPPASSplitterVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget);

    QString text = "Split";
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    if (round_total_ != -1) // draw round number
    {
      text = QString::number(round_counter_) + " / " + QString::number(round_total_);
      text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
      painter->drawText(-(int)(text_boundings.width() / 2.0), 31, text);
    }
  }

  QRectF TOPPASSplitterVertex::boundingRect() const
  {
    return QRectF(-41, -41, 82, 82);
  }

  void TOPPASSplitterVertex::markUnreachable()
  {
    // only mark as unreachable if all inputs are unreachable. otherwise the dead inputs will just be ignored.
    bool some_input_reachable_ = false;
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      TOPPASVertex* tv = (*it)->getSourceVertex();
      if (tv->isReachable())
      {
        some_input_reachable_ = true;
        break;
      }
    }
    if (!some_input_reachable_)
    {
      TOPPASVertex::markUnreachable();
    }
  }

  void TOPPASSplitterVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
  {
  }

}
