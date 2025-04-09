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

#include <iostream>

namespace OpenMS
{

  TOPPASSplitterVertex::TOPPASSplitterVertex(const TOPPASSplitterVertex& rhs) = default;

  TOPPASSplitterVertex& TOPPASSplitterVertex::operator=(const TOPPASSplitterVertex& rhs) = default;

  std::unique_ptr<TOPPASVertex> TOPPASSplitterVertex::clone() const
  {
    return std::make_unique<TOPPASSplitterVertex>(*this);
  }

  String TOPPASSplitterVertex::getName() const
  {
    return "SplitterVertex";
  }

  void TOPPASSplitterVertex::run()
  {
    // check if everything ready
    if (!isUpstreamFinished()) 
    {
      return;
    }
    RoundPackages pkg;
    String error_msg("");
    bool success = buildRoundPackages(pkg, error_msg);
    if (!success)
    {
      std::cerr << "Could not retrieve input files from upstream nodes...\n";
      // emit mergeFailed((String("Splitter #") + this->getTopoNr() + " failed. " + error_msg).toQString());
      return;
    }

    output_files_.clear();
    round_counter_ = 0;

    // do the virtual splitting (1 round of N files becomes N rounds of 1 file):
    for (RoundPackages::iterator pkg_it = pkg.begin(); pkg_it != pkg.end();
      ++pkg_it)
    {
      // there can only be one upstream (input) node:
      QStringList files = pkg_it->begin()->second.filenames.get();
      for (QStringList::iterator file_it = files.begin();
        file_it != files.end(); ++file_it)
      {
        RoundPackage new_pkg;
        new_pkg[-1].filenames.push_back(*file_it);
        output_files_.push_back(new_pkg);
        ++round_counter_;
      }
    }

    round_total_ = round_counter_;
    finished_ = true;

    // call all children, proceed in pipeline
    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
    {
      TOPPASVertex* tv = (*it)->getTargetVertex();
      debugOut_(String("Starting child ") + tv->getTopoNr());
      tv->run();
    }
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
