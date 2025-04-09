// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASMergerVertex.h>

#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

#include <iostream>

namespace OpenMS
{
  TOPPASMergerVertex::TOPPASMergerVertex(bool round_based) :
    round_based_mode_(round_based)
  {
  }

  std::unique_ptr<TOPPASVertex> TOPPASMergerVertex::clone() const
  {
    return std::make_unique<TOPPASMergerVertex>(*this);
  }

  String TOPPASMergerVertex::getName() const
  {
    return "MergerVertex";
  }

  void TOPPASMergerVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget)
  {
    TOPPASVertex::paint(painter, option, widget);

    QString text = round_based_mode_ ? "Merge" : "Collect";
    QRectF text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
    painter->drawText(-(int)(text_boundings.width() / 2.0), (int)(text_boundings.height() / 4.0), text);

    if (round_total_ != -1) // draw round number
    {
      text = QString::number(round_counter_) + " / " + QString::number(round_total_);
      text_boundings = painter->boundingRect(QRectF(0, 0, 0, 0), Qt::AlignCenter, text);
      painter->drawText(-(int)(text_boundings.width() / 2.0), 31, text);
    }
  }

  QRectF TOPPASMergerVertex::boundingRect() const
  {
    return QRectF(-41, -41, 82, 82);
  }

  bool TOPPASMergerVertex::roundBasedMode() const
  {
    return round_based_mode_;
  }

  void TOPPASMergerVertex::markUnreachable()
  {
    //only mark as unreachable if all inputs are unreachable. otherwise the dead inputs will just be ignored.
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

  void TOPPASMergerVertex::run()
  {
    //check if everything ready
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
      emit mergeFailed((String("Merger #") + this->getTopoNr() + " failed. " + error_msg).toQString());
      return;
    }

    /// update round status
    Size input_rounds = pkg.size();
    round_total_ = (round_based_mode_ ? (int) input_rounds : 1);  // for round based: take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0; // once round_counter_ reaches round_total_, we are done

    // clear output file list
    output_files_.clear();
    output_files_.resize(round_total_); // #rounds

    // Do the virtual merging (nothing more than reorganizing filenames)
    for (Size round = 0; round < input_rounds; ++round)
    {
      QStringList files;
      // warning: ite->first (i.e. target-in param could be -1,-2,... etc to cover all incoming edges (they all have -1 theoretically - see buildRoundPackages())
      for (RoundPackageConstIt ite = pkg[round].begin();
           ite != pkg[round].end(); ++ite)
      {
        files.append(ite->second.filenames.get()); // concat filenames from all incoming edges
      }
      Size round_index = (round_based_mode_ ? round : 0);
      output_files_[round_index][-1].filenames.append(files); // concat over all rounds (if required)
    }

    round_counter_ = round_total_;
    finished_ = true;

    // call all children, proceed in pipeline
    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
    {
      TOPPASVertex* tv = (*it)->getTargetVertex();
      debugOut_(String("Starting child ") + tv->getTopoNr());
      tv->run();
    }

  }

}
