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

  bool TOPPASMergerVertex::roundBasedMode()
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
    if (!isUpstreamFinished())  return;

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
