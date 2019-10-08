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

  TOPPASSplitterVertex::TOPPASSplitterVertex(const TOPPASSplitterVertex& rhs) :
    TOPPASVertex(rhs)
  {
  }

  TOPPASSplitterVertex & TOPPASSplitterVertex::operator=(const TOPPASSplitterVertex& rhs)
  {
    TOPPASVertex::operator=(rhs);

    return *this;
  }

  String TOPPASSplitterVertex::getName() const
  {
    return "SplitterVertex";
  }

  void TOPPASSplitterVertex::run()
  {
    // check if everything ready
    if (!isUpstreamFinished())  return;

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
