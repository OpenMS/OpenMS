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

#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <iostream>

namespace OpenMS
{
#ifdef TOPPAS_DEBUG
  int TOPPASVertex::global_debug_indent_ = 0;
#endif

    int TOPPASVertex::TOPPASFilenames::size() const
  {
    return filenames_.size();
  }

  const QStringList& TOPPASVertex::TOPPASFilenames::get() const
  {
    return filenames_;
  }
  const QString& TOPPASVertex::TOPPASFilenames::operator[](int i) const
  {
    return filenames_[i];
  }

  void TOPPASVertex::TOPPASFilenames::set(const QStringList& filenames)
  {
    filenames_.clear();
    this->append(filenames);
  }

  void TOPPASVertex::TOPPASFilenames::set(const QString& filename, int i)
  {
    check_(filename);
    filenames_[i] = filename;
  }

  void TOPPASVertex::TOPPASFilenames::push_back(const QString& filename)
  {
    check_(filename);
    filenames_.push_back(filename);
  }

  void TOPPASVertex::TOPPASFilenames::append(const QStringList& filenames)
  {
    foreach(const QString& fn, filenames)
    {
      check_(fn);
      push_back(fn);
    }
  }

  void TOPPASVertex::TOPPASFilenames::check_(const QString& filename)
  {
    const Size max_filename_length = 255;
#ifdef OPENMS_WINDOWSPLATFORM
    // on Windows, constraint applies to the whole path:
    if (filename.count() > max_filename_length)
    {
      throw Exception::FileNameTooLong(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename.toStdString(), max_filename_length);
    }
#else
    // on Unix (Linux/Mac), constraint applies to the file name:
    if (File::basename(filename).size() > max_filename_length)
    {
      throw Exception::FileNameTooLong(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, File::basename(filename), max_filename_length);
    }
#endif
  }


  TOPPASVertex::TOPPASVertex() :
    QObject(),
    QGraphicsItem(),
    in_edges_(),
    out_edges_(),
    edge_being_created_(false),
    pen_color_(),
    brush_color_(),
    dfs_color_(DFS_WHITE),
    dfs_parent_(nullptr),
    topo_sort_marked_(false),
    topo_nr_(0),
    round_total_(-1),
    round_counter_(0),
    finished_(false),
    reachable_(true),
    allow_output_recycling_(false)
  {
    setFlag(QGraphicsItem::ItemIsSelectable, true);
    setZValue(42);
  }

  TOPPASVertex::TOPPASVertex(const TOPPASVertex& rhs) :
    QObject(),
    QGraphicsItem(),
    // do not copy pointers to edges
    in_edges_(/*rhs.in_edges_*/),
    out_edges_(/*rhs.out_edges_*/),
    edge_being_created_(rhs.edge_being_created_),
    pen_color_(rhs.pen_color_),
    brush_color_(rhs.brush_color_),
    dfs_color_(rhs.dfs_color_),
    dfs_parent_(rhs.dfs_parent_),
    topo_sort_marked_(rhs.topo_sort_marked_),
    topo_nr_(rhs.topo_nr_),
    round_total_(rhs.round_total_),
    round_counter_(rhs.round_counter_),
    finished_(rhs.finished_),
    reachable_(rhs.reachable_),
    allow_output_recycling_(rhs.allow_output_recycling_)
  {
    setFlag(QGraphicsItem::ItemIsSelectable, true);
    setZValue(42);

    setPos(rhs.pos());
  }

  TOPPASVertex::~TOPPASVertex()
  {
  }

  TOPPASVertex& TOPPASVertex::operator=(const TOPPASVertex& rhs)
  {
    in_edges_ = rhs.in_edges_;
    out_edges_ = rhs.out_edges_;
    edge_being_created_ = rhs.edge_being_created_;
    pen_color_ = rhs.pen_color_;
    brush_color_ = rhs.brush_color_;
    dfs_color_ = rhs.dfs_color_;
    dfs_parent_ = rhs.dfs_parent_;
    topo_sort_marked_ = rhs.topo_sort_marked_;
    topo_nr_ = rhs.topo_nr_;

    round_total_ = rhs.round_total_;
    round_counter_ = rhs.round_counter_;
    finished_ = rhs.finished_;
    reachable_ = rhs.reachable_;
    allow_output_recycling_ = rhs.allow_output_recycling_;

    setPos(rhs.pos());

    return *this;
  }

  bool TOPPASVertex::isUpstreamFinished() const
  {
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      TOPPASVertex * tv = (*it)->getSourceVertex();
      if (!tv->isFinished())
      {
        // some tool that we depend on has not finished execution yet --> do not start yet
        debugOut_("Not run (parent not finished)");

        __DEBUG_END_METHOD__
        return false;
      }
    }
    //std::cerr << "upstream of " << this->getTopoNr() << " is ready!\n";
    return true;
  }

  bool TOPPASVertex::buildRoundPackages(RoundPackages& pkg, String& error_msg) // check all incoming edges for this node and construct the package
  {
    if (inEdgesBegin() == inEdgesEnd())
    {
      error_msg = "buildRoundPackages() called on vertex with no input edges!\n";
      LOG_ERROR << error_msg;
      return false;
    }

    // -- determine number of rounds from incoming edges
    int round_common = -1; // number of rounds common to all
    int no_recycle_count = 0; // number of edges that do NOT do recycling (there needs to be at least one)
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it) // all incoming edges should have the same number of rounds (or should be set to 'recycle') !
    {
      TOPPASVertex* tv = (*it)->getSourceVertex();
      if (tv->allow_output_recycling_)
      {
        continue;
      }

      ++no_recycle_count;
      if (round_common == -1)
      {
        round_common = tv->round_total_; // first non-recycler sets the pace
      }

      if (round_common != tv->round_total_)
      {
        error_msg = String("Number of rounds for incoming edges of node #") + this->getTopoNr() + " are not equal. No idea on how to combine them! Did you want to recycle its input?\n";
        std::cerr << error_msg;
        return false;
      }
    }

    // -- we demand at least one node with no recycling to allow to determine number of rounds
    if (no_recycle_count == 0)
    {
      error_msg = String("Number of rounds of node #") + this->getTopoNr() + " cannot be determined since all input nodes have recycling enabled. Disable for at least one input!\n";
      std::cerr << error_msg;
      return false;
    }

    // -- check if rounds from recyling nodes are an integer part of total rounds, i.e. total_rounds = X * node_rounds, X from N+
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it) // look at all all recycling edges
    {
      TOPPASVertex * tv = (*it)->getSourceVertex();
      if (!tv->allow_output_recycling_)
        continue;

      if (round_common % tv->round_total_ != 0) // modulo should be 0, if not ...
      {
        error_msg = String(tv->round_total_) + " rounds for incoming edges of node #" + this->getTopoNr() + " are recycled to meet a total of " + round_common + " rounds. But modulo is not 0. No idea on how to combine them! Adapt the number of input files?\n";
        std::cerr << error_msg;
        return false;
      }
    }

    if (round_common <= 0)
    {
      error_msg =  "Number of input rounds is 0 or negative. This cannot be! Aborting!\n";
      std::cerr << error_msg;
      return false;
    }
    pkg.clear();
    pkg.resize(round_common);

    // all incoming edges should have the same number of rounds!
    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      TOPPASVertex* tv_upstream = (*it)->getSourceVertex();

      // fill files for each round
      int param_index_src_out = (*it)->getSourceOutParam();
      int param_index_tgt_in = (*it)->getTargetInParam();
      for (int round = 0; round < round_common; ++round)
      {
        VertexRoundPackage rpg;
        rpg.edge = *it;
        int upstream_round = round;
        if (tv_upstream->allow_output_recycling_ && upstream_round >= tv_upstream->round_total_)
        {
          upstream_round %= tv_upstream->round_total_;
        }
        rpg.filenames.set(tv_upstream->getFileNames(param_index_src_out, upstream_round));

        // hack for merger vertices, as they have multiple incoming edges with -1 as index
        while (pkg[round].count(param_index_tgt_in))
        {
          --param_index_tgt_in; // find free slot, i.e. -2, -3 ....
        }

        pkg[round][param_index_tgt_in] = rpg; // index by incoming edge number
      }
    }
    return true;
  }

  QStringList TOPPASVertex::getFileNames(int param_index, int round) const
  {
    if ((Size)round >= output_files_.size())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, round, output_files_.size());
    }
    RoundPackage rp = output_files_[round];
    if (rp.find(param_index) == rp.end())
    {
      throw Exception::IndexOverflow(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, param_index, rp.size());  // index could be larger (its a map, but nevertheless)
    }
    //String s = String(rp[param_index].filenames.join("\" \""));
    return rp[param_index].filenames.get();
  }

  QStringList TOPPASVertex::getFileNames() const
  {
    // concatenate over all rounds
    QStringList fl;

    for (Size r = 0; r < output_files_.size(); ++r)
    {
      for (RoundPackage::const_iterator it  = output_files_[r].begin();
           it != output_files_[r].end();
           ++it)
      {
        fl.append(it->second.filenames.get());
      }
    }
    return fl;
  }

  TOPPASVertex::SUBSTREESTATUS TOPPASVertex::getSubtreeStatus() const
  {
    if (!this->isFinished()) return TV_UNFINISHED;

    if (!this->isUpstreamFinished()) return TV_UNFINISHED_INBRANCH; // only looks for immediate predecessors!

    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
    {
      SUBSTREESTATUS status = (*it)->getTargetVertex()->getSubtreeStatus();
      if (status != TV_ALLFINISHED) return status;
    }
    return TV_ALLFINISHED;
  }

  const TOPPASVertex::RoundPackages& TOPPASVertex::getOutputFiles() const
  {
    return output_files_;
  }

  void TOPPASVertex::mousePressEvent(QGraphicsSceneMouseEvent* e)
  {
    if (!(e->modifiers() & Qt::ControlModifier))
    {
      emit clicked();
    }
  }

  void TOPPASVertex::mouseReleaseEvent(QGraphicsSceneMouseEvent* e)
  {
    if (edge_being_created_)
    {
      emit finishHoveringEdge();
      edge_being_created_ = false;
    }
    else if (e->modifiers() & Qt::ControlModifier)
    {
      QGraphicsItem::mouseReleaseEvent(e);
    }
    else
    {
      emit released();
      // resize scene rect in case item has been moved outside
      const QRectF& scene_rect = scene()->sceneRect();
      const QRectF& items_bounding = scene()->itemsBoundingRect();
      scene()->setSceneRect(scene_rect.united(items_bounding));
    }
  }

  void TOPPASVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e)
  {
    e->ignore();
  }

  void TOPPASVertex::contextMenuEvent(QGraphicsSceneContextMenuEvent* e)
  {
    e->ignore();
  }

  void TOPPASVertex::mouseMoveEvent(QGraphicsSceneMouseEvent* e)
  {
    TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());

    if (isSelected())
    {
      QPointF delta = e->pos() - e->lastPos();
      emit itemDragged(delta.x(), delta.y());
    }
    else
    {
      ts->setActionMode(TOPPASScene::AM_NEW_EDGE);
      moveNewEdgeTo_(e->pos());
    }
  }

  void TOPPASVertex::moveNewEdgeTo_(const QPointF& pos)
  {
    if (!edge_being_created_)
    {
      emit newHoveringEdge(mapToScene(pos));
      edge_being_created_ = true;
    }

    emit hoveringEdgePosChanged(mapToScene(pos));
  }

  TOPPASVertex::ConstEdgeIterator TOPPASVertex::outEdgesBegin() const
  {
    return out_edges_.begin();
  }

  TOPPASVertex::ConstEdgeIterator TOPPASVertex::outEdgesEnd() const
  {
    return out_edges_.end();
  }

  TOPPASVertex::ConstEdgeIterator TOPPASVertex::inEdgesBegin() const
  {
    return in_edges_.begin();
  }

  TOPPASVertex::ConstEdgeIterator TOPPASVertex::inEdgesEnd() const
  {
    return in_edges_.end();
  }

  void TOPPASVertex::addInEdge(TOPPASEdge* edge)
  {
    in_edges_.push_back(edge);
  }

  void TOPPASVertex::addOutEdge(TOPPASEdge* edge)
  {
    out_edges_.push_back(edge);
  }

  void TOPPASVertex::removeInEdge(TOPPASEdge* edge)
  {
    in_edges_.removeAll(edge);
  }

  void TOPPASVertex::removeOutEdge(TOPPASEdge* edge)
  {
    out_edges_.removeAll(edge);
  }

  TOPPASVertex::DFS_COLOR TOPPASVertex::getDFSColor()
  {
    return dfs_color_;
  }

  void TOPPASVertex::setDFSColor(DFS_COLOR color)
  {
    dfs_color_ = color;
  }

  TOPPASVertex * TOPPASVertex::getDFSParent()
  {
    return dfs_parent_;
  }

  void TOPPASVertex::setDFSParent(TOPPASVertex* parent)
  {
    dfs_parent_ = parent;
  }

  Size TOPPASVertex::incomingEdgesCount()
  {
    return in_edges_.size();
  }

  Size TOPPASVertex::outgoingEdgesCount()
  {
    return out_edges_.size();
  }

  void TOPPASVertex::inEdgeHasChanged()
  {
    // (overridden behavior in output and tool vertices)
    qobject_cast<TOPPASScene*>(scene())->setChanged(true);
    emit somethingHasChanged();
  }

  void TOPPASVertex::outEdgeHasChanged()
  {
    // (overridden behavior in input and tool vertices)
    qobject_cast<TOPPASScene*>(scene())->setChanged(true);
    emit somethingHasChanged();
  }

  bool TOPPASVertex::isTopoSortMarked()
  {
    return topo_sort_marked_;
  }

  void TOPPASVertex::setTopoSortMarked(bool b)
  {
    topo_sort_marked_ = b;
  }

  UInt TOPPASVertex::getTopoNr()
  {
    return topo_nr_;
  }

  void TOPPASVertex::setTopoNr(UInt nr)
  {
    // (overridden in tool and output vertices)
    topo_nr_ = nr;
  }

  String TOPPASVertex::get3CharsNumber_(UInt number) const
  {
    String num_str(number);
    num_str.fillLeft('0', 3);
    return num_str;
  }

  void TOPPASVertex::reset(bool /*reset_all_files*/)
  {
    __DEBUG_BEGIN_METHOD__

    round_total_ = -1;
    round_counter_ = 0;

    finished_ = false;
    reachable_ = true;

    update(boundingRect());

    __DEBUG_END_METHOD__
  }

  bool TOPPASVertex::isFinished() const
  {
    return finished_;
  }

  void TOPPASVertex::run()
  {
    throw Exception::NotImplemented(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
  }

  bool TOPPASVertex::invertRecylingMode()
  {
    allow_output_recycling_ = !allow_output_recycling_;
    emit parameterChanged(true); // using 'true' is very conservative but safe. One could override this in child classes.
    return allow_output_recycling_;
  }

  bool TOPPASVertex::isRecyclingEnabled() const
  {
    return allow_output_recycling_;
  }

  void TOPPASVertex::setRecycling(const bool is_enabled)
  {
    if (allow_output_recycling_ == is_enabled) return; // nothing changed

    invertRecylingMode();
  }

  bool TOPPASVertex::allInputsReady()
  {
    __DEBUG_BEGIN_METHOD__

    for (ConstEdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
    {
      TOPPASVertex* tv = qobject_cast<TOPPASVertex*>((*it)->getSourceVertex());
      if (tv && !tv->isFinished())
      {
        // some (reachable) tool that we depend on has not finished execution yet --> do not start yet
        __DEBUG_END_METHOD__
        return false;
      }
    }

    __DEBUG_END_METHOD__
    return true;
  }

  void TOPPASVertex::markUnreachable()
  {
    reachable_ = false;
    for (ConstEdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
    {
      TOPPASVertex* tv = (*it)->getTargetVertex();
      if (tv->reachable_)
      {
        tv->markUnreachable();
      }
    }
  }

  bool TOPPASVertex::isReachable()
  {
    return reachable_;
  }

}
