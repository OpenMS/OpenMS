// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/PIPELINE/PipelineExecutor.h>
#include <OpenMS/APPLICATIONS/PIPELINE/PipelineFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/SYSTEM/SystemSettings.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASVertexNameDialog.h>
#include <OpenMS/VISUAL/MISC/Qt5Port.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <QApplication>
#include <QtCore/QDir>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QSet>
#include <QtCore/QTemporaryDir>
#include <QtCore/QTextStream>
#include <QtWidgets/QMessageBox>
#include <map>
#include <optional>
#include <set>
#include <thread>
#include <vector>

namespace OpenMS
{


struct TOPPASScene::ExecutionState
{
  std::unique_ptr<PipelineExecutor> executor;
  std::thread worker;
  PipelineExecutor::Result previous_result;
  PipelineExecutor::Result worker_result;
  std::map<TOPPASVertex*, Size> ids;
  // Each run owns an exclusively created parent directory. Keeping these RAII
  // owners separate from the engine allows retained ancestors to outlive reruns.
  std::vector<std::unique_ptr<QTemporaryDir>> temporary_runs;
  Size next_id {0};
  Size generation {0};

  void acceptResult()
  {
    // Cancellation may happen before run() initializes its result. Preserve
    // unaffected cached ancestors; invalidated descendants were already erased.
    previous_result.exit_code = worker_result.exit_code;
    previous_result.error_message = std::move(worker_result.error_message);
    for (auto& [id, result] : worker_result.nodes)
    {
      previous_result.nodes[id] = std::move(result);
    }
    worker_result = {};
    discardUnusedTemporaryRuns();
  }

  void discardUnusedTemporaryRuns()
  {
    // Called only after joining the worker. Never infer ownership from a name
    // prefix or delete directories belonging to another scene or application.
    std::erase_if(temporary_runs, [this](const auto& run) {
      const QDir directory(run->path());
      for (const auto& [id, result] : previous_result.nodes)
      {
        for (const auto& round : result.outputs)
        {
          for (const auto& [port, files] : round)
          {
            for (const auto& file : files)
            {
              const QString relative = directory.relativeFilePath(QFileInfo(toQString(file)).absoluteFilePath());
              if (! QDir::isAbsolutePath(relative) && relative != ".." && ! relative.startsWith("../")) { return false; }
            }
          }
        }
      }
      return true;
    });
  }
};

namespace
{
  // Apply a result on the GUI thread. Graphics items never schedule descendants.
  void presentNodeResult(TOPPASVertex* vertex, const PipelineExecutor::NodeResult& result)
  {
    TOPPASVertex::RoundPackages outputs(result.outputs.size());
    auto* tool = qobject_cast<TOPPASToolVertex*>(vertex);
    const auto ports = tool ? tool->getOutputParameters() : QVector<TOPPASToolVertex::IOInfo> {};
    for (Size round = 0; round < result.outputs.size(); ++round)
    {
      for (const auto& [name, files] : result.outputs[round])
      {
        Int index = -1;
        for (int i = 0; i < ports.size(); ++i)
        {
          if (ports[i].param_name == name)
          {
            index = i;
            break;
          }
        }
        for (const auto& file : files)
        {
          outputs[round][index].filenames.push_back(toQString(file));
        }
      }
    }
    const bool finished = result.state == PipelineExecutor::State::SUCCEEDED;
    // Input filenames are editable definition data. A completed input snapshot
    // may reach this adapter while an edit is cancelling the old run; preserve
    // the current selection instead of restoring the snapshot's filenames.
    if (! qobject_cast<TOPPASInputFileListVertex*>(vertex) || finished)
    {
      vertex->setExecutionResult(qobject_cast<TOPPASInputFileListVertex*>(vertex) ? vertex->getOutputFiles() : outputs,
                                 result.completed_rounds, result.total_rounds, finished);
    }
    if (tool)
    {
      switch (result.state)
      {
        case PipelineExecutor::State::RUNNING:
          tool->toolStartedSlot();
          break;
        case PipelineExecutor::State::SUCCEEDED:
          tool->toolFinishedSlot();
          break;
        case PipelineExecutor::State::FAILED:
        case PipelineExecutor::State::CANCELLED:
          tool->toolFailedSlot();
          break;
        case PipelineExecutor::State::PENDING:
        case PipelineExecutor::State::BLOCKED:
          break;
      }
    }
    if (auto* output = qobject_cast<TOPPASOutputVertex*>(vertex))
    {
      Size files = 0;
      for (const auto& round : result.outputs)
      {
        for (const auto& bundle : round)
        {
          files += bundle.second.size();
        }
      }
      output->setOutputProgress(finished ? files : 0, files);
    }
  }
} // namespace

TOPPASScene::TOPPASScene(QObject* parent, const QString& tmp_path, bool gui):
    QGraphicsScene(parent),
    action_mode_(AM_NEW_EDGE),
    vertices_(),
    edges_(),
    hover_edge_(nullptr),
    potential_target_(nullptr),
    file_name_(),
    tmp_path_(tmp_path + "/" + toQString(File::getUniqueName(false))),
    gui_(gui),
    out_dir_(toQString(SystemSettings::getUserDirectory())),
    changed_(false),
    running_(false),
    error_occured_(false),
    user_specified_out_dir_(false),
    clipboard_(nullptr),
    allowed_threads_(1),
    execution_(std::make_unique<ExecutionState>())
{
  /*	ATTENTION!

          The following line is important! Without it, we get
          hard-to-reproduce segmentation faults and
          "pure virtual method calls" due to a bug in Qt!

          (http://lists.trolltech.com/qt4-preview-feedback/2006-09/thread00124-0.html)
  */
  setItemIndexMethod(QGraphicsScene::NoIndex);
  // Parameter refresh writes an INI here before the first execution.
  if (! QDir().mkpath(tmp_path_)) { throw Exception::FileNotWritable(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, fromQString(tmp_path_)); }
}

  TOPPASScene::~TOPPASScene()
  {
    // Do not propagate scene changes while tearing down: setChanged() below would emit
    // mainWindowNeedsUpdate(), whose receiver (the main window) may already be partially
    // destroyed during application shutdown. This mirrors the per-item blockSignals() below.
    blockSignals(true);
    abortPipeline();
    // Delete all items in a controlled way:
    for (TOPPASVertex* vertex : vertices_)
    {
      vertex->blockSignals(true); // do not propagate changes, remove output files, etc..
      vertex->setSelected(true);
    }
    for (TOPPASEdge* edge : edges_)
    {
      edge->blockSignals(true); // do not propagate changes, remove output files, etc..
      edge->setSelected(true);
    }
    removeSelected();
    execution_->temporary_runs.clear();
    if (File::exists(fromQString(tmp_path_))) { File::removeDirRecursively(fromQString(tmp_path_)); }
  }

  void TOPPASScene::setActionMode(ActionMode mode)
  {
    action_mode_ = mode;
  }

  TOPPASScene::ActionMode TOPPASScene::getActionMode()
  {
    return action_mode_;
  }

  TOPPASScene::VertexIterator TOPPASScene::verticesBegin()
  {
    return vertices_.begin();
  }

  TOPPASScene::VertexIterator TOPPASScene::verticesEnd()
  {
    return vertices_.end();
  }

  TOPPASScene::EdgeIterator TOPPASScene::edgesBegin()
  {
    return edges_.begin();
  }

  TOPPASScene::EdgeIterator TOPPASScene::edgesEnd()
  {
    return edges_.end();
  }

  void TOPPASScene::addVertex(TOPPASVertex* tv)
  {
    abortPipeline();
    execution_->ids.emplace(tv, execution_->next_id++);
    vertices_.push_back(tv);
    addItem(tv);
  }

  void TOPPASScene::addEdge(TOPPASEdge* te)
  {
    abortPipeline();
    edges_.push_back(te);
    addItem(te);
  }

  void TOPPASScene::itemClicked()
  {

  }

  void TOPPASScene::itemReleased()
  {
    TOPPASVertex* sender = qobject_cast<TOPPASVertex*>(QObject::sender());
    if (!sender)
    {
      return;
    }

    // deselect all items except for the one under the cursor, but only if no multiple selection
    if (selectedItems().size() <= 1)
    {
      unselectAll();
      sender->setSelected(true);
    }

    snapToGrid();
  }

  void TOPPASScene::updateHoveringEdgePos(const QPointF& new_pos)
  {
    if (!hover_edge_)
    {
      return;
    }

    hover_edge_->setHoverPos(new_pos);

    TOPPASVertex* target = getVertexAt_(new_pos);
    if (target)
    {
      if (target != potential_target_)
      {
        potential_target_ = target;
        bool ev = isEdgeAllowed_(hover_edge_->getSourceVertex(), target);
        if (ev)
        {
          hover_edge_->setColor(Qt::darkGreen);
        }
        else
        {
          hover_edge_->setColor(Qt::red);
        }
      }
    }
    else
    {
      hover_edge_->setColor(Qt::black);
      potential_target_ = nullptr;
    }
  }

  void TOPPASScene::addHoveringEdge(const QPointF& pos)
  {
    TOPPASVertex* sender = qobject_cast<TOPPASVertex*>(QObject::sender());
    if (!sender)
    {
      return;
    }
    TOPPASEdge* new_edge = new TOPPASEdge(sender, pos);
    hover_edge_ = new_edge;
    addEdge(new_edge);
  }

  void TOPPASScene::finishHoveringEdge()
  {
    TOPPASVertex* target = getVertexAt_(hover_edge_->endPos());
    bool remove_edge = false;

    if (target && target != hover_edge_->getSourceVertex())
    {
      hover_edge_->setTargetVertex(target);
      TOPPASVertex* source = hover_edge_->getSourceVertex();

      // check for parameter copy action (only if source is a tool node (--> edge is purple already, user expects this to happen))
      TOPPASToolVertex* tv_source = qobject_cast<TOPPASToolVertex*>(source);
      if ((QGuiApplication::keyboardModifiers() & Qt::ControlModifier) && tv_source)
      {
        TOPPASToolVertex* tv_target = qobject_cast<TOPPASToolVertex*>(target);
        if (!(tv_source && tv_target))
        {
          emit messageReady("Copying parameters is only allowed between Tool nodes! No copy was performed!\n");
        }
        else
        {
          emit messageReady("Transferring parameters between nodes ...\n");
          Param from = tv_source->getParam();
          Param to = tv_target->getParam();
          Param to_old = to; // backup, to compare

          std::stringstream ss;
          Logger::LogStream my_log(new Logger::LogStreamBuf("Transfer", nullptr));
          my_log.insert(ss);
          to.update(from, false, my_log);
          if (to == to_old)
          {
            my_log << "All parameters are up to date! Nothing happened!\n";
          }
          else // update the target parameters
          {
            tv_target->setParam(to);
            abortPipeline();
            setChanged(true); // to allow "Store" of pipeline
            resetDownstream(target);
          }
          //ss << "test test";
          my_log << " ---------------------------------- " << std::endl; // this will cause a flush... removing this line might cause loss(!) of log content!
          my_log.flush(); // bug! this sometimes does not cause the content to be flushed to the stringstream; the cache seems to be inactive as well. also std::endl does not help
          emit messageReady(toQString(ss.str()));
          //std::cerr << ss.str();
        }
        remove_edge = true;
      }
      else if (isEdgeAllowed_(hover_edge_->getSourceVertex(), target))
      {
        source->addOutEdge(hover_edge_);
        target->addInEdge(hover_edge_);
        hover_edge_->setColor(QColor(255, 165, 0));

        connectEdgeSignals(hover_edge_);

        TOPPASIOMappingDialog dialog(hover_edge_);
        if (dialog.firstExec())
        {
          hover_edge_->emitChanged();
        }
        else
        {
          remove_edge = true;
        }
      }
      else
      {
        remove_edge = true;
      }
    }
    else
    {
      remove_edge = true;
    }

    if (remove_edge)
    {
      edges_.removeAll(hover_edge_);
      removeItem(hover_edge_);
      delete hover_edge_;
      hover_edge_ = nullptr;
    }
    else
    {  // edge was added ...
      topoSort();
      updateEdgeColors();
    }
  }

  TOPPASVertex* TOPPASScene::getVertexAt_(const QPointF& pos)
  {
    QList<QGraphicsItem*> target_list = items(pos);

    // return first item that is a vertex
    TOPPASVertex* target = nullptr;
    for (QList<QGraphicsItem*>::iterator it = target_list.begin(); it != target_list.end(); ++it)
    {
      target = dynamic_cast<TOPPASVertex*>(*it);
      if (target)
      {
        break;
      }
    }

    return target;
  }

  void TOPPASScene::copySelected()
  {
    TOPPASScene* tmp_scene = new TOPPASScene(nullptr, this->getTempDir(), false);
    std::map<TOPPASVertex*, TOPPASVertex*> vertex_map;

    for (TOPPASVertex* v : vertices_)
    {
      if (!v->isSelected())
      {
        continue;
      }

      TOPPASVertex* new_v = v->clone().release();

      vertex_map[v] = new_v;
      tmp_scene->addVertex(new_v);
    }

    for (TOPPASEdge* e : edges_)
    {
      if (!e->isSelected())
      {
        continue;
      }

      //check if both source and target node were also selected (otherwise don't copy)
      TOPPASVertex* old_source = e->getSourceVertex();
      TOPPASVertex* old_target = e->getTargetVertex();
      if (vertex_map.find(old_source) == vertex_map.end() || vertex_map.find(old_target) == vertex_map.end()) { continue; }

      TOPPASEdge* new_e = new TOPPASEdge();
      TOPPASVertex* new_source = vertex_map[old_source];
      TOPPASVertex* new_target = vertex_map[old_target];
      new_e->setSourceVertex(new_source);
      new_e->setTargetVertex(new_target);
      new_e->setSourceOutParam(e->getSourceOutParam());
      new_e->setTargetInParam(e->getTargetInParam());
      new_source->addOutEdge(new_e);
      new_target->addInEdge(new_e);

      tmp_scene->addEdge(new_e);
    }

    emit selectionCopied(tmp_scene);
  }

  void TOPPASScene::paste(QPointF pos)
  {
    emit requestClipboardContent();

    if (clipboard_ != nullptr)
    {
      include(clipboard_, pos);
    }
  }

  void TOPPASScene::setClipboard(TOPPASScene* clipboard)
  {
    clipboard_ = clipboard;
  }

  void TOPPASScene::removeSelected()
  {
    abortPipeline();
    QList<TOPPASVertex*> vertices_to_be_removed;
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      if ((*it)->isSelected())
      {
        // also select all in and out edges (will be deleted below)
        for (TOPPASVertex::ConstEdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
        {
          (*e_it)->setSelected(true);
        }
        for (TOPPASVertex::ConstEdgeIterator e_it = (*it)->outEdgesBegin(); e_it != (*it)->outEdgesEnd(); ++e_it)
        {
          (*e_it)->setSelected(true);
        }
        vertices_to_be_removed.push_back(*it);
      }
    }
    QList<TOPPASEdge*> edges_to_be_removed;
    for (EdgeIterator it = edgesBegin(); it != edgesEnd(); ++it)
    {
      if ((*it)->isSelected())
      {
        edges_to_be_removed.push_back(*it);
      }
    }

    for (TOPPASEdge* edge : edges_to_be_removed)
    {
      edges_.removeAll(edge);
      removeItem(edge); // remove from scene
      delete edge;
    }
    for (TOPPASVertex* vertex : vertices_to_be_removed)
    {
      execution_->previous_result.nodes.erase(execution_->ids.at(vertex));
      execution_->ids.erase(vertex);
      vertices_.removeAll(vertex);
      removeItem(vertex); // remove from scene
      delete vertex;
    }

    execution_->discardUnusedTemporaryRuns();

    topoSort();
    updateEdgeColors();
    setChanged(true);
  }

  bool TOPPASScene::isEdgeAllowed_(TOPPASVertex* u, TOPPASVertex* v)
  {
    if (u == nullptr || v == nullptr || u == v ||
        // edges leading to input files make no sense:
        qobject_cast<TOPPASInputFileListVertex*>(v) ||
        // neither do edges coming from output files:
        qobject_cast<TOPPASOutputFileListVertex*>(u) ||
        // or edges coming from output directories:
        qobject_cast<TOPPASOutputFolderVertex*>(u) ||
        // nor edges from input/merger/splitter directly to output:
        ((qobject_cast<TOPPASInputFileListVertex*>(u) || qobject_cast<TOPPASMergerVertex*>(u) || qobject_cast<TOPPASSplitterVertex*>(u))
           && (qobject_cast<TOPPASOutputFileListVertex*>(v) || qobject_cast<TOPPASOutputFolderVertex*>(v)))
        ||
        // nor multiple incoming edges for an output or splitter node:
        ((qobject_cast<TOPPASOutputFileListVertex*>(v) || qobject_cast<TOPPASOutputFolderVertex*>(v) || qobject_cast<TOPPASSplitterVertex*>(v))
         && (v->inEdgesBegin() != v->inEdgesEnd())))
    {
      return false;
    }

    // can't have more incoming edges than a tool has inputs:
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(v);
    if (tv)
    {
      QVector<TOPPASToolVertex::IOInfo> input_infos = tv->getInputParameters();
      if (tv->incomingEdgesCount() >= Size(input_infos.size()))
      {
        return false;
      }
      // also, no edges from collectors to tools without input file lists:
      // @TODO: what if the input file list is already occupied by an edge?
      TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(u);
      if (mv && !mv->roundBasedMode())
      {      
        bool any_list = TOPPASToolVertex::IOInfo::isAnyList(input_infos);
        if (!any_list)
        {
          return false;
        }
      }
    }
    // no edges to splitters from tools without output file lists:
    if (qobject_cast<TOPPASSplitterVertex*>(v))
    {
      TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(u);
      if (tv)
      {
        QVector<TOPPASToolVertex::IOInfo> output_infos = tv->getOutputParameters();
        bool any_list = TOPPASToolVertex::IOInfo::isAnyList(output_infos);
        if (!any_list)
        {
          return false;
        }
      }
    }

    // does this edge already exist?
    for (TOPPASVertex::ConstEdgeIterator it = u->outEdgesBegin(); it != u->outEdgesEnd(); ++it)
    {
      if ((*it)->getTargetVertex() == v)
      {
        return false;
      }
    }

    // insert edge between u and v for testing, is removed afterwards
    TOPPASEdge* test_edge = new TOPPASEdge(u, QPointF());
    test_edge->setTargetVertex(v);
    u->addOutEdge(test_edge);
    v->addInEdge(test_edge);
    addEdge(test_edge);

    bool graph_has_cycles = false;
    // find back edges via DFS
    for (TOPPASVertex* vertex : vertices_)
    {
      vertex->setDFSColor(TOPPASVertex::DFS_WHITE);
    }
    for (TOPPASVertex* vertex : vertices_)
    {
      if (vertex->getDFSColor() == TOPPASVertex::DFS_WHITE)
      {
        graph_has_cycles = dfsVisit_(vertex);
        if (graph_has_cycles)
        {
          break;
        }
      }
    }

    // remove previously inserted edge
    edges_.removeAll(test_edge);
    removeItem(test_edge);
    delete test_edge;

    return !graph_has_cycles;
  }

  void TOPPASScene::updateEdgeColors()
  {
    for (TOPPASEdge* edge : edges_)
    {
      edge->updateColor();
    }
    update(sceneRect());
  }

  bool TOPPASScene::dfsVisit_(TOPPASVertex* vertex)
  {
    vertex->setDFSColor(TOPPASVertex::DFS_GRAY);
    for (TOPPASVertex::ConstEdgeIterator it = vertex->outEdgesBegin(); it != vertex->outEdgesEnd(); ++it)
    {
      TOPPASVertex* target = (*it)->getTargetVertex();
      if (target->getDFSColor() == TOPPASVertex::DFS_WHITE)
      {
        if (dfsVisit_(target))
        {
          // back edge found
          return true;
        }
      }
      else if (target->getDFSColor() == TOPPASVertex::DFS_GRAY)
      {
        // back edge found
        return true;
      }
    }
    vertex->setDFSColor(TOPPASVertex::DFS_BLACK);
    return false;
  }

  void TOPPASScene::resetDownstream(TOPPASVertex* vertex)
  {
    abortPipeline();
    std::set<TOPPASVertex*> visited;
    std::vector<TOPPASVertex*> pending {vertex};
    while (! pending.empty())
    {
      auto* current = pending.back();
      pending.pop_back();
      if (! current || ! visited.insert(current).second) { continue; }
      if (auto it = execution_->ids.find(current); it != execution_->ids.end()) { execution_->previous_result.nodes.erase(it->second); }
      current->reset(false);
      for (auto edge = current->outEdgesBegin(); edge != current->outEdgesEnd(); ++edge)
      {
        pending.push_back((*edge)->getTargetVertex());
      }
    }
    execution_->discardUnusedTemporaryRuns();
  }

  PipelineGraph TOPPASScene::pipelineGraph_() const
  {
    PipelineGraph graph;
    graph.version = VersionInfo::getVersion();
    graph.filename = file_name_;
    graph.description = fromQString(description_text_);
    for (auto* vertex : vertices_)
    {
      PipelineGraph::Node node;
      node.id = execution_->ids.at(vertex);
      node.topo_number = vertex->getTopoNr();
      node.x = vertex->x();
      node.y = vertex->y();
      node.recycle_output = vertex->isRecyclingEnabled();
      if (auto* input = qobject_cast<TOPPASInputFileListVertex*>(vertex))
      {
        node.kind = PipelineGraph::Kind::INPUT;
        node.resource_key = fromQString(input->getKey());
        for (const auto& file : input->getFileNames())
        {
          node.files.push_back(fromQString(file));
        }
      }
      else if (auto* tool = qobject_cast<TOPPASToolVertex*>(vertex))
      {
        node.kind = PipelineGraph::Kind::TOOL;
        node.tool_name = tool->getName();
        node.tool_type = tool->getType();
        node.parameters = tool->getParam();
      }
      else if (auto* merger = qobject_cast<TOPPASMergerVertex*>(vertex))
      {
        node.kind = PipelineGraph::Kind::MERGER;
        node.round_based = merger->roundBasedMode();
      }
      else if (qobject_cast<TOPPASSplitterVertex*>(vertex)) { node.kind = PipelineGraph::Kind::SPLITTER; }
      else if (auto* output = qobject_cast<TOPPASOutputVertex*>(vertex))
      {
        node.kind = qobject_cast<TOPPASOutputFolderVertex*>(vertex) ? PipelineGraph::Kind::OUTPUT_DIRECTORY : PipelineGraph::Kind::OUTPUT;
        node.output_folder = fromQString(output->getOutputFolderName());
      }
      graph.nodes.push_back(std::move(node));
    }
    for (auto* edge : edges_)
    {
      if (! edge->getSourceVertex() || ! edge->getTargetVertex()) { continue; }
      PipelineGraph::Edge binding;
      binding.source = execution_->ids.at(edge->getSourceVertex());
      binding.target = execution_->ids.at(edge->getTargetVertex());
      auto portName = [](TOPPASVertex* vertex, Int index, bool input) {
        auto* tool = qobject_cast<TOPPASToolVertex*>(vertex);
        if (! tool || index < 0) { return std::string {}; }
        const auto ports = input ? tool->getInputParameters() : tool->getOutputParameters();
        if (index >= ports.size()) { return std::string {}; }
        return ports[index].param_name;
      };
      binding.source_port = portName(edge->getSourceVertex(), edge->getSourceOutParam(), false);
      binding.target_port = portName(edge->getTargetVertex(), edge->getTargetInParam(), true);
      graph.edges.push_back(std::move(binding));
    }
    return graph;
  }

  void TOPPASScene::runPipeline()
  { startExecution_(nullptr); }

  void TOPPASScene::resumePipeline(TOPPASToolVertex* vertex)
  { startExecution_(vertex); }

  void TOPPASScene::startExecution_(TOPPASToolVertex* resume_vertex)
  {
    abortPipeline();
    if (! sanityCheck_(gui_) || ! askForOutputDir(resume_vertex == nullptr)) { return; }
    auto graph = pipelineGraph_();
    std::optional<Size> resume;
    if (resume_vertex) { resume = execution_->ids.at(resume_vertex); }
    if (! resume)
    {
      execution_->previous_result = {};
      execution_->discardUnusedTemporaryRuns();
      for (auto* vertex : vertices_)
      {
        vertex->reset(false);
      }
    }
    else
    {
      resetDownstream(resume_vertex);
    }
    auto run_directory = std::make_unique<QTemporaryDir>(tmp_path_ + "/TOPPAS_run_XXXXXX");
    if (! run_directory->isValid())
    {
      error_occured_ = true;
      emit messageReady("Could not create a temporary workflow directory in '" + tmp_path_ + "'.");
      emit pipelineExecutionFailed(TOPPBase::CANNOT_WRITE_OUTPUT_FILE);
      return;
    }
    PipelineExecutor::Options options;
    options.output_directory = fromQString(out_dir_);
    options.temp_directory = fromQString(run_directory->path());
    execution_->temporary_runs.push_back(std::move(run_directory));
    options.num_jobs = static_cast<Size>(allowed_threads_);
    options.keep_temporary_files = true;
    const auto previous = execution_->previous_result;
    execution_->worker_result = {};
    execution_->executor = std::make_unique<PipelineExecutor>();
    auto* executor = execution_->executor.get();
    const Size generation = ++execution_->generation;
    error_occured_ = false;
    setPipelineRunning(true);
    // The worker owns a graph snapshot and sends value events only. The Qt context
    // drops queued calls after destruction; generation also rejects calls after edits.
    execution_->worker = std::thread([this, executor, graph = std::move(graph), options, previous, resume, generation]() {
      auto callback = [this, generation](const PipelineExecutor::Event& event) {
        QMetaObject::invokeMethod(
          this,
          [this, generation, event]() {
            if (execution_->generation != generation) { return; }
            if (event.type == PipelineExecutor::Event::Type::LOG || event.type == PipelineExecutor::Event::Type::OUTPUT_WRITTEN)
            {
              emit messageReady(toQString(event.text));
            }
            QString node_label;
            for (const auto& [vertex, id] : execution_->ids)
            {
              if (id != event.node_id) { continue; }
              node_label = toQString(vertex->getName()) + " (#" + QString::number(vertex->getTopoNr()) + ")";
              if (event.type == PipelineExecutor::Event::Type::NODE_SCHEDULED || event.type == PipelineExecutor::Event::Type::NODE_STARTED
                  || event.type == PipelineExecutor::Event::Type::ROUND_COMPLETED || event.type == PipelineExecutor::Event::Type::NODE_FINISHED
                  || event.type == PipelineExecutor::Event::Type::NODE_FAILED)
              {
                presentNodeResult(vertex, event.result);
                if (event.type == PipelineExecutor::Event::Type::NODE_SCHEDULED)
                {
                  if (auto* tool = qobject_cast<TOPPASToolVertex*>(vertex)) { tool->toolScheduledSlot(); }
                }
              }
              break;
            }
            if (event.type == PipelineExecutor::Event::Type::NODE_STARTED || event.type == PipelineExecutor::Event::Type::NODE_FINISHED
                || event.type == PipelineExecutor::Event::Type::NODE_FAILED)
            {
              const QString state = event.type == PipelineExecutor::Event::Type::NODE_STARTED    ? " started."
                                    : event.type == PipelineExecutor::Event::Type::NODE_FINISHED ? " finished."
                                                                                                 : " failed.";
              emit messageReady(node_label + state);
            }
            update(sceneRect());
          },
          Qt::QueuedConnection);
      };
      try
      {
        execution_->worker_result = executor->run(graph, options, callback, resume ? &previous : nullptr, resume);
      }
      catch (const std::exception& error)
      {
        execution_->worker_result.exit_code = 1;
        execution_->worker_result.error_message = error.what();
      }
      catch (...)
      {
        execution_->worker_result.exit_code = 1;
        execution_->worker_result.error_message = "Unexpected error while executing the workflow.";
      }
      QMetaObject::invokeMethod(
        this,
        [this, generation]() {
          if (execution_->generation != generation) { return; }
          execution_->worker.join();
          execution_->acceptResult();
          execution_->executor.reset();
          for (const auto& [vertex, id] : execution_->ids)
          {
            auto result = execution_->previous_result.nodes.find(id);
            if (result != execution_->previous_result.nodes.end()) { presentNodeResult(vertex, result->second); }
          }
          setPipelineRunning(false);
          const auto& result = execution_->previous_result;
          if (result.exit_code != 0)
          {
            error_occured_ = true;
            emit messageReady(toQString(result.error_message));
            emit pipelineExecutionFailed(result.exit_code);
          }
          else
            emit entirePipelineFinished();
        },
        Qt::QueuedConnection);
    });
  }

  bool TOPPASScene::store(const std::string& file)
  {
    for (auto* edge : edges_)
    {
      if (edge->getEdgeStatus() != TOPPASEdge::ES_VALID && edge->getEdgeStatus() != TOPPASEdge::ES_NOT_READY_YET) { return false; }
    }
    try
    {
      PipelineFile().store(file, pipelineGraph_());
    }
    catch (const std::exception& error)
    {
      emit messageReady(toQString(error.what()));
      return false;
    }
    setChanged(false);
    file_name_ = file;
    return true;
  }

  QString TOPPASScene::getDescription() const
  {
    return description_text_;
  }

  ///
  void TOPPASScene::setDescription(const QString& desc)
  {
    description_text_ = desc;
  }

  void TOPPASScene::load(const std::string& file)
  {
    file_name_ = file;

    if (File::empty(file)) // allow opening of 0-byte files as pretend they are empty, new TOPPAS files
    {
      return;
    }

    Param load_param;
    ParamXMLFile paramFile;
    paramFile.load(file, load_param);

    // check for TOPPAS file version. Deny loading if too old or too new
    // get version of TOPPAS file
    std::string file_version = "1.8.0"; // default (were we did not have the tag)
    if (load_param.exists("info:version"))
    {
      file_version = load_param.getValue("info:version").toString();
    }
    VersionInfo::VersionDetails v_file = VersionInfo::VersionDetails::create(file_version);
    VersionInfo::VersionDetails v_this_low = VersionInfo::VersionDetails::create("1.9.0"); // last compatible TOPPAS file version
    VersionInfo::VersionDetails v_this_high = VersionInfo::VersionDetails::create(VersionInfo::getVersion()); // last compatible TOPPAS file version
    if (v_file < v_this_low)
    {
      if (!this->gui_)
      {
        std::cerr << "The TOPPAS file is too old! Please update the file using TOPPAS or INIUpdater!" << std::endl;
      }
      else if (this->gui_)
      {
        if (QMessageBox::warning(nullptr, tr("Old TOPPAS file -- convert and override?"), tr("The TOPPAS file you downloaded was created with an old incompatible version of TOPPAS.\nShall we try to convert the file?! The original file will be overridden, but a backup file will be saved in the same directory.\n"), QMessageBox::Yes, QMessageBox::No) == QMessageBox::No)
        {
          return;
        }
        // only update in GUI mode, as in non-GUI mode, we'd create infinite recursive calls when instantiating TOPPASScene in INIUpdater
#ifdef OPENMS_WINDOWSPLATFORM
        std::string extra_quotes = "\""; // note: double quoting required for Windows, as outer quotes are required by cmd.exe (arghh)...
#else
        std::string extra_quotes;
#endif

        std::string cmd = extra_quotes + "\"" + File::findSiblingTOPPExecutable("INIUpdater") + "\" -in \"" + file + "\" -i " + extra_quotes;
        std::cerr << cmd << "\n\n";
        if (std::system(cmd.c_str()))
        {
          QMessageBox::warning(nullptr, tr("INIUpdater failed"), tr("Updating using the INIUpdater tool failed. Please submit a bug report!\n"), QMessageBox::Ok);
          return;
        }
        // reload updated file
        ParamXMLFile paramFile;
        paramFile.load(file, load_param);
      }
    }
    else if (v_file > v_this_high)
    {
      if (this->gui_ && QMessageBox::warning(nullptr, tr("TOPPAS file too new"), tr("The TOPPAS file you downloaded was created with a more recent version of TOPPAS. Shall we will try to open it?\nIf this fails, update to the new TOPPAS version.\n"), QMessageBox::Yes, QMessageBox::No) == QMessageBox::No)
      {
        return;
      }
    }


    PipelineGraph graph;
    PipelineFile().loadParam(load_param, graph, file);
    abortPipeline();
    for (auto* vertex : vertices_)
    {
      vertex->blockSignals(true);
      vertex->setSelected(true);
    }
    for (auto* edge : edges_)
    {
      edge->blockSignals(true);
      edge->setSelected(true);
    }
    removeSelected();
    execution_->previous_result = {};
    execution_->discardUnusedTemporaryRuns();
    description_text_ = toQString(graph.description);
    std::map<Size, TOPPASVertex*> by_id;
    for (const auto& node : graph.nodes)
    {
      TOPPASVertex* vertex = nullptr;
      switch (node.kind)
      {
        case PipelineGraph::Kind::INPUT: {
          QStringList files;
          for (const auto& path : node.files)
          {
            files.push_back(toQString(path));
          }
          auto* input = new TOPPASInputFileListVertex(files);
          input->setKey(toQString(node.resource_key));
          vertex = input;
          break;
        }
        case PipelineGraph::Kind::TOOL: {
          auto tool = std::make_unique<TOPPASToolVertex>(node.tool_name, node.tool_type);
          tool->setParam(node.parameters);
          connectToolVertexSignals(tool.get());
          vertex = tool.release();
          break;
        }
        case PipelineGraph::Kind::MERGER: {
          auto* merger = new TOPPASMergerVertex(node.round_based);
          connectMergerVertexSignals(merger);
          vertex = merger;
          break;
        }
        case PipelineGraph::Kind::SPLITTER:
          vertex = new TOPPASSplitterVertex();
          break;
        case PipelineGraph::Kind::OUTPUT:
        case PipelineGraph::Kind::OUTPUT_DIRECTORY: {
          TOPPASOutputVertex* output = node.kind == PipelineGraph::Kind::OUTPUT ? static_cast<TOPPASOutputVertex*>(new TOPPASOutputFileListVertex())
                                                                                : static_cast<TOPPASOutputVertex*>(new TOPPASOutputFolderVertex());
          output->setOutputFolderName(toQString(node.output_folder));
          connectOutputVertexSignals(output);
          vertex = output;
          break;
        }
      }
      vertex->blockSignals(true);
      vertex->setPos(node.x, node.y);
      vertex->setRecycling(node.recycle_output);
      vertex->setTopoNr(static_cast<UInt>(node.topo_number));
      addVertex(vertex);
      execution_->ids[vertex] = node.id;
      execution_->next_id = std::max(execution_->next_id, node.id + 1);
      by_id[node.id] = vertex;
      connectVertexSignals(vertex);
    }
    auto portIndex = [&graph, this](TOPPASVertex* vertex, const std::string& name, bool input) {
      if (name.empty()) { return Int {-1}; }
      if (graph.legacy_port_indices) { return StringUtils::toInt32(name); }
      auto* tool = qobject_cast<TOPPASToolVertex*>(vertex);
      if (! tool) { return Int {-1}; }
      const auto ports = input ? tool->getInputParameters() : tool->getOutputParameters();
      for (int i = 0; i < ports.size(); ++i)
      {
        if (ports[i].param_name == name) return i;
      }
      logTOPPOutput(toQString("Could not find parameter '" + name + "'. Check edge!"));
      return Int {-1};
    };
    for (const auto& binding : graph.edges)
    {
      auto* source = by_id.at(binding.source);
      auto* target = by_id.at(binding.target);
      auto* edge = new TOPPASEdge();
      edge->setSourceVertex(source);
      edge->setTargetVertex(target);
      edge->setSourceOutParam(portIndex(source, binding.source_port, false));
      edge->setTargetInParam(portIndex(target, binding.target_port, true));
      source->addOutEdge(edge);
      target->addInEdge(edge);
      connectEdgeSignals(edge);
      addEdge(edge);
    }
    for (auto* vertex : vertices_)
    {
      vertex->blockSignals(false);
    }
    updateEdgeColors();
    setChanged(false);
  }


  void TOPPASScene::include(TOPPASScene* tmp_scene, QPointF pos)
  {
    qreal x_offset, y_offset;
    if (pos == QPointF()) // pasted via Ctrl-V (no mouse position given)
    {
      x_offset = 30.0; // move just a tad (in relation to old content)
      y_offset = 30.0;
    }
    else
    {
      QRectF new_bounding_rect = tmp_scene->itemsBoundingRect();
      x_offset = pos.x() - new_bounding_rect.left();
      y_offset = pos.y() - new_bounding_rect.top();
    }
    std::map<TOPPASVertex*, TOPPASVertex*> vertex_map;

    for (VertexIterator it = tmp_scene->verticesBegin(); it != tmp_scene->verticesEnd(); ++it)
    {
      TOPPASVertex* v = *it;
      TOPPASVertex* new_v = nullptr;

      TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(v);
      if (iflv)
      {
        TOPPASInputFileListVertex* new_iflv = new TOPPASInputFileListVertex(*iflv);
        std::set<QString> keys;
        for (auto* existing : vertices_)
        {
          if (auto* input = qobject_cast<TOPPASInputFileListVertex*>(existing)) { keys.insert(input->getKey()); }
        }
        const QString original_key = new_iflv->getKey();
        QString key = original_key == QString::number(iflv->getTopoNr()) ? QString {} : original_key;
        Size suffix = 2;
        while (! key.isEmpty() && keys.count(key))
        {
          key = original_key + "_" + QString::number(suffix++);
        }
        new_iflv->setKey(key);
        new_v = new_iflv;
      }

      TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(v);
      if (oflv)
      {
        TOPPASOutputFileListVertex* new_oflv = new TOPPASOutputFileListVertex(*oflv);
        new_v = new_oflv;

        connectOutputVertexSignals(new_oflv);
      }
      TOPPASOutputFolderVertex* ofv = qobject_cast<TOPPASOutputFolderVertex*>(v);
      if (ofv)
      {
        TOPPASOutputFolderVertex* new_ofv = new TOPPASOutputFolderVertex(*ofv);
        new_v = new_ofv;

        connectOutputVertexSignals(new_ofv);
      }

      TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(v);
      if (tv)
      {
        TOPPASToolVertex* new_tv = new TOPPASToolVertex(*tv);
        new_v = new_tv;

        connectToolVertexSignals(new_tv);
      }

      TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(v);
      if (mv)
      {
        TOPPASMergerVertex* new_mv = new TOPPASMergerVertex(*mv);
        new_v = new_mv;

        connectMergerVertexSignals(new_mv);
      }

      TOPPASSplitterVertex* sv = qobject_cast<TOPPASSplitterVertex*>(v);
      if (sv)
      {
        TOPPASSplitterVertex* new_sv = new TOPPASSplitterVertex(*sv);
        new_v = new_sv;
      }

      if (!new_v)
      {
        std::cerr << "Unknown vertex type! Aborting." << std::endl;
        return;
      }

      vertex_map[v] = new_v;
      new_v->moveBy(x_offset, y_offset);
      connectVertexSignals(new_v);
      addVertex(new_v);

      // temporarily block signals in order that the first topo sort does not set the changed flag
      new_v->blockSignals(true);
    }

    // add all edges (are not copied by copy constructors of vertices)
    for (EdgeIterator it = tmp_scene->edgesBegin(); it != tmp_scene->edgesEnd(); ++it)
    {
      TOPPASVertex* old_source = (*it)->getSourceVertex();
      TOPPASVertex* old_target = (*it)->getTargetVertex();
      TOPPASVertex* new_source = vertex_map[old_source];
      TOPPASVertex* new_target = vertex_map[old_target];
      TOPPASEdge* new_e = new TOPPASEdge();
      new_e->setSourceVertex(new_source);
      new_e->setTargetVertex(new_target);
      new_e->setSourceOutParam((*it)->getSourceOutParam());
      new_e->setTargetInParam((*it)->getTargetInParam());
      new_source->addOutEdge(new_e);
      new_target->addInEdge(new_e);

      connectEdgeSignals(new_e);

      addEdge(new_e);
    }

    // select new items (so the user can move them); edges do not need to be selected, only vertices
    unselectAll();
    for (std::map<TOPPASVertex*, TOPPASVertex*>::iterator it = vertex_map.begin(); it != vertex_map.end(); ++it)
    {
      it->second->setSelected(true);
    }

    topoSort();
    // unblock signals again
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      (*it)->blockSignals(false);
    }

    updateEdgeColors();
  }

  const std::string& TOPPASScene::getSaveFileName()
  {
    return file_name_;
  }

  void TOPPASScene::setSaveFileName(const std::string& name)
  {
    file_name_ = name;
  }

  void TOPPASScene::unselectAll()
  {
    const QList<QGraphicsItem*>& all_items = items();
    for (QGraphicsItem * item : all_items)
    {
      item->setSelected(false);
    }
    update(sceneRect());
  }


  void TOPPASScene::pipelineErrorSlot(int return_code, const QString& msg)
  {
    logTOPPOutput(msg); // print to log window or console
    error_occured_ = true;
    setPipelineRunning(false);
    abortPipeline();
    emit pipelineExecutionFailed(return_code);
  }

  void TOPPASScene::writeToLogFile_(const QString& text)
  {
    QFile logfile(out_dir_ + QDir::separator() + "TOPPAS.log");
    if (!logfile.open(QIODevice::Append | QIODevice::Text))
    {
      std::cerr << "Could not write to logfile '" << fromQString(logfile.fileName()) << "'" << std::endl;
      return;
    }

    QTextStream ts(&logfile);
    ts << "\n" << text << "\n";
    logfile.close();
  }

  void TOPPASScene::logTOPPOutput(const QString& out)
  {
    TOPPASToolVertex* sender = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (!sender)
    {
      //return;
    }
    std::string text = fromQString(out);

    if (!gui_)
    {
      std::cout << std::endl << text << std::endl;
    }
    emit messageReady(out); // let TOPPAS know about it

    writeToLogFile_(toQString(text));
  }

  void TOPPASScene::logToolStarted()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      std::string text = tv->getName();
      std::string type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " started. Processing ...";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(toQString(text));
    }
  }

  void TOPPASScene::logToolFinished()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      std::string text = tv->getName();
      std::string type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " finished!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(toQString(text));
    }
  }

  void TOPPASScene::logToolFailed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      std::string text = tv->getName();
      std::string type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " failed!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(toQString(text));
    }
  }

  void TOPPASScene::logToolCrashed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      std::string text = tv->getName();
      std::string type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " crashed!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(toQString(text));
    }
  }

  void TOPPASScene::logOutputFileWritten(const std::string& file)
  {
    std::string text = "Output file '" + file + "' written.";

    if (!gui_)
    {
      std::cout << std::endl << text << std::endl;
    }

    writeToLogFile_(toQString(text));
  }

  void TOPPASScene::topoSort(bool /*resort_all*/)
  {
    // Existing vertices are kept in their previous topological order. The shared
    // stable scan therefore preserves their numbers when appending new vertices,
    // including the former resort_all=false case.
    auto graph = pipelineGraph_();
    graph.assignTopologicalNumbers();
    for (auto* vertex : vertices_)
    {
      const auto& node = graph.node(execution_->ids.at(vertex));
      if (auto* input = qobject_cast<TOPPASInputFileListVertex*>(vertex)) { input->setKey(toQString(node.resource_key)); }
      vertex->setTopoNr(static_cast<UInt>(node.topo_number));
      vertex->setTopoSortMarked(true);
    }
    std::sort(vertices_.begin(), vertices_.end(), [](TOPPASVertex* a, TOPPASVertex* b) { return a->getTopoNr() < b->getTopoNr(); });
    update(sceneRect());
  }

  const QString& TOPPASScene::getOutDir() const
  {
    return out_dir_;
  }

  const QString& TOPPASScene::getTempDir() const
  {
    return tmp_path_;
  }

  void TOPPASScene::setOutDir(const QString& dir)
  {
    QDir d(dir);
    out_dir_ = d.absolutePath();
    user_specified_out_dir_ = true;
  }

  void TOPPASScene::moveSelectedItems(qreal dx, qreal dy)
  {
    setActionMode(AM_MOVE);

    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      if (!(*it)->isSelected())
      {
        continue;
      }
      for (TOPPASVertex::ConstEdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
      {
        (*e_it)->prepareResize();
      }
      for (TOPPASVertex::ConstEdgeIterator e_it = (*it)->outEdgesBegin(); e_it != (*it)->outEdgesEnd(); ++e_it)
      {
        (*e_it)->prepareResize();
      }

      (*it)->moveBy(dx, dy);
    }

    setChanged(true);
  }

  void TOPPASScene::snapToGrid()
  {
    int grid_step = 20;

    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      //only make selected nodes snap (those might have been moved)
      if (!(*it)->isSelected())
      {
        continue;
      }

      int x_int = (int)((*it)->x());
      int y_int = (int)((*it)->y());
      int prev_grid_x = x_int - (x_int % grid_step);
      int prev_grid_y = y_int - (y_int % grid_step);
      int new_x = prev_grid_x;
      int new_y = prev_grid_y;

      if (x_int - prev_grid_x > (grid_step / 2))
      {
        new_x += grid_step;
      }
      if (y_int - prev_grid_y > (grid_step / 2))
      {
        new_y += grid_step;
      }

      (*it)->setPos(QPointF(new_x, new_y));
    }

    update(sceneRect());
  }

  bool TOPPASScene::saveIfChanged()
  {
    // Save changes
    if (gui_ && changed_)
    {
      QString name = file_name_.empty() ? "Untitled" : toQString(File::basename(file_name_));
      QMessageBox::StandardButton ret;
      ret = QMessageBox::warning(views().first(), "Save changes?", "'" + name + "' has been modified.\n\nDo you want to save your changes?", QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
      if (ret == QMessageBox::Save)
      {
        emit saveMe();
        if (changed_)
        {
          //user has not saved the file (aborted save dialog)
          return false;
        }
      }
      else if (ret == QMessageBox::Cancel)
      {
        return false;
      }
    }
    return true;
  }

  void TOPPASScene::setChanged(bool b)
  {
    if (changed_ != b)
    {
      changed_ = b;
      emit mainWindowNeedsUpdate();
    }
  }

  bool TOPPASScene::wasChanged() const
  {
    return changed_;
  }

  bool TOPPASScene::isPipelineRunning() const
  {
    return running_;
  }

  void TOPPASScene::abortPipeline()
  {
    ++execution_->generation;
    if (execution_->executor) { execution_->executor->cancel(); }
    if (execution_->worker.joinable())
    {
      execution_->worker.join();
      execution_->acceptResult();
      for (const auto& [vertex, id] : execution_->ids)
      {
        const auto result = execution_->previous_result.nodes.find(id);
        if (result != execution_->previous_result.nodes.end()) { presentNodeResult(vertex, result->second); }
      }
    }
    execution_->executor.reset();
    execution_->discardUnusedTemporaryRuns();
    if (running_) { setPipelineRunning(false); }
  }


  void TOPPASScene::setPipelineRunning(bool b)
  {
    running_ = b;
    emit mainWindowNeedsUpdate();
    if (!running_) // whenever we stop the pipeline and user is not looking, the icon should flash
    {
      QApplication::alert(nullptr); // flash Taskbar || Dock
    }
  }


  bool TOPPASScene::askForOutputDir(bool always_ask)
  {
    if (gui_)
    {
      if (always_ask || !user_specified_out_dir_)
      {
        TOPPASOutputFilesDialog tofd(out_dir_, allowed_threads_);
        if (tofd.exec())
        {
          setOutDir(tofd.getDirectory());
          setAllowedThreads(tofd.getNumJobs());
        }
        else
        {
          return false;
        }
      }
    }

    return true;
  }

  void TOPPASScene::contextMenuEvent(QGraphicsSceneContextMenuEvent* event)
  {
    QPointF scene_pos = event->scenePos();
    QGraphicsItem* clicked_item = itemAt(scene_pos, QTransform());
    QMenu menu;

    if (clicked_item == nullptr)
    {
      QAction* new_action = menu.addAction("Paste");
      emit requestClipboardContent();
      if (clipboard_ == nullptr)
      {
        new_action->setEnabled(false);
      }
    }
    else
    {
      if (!clicked_item->isSelected())
      {
        unselectAll();
      }

      clicked_item->setSelected(true);

      // check which kinds of items are selected and display a context menu containing only actions compatible with all of them
      bool found_tool = false;
      bool found_input = false;
      bool found_output = false, found_output_files = false;
      bool found_merger = false;
      bool found_splitter = false;
      bool found_edge = false;
      bool disable_resume = this->isPipelineRunning();
      //bool disable_toppview = true;

      for (TOPPASEdge* edge : edges_)
      {
        if (edge->isSelected())
        {
          found_edge = true;
          break;
        }
      }

      for (TOPPASVertex* tv : vertices_)
      {
        if (!tv->isSelected())
        {
          continue;
        }

        if (qobject_cast<TOPPASToolVertex*>(tv))
        {
          found_tool = true;
          // all predecessor nodes finished successfully? if not, disable resuming
          for (ConstEdgeIterator it = tv->inEdgesBegin(); it != tv->inEdgesEnd(); ++it)
          {
            TOPPASToolVertex* pred_ttv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
            if (pred_ttv && (pred_ttv->getStatus() != TOPPASToolVertex::TOOL_SUCCESS))
            {
              disable_resume = true;
              break;
            }
          }
          continue;
        }
        if (qobject_cast<TOPPASInputFileListVertex*>(tv))
        {
          found_input = true;
          continue;
        }
        if (qobject_cast<TOPPASOutputVertex*>(tv))
        {
          found_output = true;
          // no continue here; derived classes below
        }
        if (qobject_cast<TOPPASOutputFileListVertex*>(tv))
        {
          found_output_files = true;
          continue;
        }
        if (qobject_cast<TOPPASMergerVertex*>(tv))
        {
          found_merger = true;
          continue;
        }
        if (qobject_cast<TOPPASSplitterVertex*>(tv))
        {
          found_splitter = true;
          continue;
        }
      }

      QSet<QString> action;

      if (found_tool)
      {
        action.insert("Edit parameters");
        action.insert("Resume");
        action.insert("Open files in TOPPView");
        action.insert("Open containing folder");
        //action.insert("Toggle breakpoint");
      }

      if (found_input)
      {
        action.insert("Change name");
        action.insert("Change files");
        action.insert("Open files in TOPPView");
        action.insert("Open containing folder");
      }

      if (found_output)
      {
        action.insert("Set output folder name");
        action.insert("Open containing folder");
      }
      if (found_output_files)
      {
        action.insert("Open files in TOPPView");
      }

      if (found_edge)
      {
        action.insert("Edit I/O mapping");
      }

      if (found_input || found_tool || found_merger || found_splitter)
      {
        action.insert("Toggle recycling mode");
      }

      QList<QSet<QString> > all_actions;
      all_actions.push_back(action);

      QSet<QString> supported_actions_set = all_actions.first();
      for (const QSet<QString>&action_set : all_actions)
      {
        supported_actions_set.intersect(action_set);
      }

      QList<QString> supported_actions = supported_actions_set.values();
      supported_actions << "Copy" << "Cut" << "Remove";
      for (const QString &supported_action : supported_actions)
      {
        QAction* new_action = menu.addAction(supported_action);
        if (supported_action == "Resume" && disable_resume)
        {
          new_action->setEnabled(false);
        }
      }
    }

    // ------ execute action  ------

    QAction* selected_action = menu.exec(event->screenPos());
    if (selected_action)
    {
      QString text = selected_action->text();

      if (text == "Remove")
      {
        removeSelected();
        event->accept();
        return;
      }

      if (text == "Copy")
      {
        copySelected();
        event->accept();
        return;
      }

      if (text == "Cut")
      {
        copySelected();
        removeSelected();
        event->accept();
        return;
      }

      if (text == "Paste")
      {
        paste(event->scenePos());
        event->accept();
        return;
      }

      for (QGraphicsItem* gi : selectedItems())
      {

        if (text == "Toggle recycling mode")
        {
          if (auto* tv = dynamic_cast<TOPPASVertex*>(gi); tv)
          {
            tv->invertRecylingMode();
            tv->update(tv->boundingRect());
          }
          continue;
        }

        if (auto* edge = dynamic_cast<TOPPASEdge*>(gi); edge)
        {
          if (text == "Edit I/O mapping")
          {
            edge->showIOMappingDialog();
          }

          continue;
        }

        if (auto* ttv = dynamic_cast<TOPPASToolVertex*>(gi); ttv)
        {
          if (text == "Edit parameters")
          {
            ttv->editParam();
          }
          else if (text == "Resume")
          {
            if (askForOutputDir(false)) { resumePipeline(ttv); }
          }
          else if (text == "Toggle breakpoint")
          {
            ttv->toggleBreakpoint();
            ttv->update(ttv->boundingRect());
          }
          else if (text == "Open files in TOPPView")
          {
            QStringList all_out_files = ttv->getFileNames();
            emit openInTOPPView(all_out_files);
          }
          else if (text == "Open containing folder")
          {
            ttv->openContainingFolder();
          }

          continue;
        }

        if (auto* ifv = dynamic_cast<TOPPASInputFileListVertex*>(gi); ifv)
        {
          if (text == "Open files in TOPPView")
          {
            QStringList in_files = ifv->getFileNames();
            emit openInTOPPView(in_files);
          }
          else if (text == "Open containing folder")
          {
            ifv->openContainingFolder();
          }
          else if (text == "Change files")
          {
            ifv->showFilesDialog();
          }
          else if (text == "Change name")
          {
            TOPPASVertexNameDialog dlg(ifv->getKey());
            if (dlg.exec())
            {
              ifv->setKey(dlg.getName());
            }
          }
          continue;
        }

        
        if (auto* ov = dynamic_cast<TOPPASOutputVertex*>(gi); ov)
        {
          if (text == "Open containing folder")
          {
            ov->openContainingFolder();
          }
          else if (text == "Set output folder name")
          {
            TOPPASVertexNameDialog dlg(ov->getOutputFolderName(), "[a-zA-Z0-9_-]*");
            if (dlg.exec())
            { ov->setOutputFolderName(dlg.getName());
            }
          }
          // no continue - derived classes below
        }
        
        if (auto* ofv = dynamic_cast<TOPPASOutputFileListVertex*>(gi); ofv)
        {
          if (text == "Open files in TOPPView")
          {
            QStringList out_files = ofv->getFileNames();
            emit openInTOPPView(out_files);
          }
          continue;
        }
      }
    }

    event->accept();
  }


  bool TOPPASScene::sanityCheck_(bool allowUserOverride)
  {
    QStringList strange_vertices;

    // ----- are there any input nodes and are files specified? ----

    /// check if we have any input nodes
    QVector<TOPPASInputFileListVertex*> input_nodes;
    for (TOPPASVertex* tv : vertices_)
    {
      TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
      if (iflv)
      {
        input_nodes.push_back(iflv);
      }
    }
    if (input_nodes.empty())
    {
      if (allowUserOverride)
      {
        QMessageBox::warning(nullptr, "No input files", "The pipeline does not contain any input file nodes!");
      }
      else
      {
        std::cerr << "The pipeline does not contain any input file nodes!" << std::endl;
      }
      return false;
    }

    /// warn about empty input nodes
    for (TOPPASInputFileListVertex* iflv : input_nodes)
    {
      if ((iflv->outgoingEdgesCount() > 0) && (iflv->getFileNames().empty()))  // allow disconnected input node with empty file list
      {
        strange_vertices.push_back(QString::number(iflv->getTopoNr()));
      }
    }
    if (!strange_vertices.empty())
    {
      if (allowUserOverride)
      {
        QMessageBox::warning(views().first(), "Empty input file nodes",
                             QString("Node")
                             + (strange_vertices.size() > 1 ? "s " : " ")
                             + strange_vertices.join(", ")
                             + (strange_vertices.size() > 1 ? " have " : " has ")
                             + " an empty input file list!");
      }
      else
      {
        std::cerr << "Pipeline contains input file nodes without specified files!" << std::endl;
      }
      return false;
    }

    /// check if input files exist
    strange_vertices.clear();
    for (TOPPASInputFileListVertex* iflv : input_nodes)
    {
      if ((iflv->outgoingEdgesCount() > 0) && (!iflv->fileNamesValid()))  // allow disconnected input node with invalid files
      {
        strange_vertices.push_back(QString::number(iflv->getTopoNr()));
      }
    }
    if (!strange_vertices.empty())
    {
      if (allowUserOverride)
      {
        QMessageBox::warning(views().first(), "Input file names wrong",
                             QString("Node")
                             + (strange_vertices.size() > 1 ? "s " : " ")
                             + strange_vertices.join(", ")
                             + (strange_vertices.size() > 1 ? " have " : " has ")
                             + " invalid (non-existing or duplicate) input files!");
      }
      else
      {
        std::cerr << "Pipeline contains input file nodes with invalid (non-existing or duplicate) input files!" << std::endl;
      }
      return false;
    }

    // ----- are there nodes without parents (besides input nodes)? -----
    strange_vertices.clear();
    for (TOPPASVertex* tv : vertices_)
    {
      if (qobject_cast<TOPPASInputFileListVertex*>(tv)) // input nodes don't need a parent
      {
        continue;
      }
      if (tv->inEdgesBegin() == tv->inEdgesEnd())
      {
        strange_vertices << QString::number(tv->getTopoNr());
        tv->markUnreachable();
      }
    }
    if (!strange_vertices.empty())
    {
      if (allowUserOverride)
      {
        QMessageBox::StandardButton ret;
        ret = QMessageBox::warning(views().first(), "Nodes without incoming edges", QString("Node") + (strange_vertices.size() > 1 ? "s " : " ") + strange_vertices.join(", ") + " will never be reached.\n\nDo you still want to run the pipeline?", QMessageBox::Yes | QMessageBox::No);
        if (ret == QMessageBox::No)
        {
          return false;
        }
      }
      //else
      //{
      // assume the pipeline was tested in the gui, continue
      //}
    }

    // ----- are there nodes without children (besides output nodes)? -----
    strange_vertices.clear();
    for (TOPPASVertex* tv : vertices_)
    {
      if (qobject_cast<TOPPASOutputVertex*>(tv))
      {
        continue;
      }
      if (tv->outEdgesBegin() == tv->outEdgesEnd())
      {
        strange_vertices << QString::number(tv->getTopoNr());
      }
    }
    if (!strange_vertices.empty())
    {
      if (allowUserOverride)
      {
        QMessageBox::StandardButton ret;
        ret = QMessageBox::warning(views().first(), "Nodes without outgoing edges", QString("Node") +
                                   (strange_vertices.size() > 1 ? "s " : " ") + strange_vertices.join(", ") +
                                   (strange_vertices.size() > 1 ? " have " : " has ") +
                                   "no outgoing edges.\n\nDo you still want to run the pipeline?", QMessageBox::Yes | QMessageBox::No);
        if (ret == QMessageBox::No)
        {
          return false;
        }
      }
      //else
      //{
      // assume the pipeline was tested in the gui, continue
      //}
    }

    // check edges
    bool edges_ok = true;
    for (TOPPASEdge* edge : edges_)
    {
      if (edge->getEdgeStatus() != TOPPASEdge::ES_VALID)
      {
        edges_ok = false;
        break;
      }
    }
    if (!edges_ok)
    {
      if (allowUserOverride) 
      {
          QMessageBox::StandardButton ret;
          ret = QMessageBox::warning(views().first(), "Invalid edges detected", "Invalid edges detected. Do you still want to run the pipeline?",
                                      QMessageBox::Yes | QMessageBox::No);
          if (ret == QMessageBox::No)
          {
            return false;
          }
      }
      else 
      { // do not allow silent execution with invalid edges
        return false;
      }
    }

    return true;
  }

  void TOPPASScene::connectVertexSignals(TOPPASVertex* tv)
  {
    connect(tv, &TOPPASVertex::clicked, this, &TOPPASScene::itemClicked);
    connect(tv, &TOPPASVertex::released, this, &TOPPASScene::itemReleased);
    connect(tv, &TOPPASVertex::hoveringEdgePosChanged, this, &TOPPASScene::updateHoveringEdgePos);
    connect(tv, &TOPPASVertex::newHoveringEdge, this, &TOPPASScene::addHoveringEdge);
    connect(tv, &TOPPASVertex::finishHoveringEdge, this, &TOPPASScene::finishHoveringEdge);
    connect(tv, &TOPPASVertex::itemDragged, this, &TOPPASScene::moveSelectedItems);
    connect(tv, &TOPPASVertex::parameterChanged, this, &TOPPASScene::changedParameter);
    connect(tv, &TOPPASVertex::somethingHasChanged, this, [this, tv]() {
      resetDownstream(tv);
      setChanged(true);
    });
  }

  void TOPPASScene::connectToolVertexSignals(TOPPASToolVertex* ttv)
  {
    connect(ttv, &TOPPASToolVertex::toppOutputReady, this, &TOPPASScene::logTOPPOutput);
    connect(ttv, &TOPPASToolVertex::toolStarted,  this, &TOPPASScene::logToolStarted);
    connect(ttv, &TOPPASToolVertex::toolFinished, this, &TOPPASScene::logToolFinished);
    connect(ttv, &TOPPASToolVertex::toolFailed,   this, &TOPPASScene::logToolFailed);
    connect(ttv, &TOPPASToolVertex::toolCrashed,  this, &TOPPASScene::logToolCrashed);

    connect(ttv, &TOPPASToolVertex::toolFailed,          this, &TOPPASScene::pipelineErrorSlot);
    connect(ttv, &TOPPASToolVertex::toolCrashed, [&]() { this->pipelineErrorSlot(); });
  }

  void TOPPASScene::connectMergerVertexSignals(TOPPASMergerVertex* tmv)
  {
    // mergeFailed carries only a message; forward it to pipelineErrorSlot's 'msg' parameter
    // (the old SLOT(pipelineErrorSlot(QString)) matched no slot, as the slot is pipelineErrorSlot(int, QString))
    connect(tmv, &TOPPASMergerVertex::mergeFailed, this, [this](const QString& msg) { pipelineErrorSlot(-1, msg); });
  }

  void TOPPASScene::connectOutputVertexSignals(TOPPASOutputVertex* oflv)
  {
    connect(oflv, &TOPPASOutputVertex::outputFileWritten, this, &TOPPASScene::logOutputFileWritten);
    connect(oflv, &TOPPASOutputVertex::outputFolderNameChanged, this, &TOPPASScene::changedOutputFolder);
  }

  void TOPPASScene::connectEdgeSignals(TOPPASEdge* e)
  {
    TOPPASVertex* source = e->getSourceVertex();
    TOPPASVertex* target = e->getTargetVertex();
    connect(e, &TOPPASEdge::somethingHasChanged, source, &TOPPASVertex::outEdgeHasChanged);
    connect(e, &TOPPASEdge::somethingHasChanged, target, &TOPPASVertex::inEdgeHasChanged);
    connect(e, &TOPPASEdge::somethingHasChanged, this, &TOPPASScene::abortPipeline);
  }

  void TOPPASScene::changedOutputFolder()
  {
    resetDownstream(qobject_cast<TOPPASVertex*>(sender()));
    setChanged(true); // to allow "Store" of pipeline
  }

  void TOPPASScene::changedParameter(const bool /*invalidates_running_pipeline*/)
  {
    // A run uses an immutable snapshot. Any execution-affecting edit cancels it
    // before retained upstream results are considered for a subsequent rerun.
    resetDownstream(qobject_cast<TOPPASVertex*>(sender()));
    setChanged(true);
  }

  void TOPPASScene::loadResources(const TOPPASResources& resources)
  {
    abortPipeline();
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
      if (iflv)
      {
        const QString& key = iflv->getKey();
        const QList<TOPPASResource>& resource_list = resources.get(key);
        QStringList files;
        for (const TOPPASResource& res : resource_list)
        {
          files << res.getLocalFile();
        }
        iflv->setFilenames(files);
        resetDownstream(iflv);
        setChanged(true);
      }
    }
  }

  void TOPPASScene::createResources(TOPPASResources& resources)
  {
    resources.clear();
    QStringList used_keys;
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
      if (iflv)
      {
        QString key = iflv->getKey();
        if (used_keys.contains(key))
        {
          if (gui_)
          {
            QMessageBox::warning(nullptr, "Non-unique input node names", "Some of the input nodes have the same names. Cannot create resource file.");
          }
          else
          {
            std::cerr << "Some of the input nodes have the same names. Cannot create resource file." << std::endl;
          }
          return;
        }
        used_keys << key;
        QList<TOPPASResource> resource_list;
        QStringList files = iflv->getFileNames();
        for (const QString& file : files)
        {
          resource_list << TOPPASResource(file);
        }
        resources.add(key, resource_list);
      }
    }
  }

  TOPPASScene::RefreshStatus TOPPASScene::refreshParameters()
  {
    abortPipeline();
    execution_->previous_result = {};
    execution_->discardUnusedTemporaryRuns();
    bool sane_before = sanityCheck_(false);
    bool change = false;
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(*it);
      if (ttv && ttv->refreshParameters())
      {
        change = true;
      }
    }

    TOPPASScene::RefreshStatus result;
    if (!change)
    {
      result = ST_REFRESH_NOCHANGE;
    }
    else if (!sanityCheck_(false)) 
    {
      if (sane_before)
      {
        result = ST_REFRESH_CHANGEINVALID;
      }
      else
      {
        result = ST_REFRESH_REMAINSINVALID;
      }
    }
    else result = ST_REFRESH_CHANGED;
    
    return result;
  }

  void TOPPASScene::setAllowedThreads(int num_jobs)
  {
    if (num_jobs < 1)
    {
      return;
    }
    allowed_threads_ = num_jobs;
  }

  bool TOPPASScene::isGUIMode() const
  {
    return gui_;
  }


  TOPPASEdge* TOPPASScene::getHoveringEdge()
  {
    return hover_edge_;
  }

} //namespace OpenMS
