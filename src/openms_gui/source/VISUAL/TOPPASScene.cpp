// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Veit $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASResources.h>
#include <OpenMS/VISUAL/TOPPASSplitterVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASVertexNameDialog.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/FORMAT/ParamXMLFile.h>

#include <QApplication>
#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QSet>
#include <QtCore/QTextStream>
#include <QtWidgets/QMessageBox>

#include <map>
#include <OpenMS/VISUAL/TOPPASOutputFolderVertex.h>

namespace OpenMS
{


  void FakeProcess::start(const QString& /*program*/, const QStringList& /*arguments*/, OpenMode /*mode = ReadWrite*/)
  {
    // don't do anything...
    //std::cout << "fake process " << program.toStdString() << " called.\n";
    emit finished(0, QProcess::NormalExit);
  }

  TOPPASScene::TOPPASScene(QObject* parent, const QString& tmp_path, bool gui) :
    QGraphicsScene(parent),
    action_mode_(AM_NEW_EDGE),
    vertices_(),
    edges_(),
    hover_edge_(nullptr),
    potential_target_(nullptr),
    file_name_(),
    tmp_path_(tmp_path),
    gui_(gui),
    out_dir_(File::getUserDirectory().toQString()),
    changed_(false),
    running_(false),
    error_occured_(false),
    user_specified_out_dir_(false),
    clipboard_(nullptr),
    dry_run_(true),
    threads_active_(0),
    allowed_threads_(1),
    resume_source_(nullptr)
  {
    /*	ATTENTION!

            The following line is important! Without it, we get
            hard-to-reproduce segmentation faults and
            "pure virtual method calls" due to a bug in Qt!

            (http://lists.trolltech.com/qt4-preview-feedback/2006-09/thread00124-0.html)
    */
    setItemIndexMethod(QGraphicsScene::NoIndex);
  }

  TOPPASScene::~TOPPASScene()
  {
    // Delete all items in a controlled way:
    foreach(TOPPASVertex* vertex, vertices_)
    {
      vertex->blockSignals(true); // do not propagate changes, remove output files, etc..
      vertex->setSelected(true);
    }
    foreach(TOPPASEdge* edge, edges_)
    {
      edge->blockSignals(true); // do not propagate changes, remove output files, etc..
      edge->setSelected(true);
    }
    removeSelected();
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
    vertices_.push_back(tv);
    addItem(tv);
  }

  void TOPPASScene::addEdge(TOPPASEdge* te)
  {
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
          emit messageReady(String(ss.str()).toQString());
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

    foreach(TOPPASVertex* v, vertices_)
    {
      if (!v->isSelected())
      {
        continue;
      }

      TOPPASVertex* new_v = v->clone().release();

      vertex_map[v] = new_v;
      tmp_scene->addVertex(new_v);
    }

    foreach(TOPPASEdge* e, edges_)
    {
      if (!e->isSelected())
      {
        continue;
      }

      //check if both source and target node were also selected (otherwise don't copy)
      TOPPASVertex* old_source = e->getSourceVertex();
      TOPPASVertex* old_target = e->getTargetVertex();
      if (vertex_map.find(old_source) == vertex_map.end())
      {
        continue;
      }

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

    TOPPASEdge* edge = nullptr;
    foreach(edge, edges_to_be_removed)
    {
      edges_.removeAll(edge);
      removeItem(edge); // remove from scene
      delete edge;
    }
    TOPPASVertex* vertex = nullptr;
    foreach(vertex, vertices_to_be_removed)
    {
      vertices_.removeAll(vertex);
      removeItem(vertex); // remove from scene
      delete vertex;
    }

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
    foreach(TOPPASVertex* vertex, vertices_)
    {
      vertex->setDFSColor(TOPPASVertex::DFS_WHITE);
    }
    foreach(TOPPASVertex* vertex, vertices_)
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
    foreach(TOPPASEdge* edge, edges_)
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
    // reset all nodes
    vertex->reset(true);
    for (TOPPASVertex::ConstEdgeIterator it = vertex->outEdgesBegin(); it != vertex->outEdgesEnd(); ++it)
    {
      TOPPASVertex* target = (*it)->getTargetVertex();
      this->resetDownstream(target);
    }
  }

  void TOPPASScene::runPipeline()
  {
    error_occured_ = false;
    resume_source_ = nullptr; // we are not resuming, so reset the resume node

    // reset all nodes
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      (*it)->reset(true);
    }
    update(sceneRect());

    // check if pipeline OK
    if (!sanityCheck_(gui_))
    {
      if (!gui_)
      {
        emit pipelineExecutionFailed(); // the user cannot interact. End processing.
      }
      return;
    }

    // ask for output directory
    if (!askForOutputDir(true))
    {
      return;
    }

    std::vector<bool> runs;
    runs.push_back(true); // iterate through dry run and normal run
    runs.push_back(false);

    foreach(bool dry_run_state, runs)
    {
      this->dry_run_ = dry_run_state;
      setPipelineRunning();

      std::cout << "current dry-run state: " << dry_run_state << "\n";

      // reset all nodes
      for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
      {
        (*it)->reset(true);
      }
      update(sceneRect());

      // reset logfile
      QFile logfile(out_dir_ + QDir::separator() + "TOPPAS.log");
      if (logfile.exists())
        logfile.remove();

      // reset processes
      topp_processes_queue_.clear();

      // start at input nodes
      for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
      {
        if (error_occured_) break; // someone raised an error

        TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
        if (iflv)
        {
          iflv->run();
        }
      }
    } // foreach
  }

  bool TOPPASScene::store(const String& file)
  {
    Param save_param;

    save_param.setValue("info:version", VersionInfo::getVersion());
    save_param.setValue("info:num_vertices", vertices_.size());
    save_param.setValue("info:num_edges", edges_.size());
    save_param.setValue("info:description", String("<![CDATA[") + String(this->description_text_) + String("]]>"));

    // lambda function to store common parameters of all vertices
    auto save_common_params =
      [&save_param](const TOPPASVertex* tv, const String& id, const String& type)
      {
        save_param.setValue("vertices:" + id + ":toppas_type", type);
        save_param.setValue("vertices:" + id + ":x_pos", tv->x());
        save_param.setValue("vertices:" + id + ":y_pos", tv->y());
        save_param.setValue("vertices:" + id + ":recycle_output", tv->isRecyclingEnabled() ? "true" : "false");
    };
      

    // store all vertices (together with all parameters)
    for (TOPPASVertex * tv : vertices_)
    {
      String id(tv->getTopoNr() - 1);

      // vertex subclasses
      if (auto* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv); iflv)
      {
        // store file names relative to toppas file
        QDir save_dir(File::path(file).toQString());
        const QStringList& files_qt = iflv->getFileNames();
        std::vector<std::string> files;
        foreach(const QString &file_qt, files_qt)
        {
          files.push_back(save_dir.relativeFilePath(file_qt).toStdString());
        }
        save_common_params(iflv, id, "input file list");
        save_param.setValue("vertices:" + id + ":file_names", files);
        continue;
      }
      
      if (auto* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv); oflv)
      {
        save_common_params(oflv, id, "output file list");
        save_param.setValue("vertices:" + id + ":output_folder_name", oflv->getOutputFolderName().toStdString());
        continue;
      }
      
      if (auto* ofv = qobject_cast<TOPPASOutputFolderVertex*>(tv); ofv)
      {
        save_common_params(ofv, id, "output folder");
        save_param.setValue("vertices:" + id + ":output_folder_name", ofv->getOutputFolderName().toStdString());
        continue;
      }

      if (auto* ttv = qobject_cast<TOPPASToolVertex*>(tv); ttv)
      {
        save_common_params(ttv, id, "tool");
        save_param.setValue("vertices:" + id + ":tool_name", ttv->getName());
        save_param.setValue("vertices:" + id + ":tool_type", ttv->getType());
        save_param.insert("vertices:" + id + ":parameters:", ttv->getParam());
        continue;
      }

      if (auto* mv = qobject_cast<TOPPASMergerVertex*>(tv); mv)
      {
        save_common_params(mv, id, "merger");
        save_param.setValue("vertices:" + id + ":round_based", mv->roundBasedMode() ? "true" : "false");
        continue;
      }

      if (auto* sv = qobject_cast<TOPPASSplitterVertex*>(tv); sv)
      {
        save_common_params(sv, id, "splitter");
        continue;
      }
    }

    // store all edges
    int counter = 0;
    for (TOPPASEdge* te : edges_)
    {
      if (!((te->getEdgeStatus() == TOPPASEdge::ES_VALID) || (te->getEdgeStatus() == TOPPASEdge::ES_NOT_READY_YET)))
      { // do not allow to store an invalid pipeline, e.g., after a "param refresh()", since this might lead to inconsistencies when storing the edge mapping parameters (segfaults even).
        // alternatively, we could discard invalid edges during loading, but then the user looses the information where edges were present (currently they become red)
        return false;
      }
      if (!(te->getSourceVertex() && te->getTargetVertex()))
      {
        continue;
      }

      save_param.setValue("edges:" + String(counter) + ":source/target:", String(te->getSourceVertex()->getTopoNr() - 1) + "/" + String(te->getTargetVertex()->getTopoNr() - 1));
      //save_param.setValue("edges:"+String(counter)+":source_out_param:", te->getSourceOutParam()));
      //save_param.setValue("edges:"+String(counter)+":target_in_param:", te->getTargetInParam()));
      String v = "__no_name__";
      if (te->getSourceOutParam() >= 0)
      {
        TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(te->getSourceVertex());
        if (tv_src)
        {
          QVector<TOPPASToolVertex::IOInfo> files = tv_src->getOutputParameters();
          //std::cout << "#p: " << files.size() << " . " << te->getSourceOutParam() << "\n";
          v = files[te->getSourceOutParam()].param_name;
        }
      }
      save_param.setValue("edges:" + String(counter) + ":source_out_param:", v);

      v = "__no_name__";
      if (te->getTargetInParam() >= 0)
      {
        TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(te->getTargetVertex());
        if (tv_src)
        {
          QVector<TOPPASToolVertex::IOInfo> files = tv_src->getInputParameters();
          //std::cout << "#p: " << files.size() << " . " << te->getTargetInParam() << "\n";
          v = files[te->getTargetInParam()].param_name;
        }
      }
      save_param.setValue("edges:" + String(counter) + ":target_in_param:", v);

      ++counter;
    }

    // save file
    ParamXMLFile paramFile;
    paramFile.store(file, save_param);
    setChanged(false);
    file_name_ = file;

    return true; // success
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

  void TOPPASScene::load(const String& file)
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
    String file_version = "1.8.0"; // default (were we did not have the tag)
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
        String extra_quotes = "\""; // note: double quoting required for Windows, as outer quotes are required by cmd.exe (arghh)...
#else
        String extra_quotes = "";
#endif

        String cmd = extra_quotes + "\"" + File::findSiblingTOPPExecutable("INIUpdater") + "\" -in \"" + file + "\" -i " + extra_quotes;
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


    Param vertices_param = load_param.copy("vertices:", true);
    Param edges_param = load_param.copy("edges:", true);

    bool pre_1_9_toppas = true;
    if (load_param.exists("info:version"))
    {
      pre_1_9_toppas = false; // using param names instead of indices for connecting edges
    }
    if (load_param.exists("info:description"))
    {
      String text = String(load_param.getValue("info:description").toString()).toQString();
      text.substitute("<![CDATA[", "");
      text.substitute("]]>", "");
      description_text_ = text.trim().toQString();
    }

    String current_type, current_id;
    TOPPASVertex* current_vertex = nullptr;
    QVector<TOPPASVertex*> vertex_vector;
    vertex_vector.resize((Size)(int)load_param.getValue("info:num_vertices"));

    // load all vertices
    for (Param::ParamIterator it = vertices_param.begin(); it != vertices_param.end(); ++it)
    {
      StringList substrings;
      String(it.getName()).split(':', substrings);
      if (substrings.back() == "toppas_type") // next node (all nodes have a "toppas_type")
      {
        current_vertex = nullptr;
        current_type = (it->value).toString();
        current_id = substrings[0];
        Int index = current_id.toInt();

        if (current_type == "input file list")
        {
          StringList file_names = ListUtils::toStringList<std::string>(vertices_param.getValue(current_id + ":file_names"));
          QStringList file_names_qt;

          for (StringList::const_iterator str_it = file_names.begin(); str_it != file_names.end(); ++str_it)
          {
            QString f = str_it->toQString();
            if (QDir::isRelativePath(f)) // prepend path of toppas file to relative path of the input files
            {
              f = File::path(file).toQString() + "/" + f;
            }
            file_names_qt.push_back(QDir::cleanPath(f));
          }
          TOPPASInputFileListVertex* iflv = new TOPPASInputFileListVertex(file_names_qt);
          current_vertex = iflv;
        }
        else if (current_type == "output file list")
        {
          TOPPASOutputFileListVertex* oflv = new TOPPASOutputFileListVertex();
          // custom output folder
          if (vertices_param.exists(current_id + ":output_folder_name"))
          {
            oflv->setOutputFolderName(String(vertices_param.getValue(current_id + ":output_folder_name").toString()).toQString());
          }
          
          connectOutputVertexSignals(oflv); // todo

          current_vertex = oflv;
        }
        else if (current_type == "output folder")
        {
          auto* ofv = new TOPPASOutputFolderVertex();
          // custom output folder
          if (vertices_param.exists(current_id + ":output_folder_name"))
          {
            ofv->setOutputFolderName(String(vertices_param.getValue(current_id + ":output_folder_name").toString()).toQString());
          }

          connectOutputVertexSignals(ofv);

          current_vertex = ofv;
        }
        else if (current_type == "tool")
        {
          String tool_name = vertices_param.getValue(current_id + ":tool_name").toString();
          String tool_type = vertices_param.getValue(current_id + ":tool_type").toString();
          Param param_param = vertices_param.copy(current_id + ":parameters:", true);
          TOPPASToolVertex* tv = new TOPPASToolVertex(tool_name, tool_type);
          tv->setParam(param_param);

          connectToolVertexSignals(tv);

          current_vertex = tv;
        }
        else if (current_type == "merger")
        {
          String rb = "true";
          if (vertices_param.exists(current_id + ":round_based"))
          {
            rb = vertices_param.getValue(current_id + ":round_based").toString();
          }
          TOPPASMergerVertex* mv = new TOPPASMergerVertex(rb == "true");

          connectMergerVertexSignals(mv);

          current_vertex = mv;
        }
        else if (current_type == "splitter")
        {
          TOPPASSplitterVertex* sv = new TOPPASSplitterVertex();

          current_vertex = sv;
        }
        else
        {
          std::cerr << "Unknown vertex type '" << current_type << "'" << std::endl;
        }

        if (current_vertex)
        {
          float x = vertices_param.getValue(current_id + ":x_pos");
          float y = vertices_param.getValue(current_id + ":y_pos");

          current_vertex->setPos(QPointF(x, y));

          // vertex parameters:
          if (vertices_param.exists(current_id + ":recycle_output")) // only since TOPPAS 1.9, so does not need to exist
          {
            String recycle = vertices_param.getValue(current_id + ":recycle_output").toString();
            current_vertex->setRecycling(recycle == "true" ? true : false);
          }

          addVertex(current_vertex);

          connectVertexSignals(current_vertex);

          // temporarily block signals in order that the first topo sort does not set the changed flag
          current_vertex->blockSignals(true);

          if (index >= vertex_vector.size())
          {
            std::cerr << "Unexpected vertex ID!" << std::endl;
          }
          else
          {
            if (vertex_vector[index] != 0)
            {
              std::cerr << "Vertex occupied!" << std::endl;
            }
            else
            {
              vertex_vector[index] = current_vertex;
            }
          }
        }
        else
        {
          std::cerr << "Current vertex not available." << std::endl;
        }
      }
    }

    // load all edges
    for (Param::ParamIterator it = edges_param.begin(); it != edges_param.end(); ++it)
    {
      const String& edge = (it->value).toString();
      StringList edge_substrings;
      edge.split('/', edge_substrings);
      if (edge_substrings.size() != 2)
      {
        std::cerr << "Invalid edge format" << std::endl;
        break;
      }
      Int index_1 = edge_substrings[0].toInt();
      Int index_2 = edge_substrings[1].toInt();

      if (index_1 >= vertex_vector.size() || index_2 >= vertex_vector.size())
      {
        std::cerr << "Invalid vertex index" << std::endl;
      }
      else
      {
        TOPPASVertex* tv_1 = vertex_vector[index_1];
        TOPPASVertex* tv_2 = vertex_vector[index_2];
        
        // future TOPPAS files may contain new nodes, which may leave `vertex_vector[i]` empty
        if (tv_1 == nullptr || tv_2 == nullptr)
        {
          std::cerr << "Invalid edge" << std::endl;
          continue;
        }

        TOPPASEdge* edge = new TOPPASEdge();
        edge->setSourceVertex(tv_1);
        edge->setTargetVertex(tv_2);
        tv_1->addOutEdge(edge);
        tv_2->addInEdge(edge);

        connectEdgeSignals(edge);

        addEdge(edge);

        String source_out_param = (++it)->value.toString();
        String target_in_param = (++it)->value.toString();
        if (pre_1_9_toppas) // just indices stored - no way we can check
        {
          edge->setSourceOutParam(source_out_param.toInt());
          edge->setTargetInParam(target_in_param.toInt());
        }
        else
        {
          Int src_index = -1;
          Int tgt_index = -1;
          TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(tv_1);
          if (source_out_param != "__no_name__" && tv_src)
          {
            QVector<TOPPASToolVertex::IOInfo> files = tv_src->getOutputParameters();
            // search for the name
            for (int i = 0; i < files.size(); ++i)
            {
              if (files[i].param_name == source_out_param)
              {
                src_index = i;
                break;
              }
            }
            if (src_index == -1)
              logTOPPOutput(String("Could not find output parameter called '" + source_out_param + "'. Check edge!").toQString());
          }

          tv_src = qobject_cast<TOPPASToolVertex*>(tv_2);
          if (target_in_param != "__no_name__" && tv_src)
          {
            QVector<TOPPASToolVertex::IOInfo> files = tv_src->getInputParameters();
            // search for the name
            for (int i = 0; i < files.size(); ++i)
            {
              if (files[i].param_name == target_in_param)
              {
                tgt_index = i;
                break;
              }
            }
            if (tgt_index == -1)
              logTOPPOutput(String("Could not find input parameter called '" + target_in_param + "'. Check edge!").toQString());
          }

          edge->setSourceOutParam(src_index);
          edge->setTargetInParam(tgt_index);
        }
      }
    }
    if (pre_1_9_toppas) // just indices stored - no way we can check
    {
      logTOPPOutput(String("Your TOPPAS file was build with an old version of TOPPAS and is susceptible to errors when used with new versions of OpenMS. Check every edge for correct input/output parameter names and store the workflow using the current version of TOPPAS (e.g using the \"Save as ...\" functionality) to make the workflow more robust to changes in future versions of TOPP tools!").toQString());
    }

/*
    if (!views().empty())
        {
            TOPPASWidget* tw = qobject_cast<TOPPASWidget*>(views().first());
            if (tw)
            {
                QRectF scene_rect = itemsBoundingRect();

                tw->fitInView(scene_rect, Qt::KeepAspectRatio);
                tw->scale(0.75, 0.75);
                setSceneRect(tw->mapToScene(tw->rect()).boundingRect());
            }
        }
*/

    topoSort();
    // unblock signals again
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      (*it)->blockSignals(false);
    }

    updateEdgeColors();
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

  const String& TOPPASScene::getSaveFileName()
  {
    return file_name_;
  }

  void TOPPASScene::setSaveFileName(const String& name)
  {
    file_name_ = name;
  }

  void TOPPASScene::unselectAll()
  {
    const QList<QGraphicsItem*>& all_items = items();
    foreach(QGraphicsItem * item, all_items)
    {
      item->setSelected(false);
    }
    update(sceneRect());
  }

  void TOPPASScene::checkIfWeAreDone()
  {
    if (dry_run_)
      return;

    if (resume_source_)
    {
      switch (resume_source_->getSubtreeStatus())
      {
      case TOPPASVertex::TV_UNFINISHED:
        return; // still processing

        break;

      case TOPPASVertex::TV_ALLFINISHED:
        break; // ok, go to bottom

      case TOPPASVertex::TV_UNFINISHED_INBRANCH:
        setPipelineRunning(false);
        emit pipelineErrorSlot(-1, "Resume cannot continue due to missing subtree.");
        break;
      }
    }
    else
    {
      for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it) // check if all nodes are done
      {
        if (!(*it)->isFinished())
        {
          return;
        }
      }
    }

    setPipelineRunning(false);
    emit entirePipelineFinished();
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
      std::cerr << "Could not write to logfile '" << String(logfile.fileName()) << "'" << std::endl;
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
    String text = String(out);

    if (!gui_)
    {
      std::cout << std::endl << text << std::endl;
    }
    emit messageReady(out); // let TOPPAS know about it

    writeToLogFile_(text.toQString());
  }

  void TOPPASScene::logToolStarted()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " started. Processing ...";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(text.toQString());
    }
  }

  void TOPPASScene::logToolFinished()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " finished!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(text.toQString());
    }
  }

  void TOPPASScene::logToolFailed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " failed!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(text.toQString());
    }
  }

  void TOPPASScene::logToolCrashed()
  {
    TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(QObject::sender());
    if (tv)
    {
      String text = tv->getName();
      String type = tv->getType();
      if (!type.empty())
      {
        text += " (" + type + ")";
      }
      text += " crashed!";

      if (!gui_)
      {
        std::cout << '\n' << text << std::endl;
      }

      writeToLogFile_(text.toQString());
    }
  }

  void TOPPASScene::logOutputFileWritten(const String& file)
  {
    String text = "Output file '" + file + "' written.";

    if (!gui_)
    {
      std::cout << std::endl << text << std::endl;
    }

    writeToLogFile_(text.toQString());
  }

  void TOPPASScene::topoSort(bool resort_all)
  {
    UInt topo_counter {1};
    for (TOPPASVertex* tv : vertices_)
    {
      if (resort_all)
      {
        tv->setTopoSortMarked(false);
      }
      else if (tv->isTopoSortMarked())
      {
        ++topo_counter; // count number of existing/sorted vertices to get correct offset for new vertices
      }
    }
  
    while (true)
    {
      bool some_vertex_not_finished = false;
      for (TOPPASVertex* tv : vertices_)
      {
        if (tv->isTopoSortMarked())
        {
          continue;
        }
        
        bool has_unmarked_predecessors = false;
        for (TOPPASVertex::ConstEdgeIterator e_it = tv->inEdgesBegin(); e_it != tv->inEdgesEnd(); ++e_it)
        {
          TOPPASVertex* v = (*e_it)->getSourceVertex();
          if (!(v->isTopoSortMarked()))
          {
            has_unmarked_predecessors = true;
            break;
          }
        }
        if (has_unmarked_predecessors)
        { // needs to be revisited in the next round (where we hopefully have found the predecessors)
          some_vertex_not_finished = true;
        }
        else
        { // mark this node
          // update name of input node
          TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
          if (iflv)
          {
            //check if key was modified by user. if yes, don't update it
            QString old_topo_nr = QString::number(tv->getTopoNr());
            if (old_topo_nr == iflv->getKey() || iflv->getKey() == "")
            {
              iflv->setKey(QString::number(topo_counter));
            }
          }

          tv->setTopoNr(topo_counter);
          tv->setTopoSortMarked(true);

          ++topo_counter;
        }
      }
      if (!some_vertex_not_finished)
      {
        break; // all sorted
      }
    }

    // sort vertices in list by their TopoNr, so that they keep their numbering when deleting edges
    std::sort(vertices_.begin(), vertices_.end(), [](TOPPASVertex* a, TOPPASVertex* b) { return a->getTopoNr() < b->getTopoNr();});

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
      QString name = file_name_.empty() ? "Untitled" : File::basename(file_name_).toQString();
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
    emit terminateCurrentPipeline();
    resetProcessesQueue();
    setPipelineRunning(false);
  }

  void TOPPASScene::resetProcessesQueue()
  {
    topp_processes_queue_.clear();
  }

  void TOPPASScene::setPipelineRunning(bool b)
  {
    running_ = b;
    if (!running_) // whenever we stop the pipeline and user is not looking, the icon should flash
    {
      resume_source_ = nullptr;
      QApplication::alert(nullptr); // flash Taskbar || Dock
    }
  }

  void TOPPASScene::processFinished()
  {
    --threads_active_;
    // try to run next in line
    runNextProcess();
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

      foreach(TOPPASEdge* edge, edges_)
      {
        if (edge->isSelected())
        {
          found_edge = true;
          break;
        }
      }

      foreach(TOPPASVertex* tv, vertices_)
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
      foreach(const QSet<QString>&action_set, all_actions)
      {
        supported_actions_set.intersect(action_set);
      }

      QList<QString> supported_actions = supported_actions_set.values();
      supported_actions << "Copy" << "Cut" << "Remove";
      foreach(const QString &supported_action, supported_actions)
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

      foreach(QGraphicsItem* gi, selectedItems())
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
            if (askForOutputDir(false))
            {
              setPipelineRunning();
              resume_source_ = ttv;
              resetDownstream(ttv);
              ttv->run();
            }
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

  void TOPPASScene::enqueueProcess(const TOPPProcess& process)
  {
    topp_processes_queue_ << process;
  }

  void TOPPASScene::runNextProcess()
  {
    static bool used = false;
    if (used)
      return;

    used = true;

    while (!topp_processes_queue_.empty() && threads_active_ < allowed_threads_)
    {
      ++threads_active_; // will be decreased, once the tool finishes
      TOPPProcess tp = topp_processes_queue_.first();
      topp_processes_queue_.pop_front();
      FakeProcess* p = qobject_cast<FakeProcess*>(tp.proc);
      if (p)
      {
        p->start(tp.command, tp.args);
      }
      else
      {
        tp.tv->emitToolStarted();
        tp.proc->start(tp.command, tp.args);
      }
    }
    used = false;

    checkIfWeAreDone();
  }

  bool TOPPASScene::sanityCheck_(bool allowUserOverride)
  {
    QStringList strange_vertices;

    // ----- are there any input nodes and are files specified? ----

    /// check if we have any input nodes
    QVector<TOPPASInputFileListVertex*> input_nodes;
    foreach(TOPPASVertex* tv, vertices_)
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
    foreach(TOPPASInputFileListVertex* iflv, input_nodes)
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
    foreach(TOPPASInputFileListVertex* iflv, input_nodes)
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
    foreach(TOPPASVertex* tv, vertices_)
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
    foreach(TOPPASVertex* tv, vertices_)
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
    foreach(TOPPASEdge* edge, edges_)
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
    connect(tv, SIGNAL(clicked()), this, SLOT(itemClicked()));
    connect(tv, SIGNAL(released()), this, SLOT(itemReleased()));
    connect(tv, SIGNAL(hoveringEdgePosChanged(const QPointF &)), this, SLOT(updateHoveringEdgePos(const QPointF &)));
    connect(tv, SIGNAL(newHoveringEdge(const QPointF &)), this, SLOT(addHoveringEdge(const QPointF &)));
    connect(tv, SIGNAL(finishHoveringEdge()), this, SLOT(finishHoveringEdge()));
    connect(tv, SIGNAL(itemDragged(qreal, qreal)), this, SLOT(moveSelectedItems(qreal, qreal)));
    connect(tv, SIGNAL(parameterChanged(const bool)), this, SLOT(changedParameter(const bool)));
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
    connect(ttv, &TOPPASToolVertex::somethingHasChanged, this, &TOPPASScene::abortPipeline);
  }

  void TOPPASScene::connectMergerVertexSignals(TOPPASMergerVertex* tmv)
  {
    connect(tmv, SIGNAL(mergeFailed(QString)), this, SLOT(pipelineErrorSlot(QString)));
    connect(tmv, SIGNAL(somethingHasChanged()), this, SLOT(abortPipeline()));
  }

  void TOPPASScene::connectOutputVertexSignals(TOPPASOutputVertex* oflv)
  {
    connect(oflv, SIGNAL(outputFileWritten(const String &)), this, SLOT(logOutputFileWritten(const String&)));
    connect(oflv, SIGNAL(outputFolderNameChanged()), this, SLOT(changedOutputFolder()));
  }

  void TOPPASScene::connectEdgeSignals(TOPPASEdge* e)
  {
    TOPPASVertex* source = e->getSourceVertex();
    TOPPASVertex* target = e->getTargetVertex();
    connect(e, SIGNAL(somethingHasChanged()), source, SLOT(outEdgeHasChanged()));
    connect(e, SIGNAL(somethingHasChanged()), target, SLOT(inEdgeHasChanged()));
    connect(e, SIGNAL(somethingHasChanged()), this, SLOT(abortPipeline()));
  }

  void TOPPASScene::changedOutputFolder()
  {
    abortPipeline();
    setChanged(true); // to allow "Store" of pipeline
  }

  void TOPPASScene::changedParameter(const bool invalidates_running_pipeline)
  {
    if (invalidates_running_pipeline) // abort only if TTV's new parameters invalidate the results
    {
      abortPipeline();
    }
    setChanged(true); // to allow "Store" of pipeline
    resetDownstream(dynamic_cast<TOPPASVertex*>(sender()));
  }

  void TOPPASScene::loadResources(const TOPPASResources& resources)
  {
    for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
    {
      TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
      if (iflv)
      {
        const QString& key = iflv->getKey();
        const QList<TOPPASResource>& resource_list = resources.get(key);
        QStringList files;
        foreach(const TOPPASResource& res, resource_list)
        {
          files << res.getLocalFile();
        }
        iflv->setFilenames(files);
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
        foreach(const QString& file, files)
        {
          resource_list << TOPPASResource(file);
        }
        resources.add(key, resource_list);
      }
    }
  }

  TOPPASScene::RefreshStatus TOPPASScene::refreshParameters()
  {
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

  bool TOPPASScene::isDryRun() const
  {
    return dry_run_;
  }
  
  void TOPPASScene::quitWithError(int exit_code)
  {
    exit(exit_code);
  }

  TOPPASEdge* TOPPASScene::getHoveringEdge()
  {
    return hover_edge_;
  }

} //namespace OpenMS
