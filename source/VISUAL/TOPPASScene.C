// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Johannes Junker $
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
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASVertexNameDialog.h>
#include <OpenMS/CONCEPT/VersionInfo.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QFile>
#include <QtCore/QFileInfo>
#include <QtCore/QDir>
#include <QtCore/QSet>
#include <QtCore/QTextStream>
#include <QtGui/QMessageBox>

namespace OpenMS
{

  
  void FakeProcess::start ( const QString & /*program*/, const QStringList & /*arguments*/, OpenMode /*mode = ReadWrite*/ )
  {
    // don't do anything...
    //std::cout << "fake process " << program.toStdString() << " called.\n";
    emit finished ( 0, QProcess::NormalExit);
  }

	TOPPASScene::TOPPASScene(QObject* parent, const QString& tmp_path, bool gui)
		:	QGraphicsScene(parent),
			action_mode_(AM_NEW_EDGE),
			vertices_(),
			edges_(),
			hover_edge_(0),
			potential_target_(0),
			file_name_(),
			tmp_path_(tmp_path),
			gui_(gui),
			out_dir_(File::getUserDirectory().toQString()),
			changed_(false),
			running_(false),
			user_specified_out_dir_(false),
			clipboard_(0),
      dry_run_(true),
      threads_active_(0)
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
		foreach (TOPPASVertex* vertex, vertices_)
		{
			vertex->blockSignals(true); // do not propagate changes, remove output files, etc..
			vertex->setSelected(true);
		}
		foreach (TOPPASEdge* edge, edges_)
		{
			edge->blockSignals(true); // do not propagate changes, remove output files, etc..
			edge->setSelected(true);
		}
		removeSelected();

    // delete temporary files (TODO: make this a user dialog and ask - for later resume)
    // safety measure: only delete if subdirectory of Temp path; we do not want to delete / or c:
    if (String(tmp_path_).substitute("\\","/").hasPrefix(File::getTempDirectory().substitute("\\","/") + "/")) 
    {
      File::removeDirRecursively(tmp_path_);
    }
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
			potential_target_ = 0;
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
		
		if (target && 
				target != hover_edge_->getSourceVertex() &&
				isEdgeAllowed_(hover_edge_->getSourceVertex(), target))
		{
			hover_edge_->setTargetVertex(target);
			TOPPASVertex* source = hover_edge_->getSourceVertex();
			source->addOutEdge(hover_edge_);
			target->addInEdge(hover_edge_);
			hover_edge_->setColor(QColor(255,165,0));
			connect (source, SIGNAL(somethingHasChanged()), hover_edge_, SLOT(sourceHasChanged()));
			connect (hover_edge_, SIGNAL(somethingHasChanged()), target, SLOT(inEdgeHasChanged()));
			
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
		
		if (remove_edge)
		{
			edges_.removeAll(hover_edge_);
			removeItem(hover_edge_);
			delete hover_edge_;
			hover_edge_ = 0;
		}
		
		topoSort();
		updateEdgeColors();
	}
	
	TOPPASVertex* TOPPASScene::getVertexAt_(const QPointF& pos)
	{
		QList<QGraphicsItem*> target_list = items(pos);
		
		// return first item that is a vertex
		TOPPASVertex* target = 0;
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
		TOPPASScene* tmp_scene = new TOPPASScene(0, this->getTempDir(), false);
		Map<TOPPASVertex*,TOPPASVertex*> vertex_map;

		foreach (TOPPASVertex* v, vertices_)
		{
			if (!v->isSelected())
			{
				continue;
			}

			TOPPASVertex* new_v = 0;
			
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
			}

			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(v);
			if (tv)
			{
				TOPPASToolVertex* new_tv = new TOPPASToolVertex(*tv);
				new_v = new_tv;
			}
			
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(v);
			if (mv)
			{
				TOPPASMergerVertex* new_mv = new TOPPASMergerVertex(*mv);
				new_v = new_mv;
			}
			
			if (!new_v)
			{
				std::cerr << "Unknown vertex type! Aborting." << std::endl;
				return;
			}
			
			vertex_map[v] = new_v;
			tmp_scene->addVertex(new_v);
		}

		foreach (TOPPASEdge* e, edges_)
		{
			if (!e->isSelected())
			{
				continue;
			}
			
			//check if both source and target node were also selected (otherwise don't copy)
			TOPPASVertex* old_source = e->getSourceVertex();
			TOPPASVertex* old_target = e->getTargetVertex();
			if (!(vertex_map.has(old_source) && vertex_map.has(old_target)))
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
		
		if (clipboard_ != 0)
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
				for (TOPPASVertex::EdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
				{
					(*e_it)->setSelected(true);
				}
				for (TOPPASVertex::EdgeIterator e_it = (*it)->outEdgesBegin(); e_it != (*it)->outEdgesEnd(); ++e_it)
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

		TOPPASEdge* edge;
		foreach (edge, edges_to_be_removed)
		{
			edges_.removeAll(edge);
			removeItem(edge); // remove from scene 
			delete edge;
		}
		TOPPASVertex* vertex;
		foreach (vertex, vertices_to_be_removed)
		{
			vertices_.removeAll(vertex);
			removeItem(vertex); // remove from scene
			delete vertex;
		}
		
		topoSort();
		updateEdgeColors();
	}
	
	bool TOPPASScene::isEdgeAllowed_(TOPPASVertex* u, TOPPASVertex* v)
	{
		if (u == 0 ||
				v == 0 ||
				u == v ||
				// edges leading to input files make no sense:
				qobject_cast<TOPPASInputFileListVertex*>(v) ||
				// neither do edges coming from output files:
				qobject_cast<TOPPASOutputFileListVertex*>(u) ||
				// nor edges from input to output without a tool in between:
				(qobject_cast<TOPPASInputFileListVertex*>(u)
					&& qobject_cast<TOPPASOutputFileListVertex*>(v)) ||
				// nor multiple incoming edges for a single output file/list node
				(qobject_cast<TOPPASOutputFileListVertex*>(v)
					&& v->inEdgesBegin() != v->inEdgesEnd()) ||
				// nor mergers connected directly to an output node
				(qobject_cast<TOPPASMergerVertex*>(u)
					&& qobject_cast<TOPPASOutputFileListVertex*>(v)))
		{
			return false;
		}
		
		// does this edge already exist?
		for (TOPPASVertex::EdgeIterator it = u->outEdgesBegin(); it != u->outEdgesEnd(); ++it)
		{
			if ((*it)->getTargetVertex() == v)
			{
				return false;
			}
		}
		
		//insert edge between u and v for testing, is removed afterwards
		TOPPASEdge* test_edge = new TOPPASEdge(u,QPointF());
		test_edge->setTargetVertex(v);
		u->addOutEdge(test_edge);
		v->addInEdge(test_edge);
		addEdge(test_edge);
		
		bool graph_has_cycles = false;
		//find back edges via DFS
		foreach (TOPPASVertex* vertex, vertices_)
		{
			vertex->setDFSColor(TOPPASVertex::DFS_WHITE);
			vertex->setDFSParent(0);
		}
		foreach (TOPPASVertex* vertex, vertices_)
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
		
		// remove priorly inserted edge
		edges_.removeAll(test_edge);
		removeItem(test_edge);
		delete test_edge;
		
		return !graph_has_cycles;
	}
	
	void TOPPASScene::updateEdgeColors()
	{
		foreach (TOPPASEdge* edge, edges_)
		{
			edge->updateColor();
		}
		update(sceneRect());
	}
	
	bool TOPPASScene::dfsVisit_(TOPPASVertex* vertex)
	{
		vertex->setDFSColor(TOPPASVertex::DFS_GRAY);
		for (TOPPASVertex::EdgeIterator it = vertex->outEdgesBegin(); it != vertex->outEdgesEnd(); ++it)
		{
			TOPPASVertex* target = (*it)->getTargetVertex();
			if (target->getDFSColor() == TOPPASVertex::DFS_WHITE)
			{
				target->setDFSParent(vertex);
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
    //reset all nodes
    vertex->reset(true);
		for (TOPPASVertex::EdgeIterator it = vertex->outEdgesBegin(); it != vertex->outEdgesEnd(); ++it)
		{
      TOPPASVertex* target = (*it)->getTargetVertex();
      this->resetDownstream(target);
    }
  }

	void TOPPASScene::runPipeline()
	{
		
    error_occured_ = false;
    
    //reset all nodes
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			(*it)->reset(true);
		}
		update(sceneRect());
		
		//check if pipeline OK
		if (!sanityCheck())
		{
      if (!gui_) emit pipelineExecutionFailed(); // the user cannot interact. End processing.
			return;
		}
		
		//ask for output directory
		if (!askForOutputDir(true))
		{
			return;
		}
		
    std::vector<bool> runs;
    runs.push_back(true);  // iterate through dry run and normal run
    runs.push_back(false);

    foreach (bool dry_run_state, runs)
    {
      this->dry_run_ = dry_run_state;
      setPipelineRunning();

      std::cout << "current dry-run state: " << dry_run_state << "\n";

      //reset all nodes
		  for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		  {
			  (*it)->reset(true);
		  }
		  update(sceneRect());

		  //reset logfile
		  QFile logfile(out_dir_ + QDir::separator() + "TOPPAS.log");
		  if (logfile.exists())	logfile.remove();
		
		  //reset processes
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

	void TOPPASScene::store(const String& file)
	{
		Param save_param;
		
    save_param.setValue("info:version", DataValue(VersionInfo::getVersion()));
		save_param.setValue("info:num_vertices", DataValue(vertices_.size()));
		save_param.setValue("info:num_edges", DataValue(edges_.size()));
    save_param.setValue("info:description", DataValue(String("<![CDATA[") + String(this->description_text_) + String("]]>")));
		
		// store all vertices (together with all parameters)
		foreach (TOPPASVertex* tv, vertices_)
		{
			String id(tv->getTopoNr()-1);
			
      // common for all vertices
      save_param.setValue("vertices:"+id+":recycle_output", DataValue(tv->isRecyclingEnabled() ? "true" : "false"));

      // vertex subclasses
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
			if (iflv)
			{
				// store file names relative to toppas file
				QDir save_dir(File::path(file).toQString());
        const QStringList& files_qt = iflv->getFileNames();
				StringList files;
				foreach (const QString& file_qt, files_qt)
				{
					files.push_back(save_dir.relativeFilePath(file_qt));
				}
				save_param.setValue("vertices:"+id+":toppas_type", DataValue("input file list"));
				save_param.setValue("vertices:"+id+":file_names", DataValue(files));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
			if (oflv)
			{
				save_param.setValue("vertices:"+id+":toppas_type", DataValue("output file list"));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
			if (ttv)
			{
				save_param.setValue("vertices:"+id+":toppas_type", DataValue("tool"));
				save_param.setValue("vertices:"+id+":tool_name", DataValue(ttv->getName()));
				save_param.setValue("vertices:"+id+":tool_type", DataValue(ttv->getType()));
				save_param.insert("vertices:"+id+":parameters:", ttv->getParam());
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}

			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(tv);
			if (mv)
			{
				save_param.setValue("vertices:"+id+":toppas_type", DataValue("merger"));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				save_param.setValue("vertices:"+id+":round_based", DataValue(mv->roundBasedMode() ? "true" : "false"));
				continue;
			}
		}
		
		//store all edges
		int counter = 0;
		foreach (TOPPASEdge* te, edges_)
		{
			if (!(te->getSourceVertex() && te->getTargetVertex()))
			{
				continue;
			}
			
			save_param.setValue("edges:"+String(counter)+":source/target:", DataValue(String(te->getSourceVertex()->getTopoNr()-1) + "/" + String(te->getTargetVertex()->getTopoNr()-1)));
			//save_param.setValue("edges:"+String(counter)+":source_out_param:", DataValue(te->getSourceOutParam()));
			//save_param.setValue("edges:"+String(counter)+":target_in_param:", DataValue(te->getTargetInParam()));
      QVector<TOPPASToolVertex::IOInfo> files;
      String v = "__no_name__";
      if (te->getSourceOutParam()>=0)
      {
        TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(te->getSourceVertex());
        if (tv_src)
        {
          tv_src->getOutputParameters(files);
          //std::cout << "#p: " << files.size() << " . " << te->getSourceOutParam() << "\n";
          v = files[te->getSourceOutParam()].param_name;
        }
      }
      save_param.setValue("edges:"+String(counter)+":source_out_param:", DataValue(v));
        
      v = "__no_name__";
      if (te->getTargetInParam()>=0)
      {
        TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(te->getTargetVertex());
        if (tv_src)
        {
          tv_src->getInputParameters(files);
          //std::cout << "#p: " << files.size() << " . " << te->getTargetInParam() << "\n";
          v = files[te->getTargetInParam()].param_name;
        }
      }
      save_param.setValue("edges:"+String(counter)+":target_in_param:", DataValue(v));

			++counter;
		}
		
		//save file
		save_param.store(file);
		setChanged(false);
		file_name_ = file;
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
    Param load_param;
    load_param.load(file);


    // check for TOPPAS file version. Deny loading if too old or too new
    // get version of TOPPAS file
    String file_version = "1.8.0"; // default (were we did not have the tag)
    if (load_param.exists("info:version")) file_version = load_param.getValue("info:version");
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
        if (QMessageBox::warning(0, tr("Old TOPPAS file -- convert and override?"), tr("The TOPPAS file you downloaded was created with an old incompatible version of TOPPAS.\n"
                                      "Shall we try to convert the file?! The original file will be overridden, but a backup file will be saved in the same directory.\n")
                                      , QMessageBox::Yes, QMessageBox::No) == QMessageBox::No)
        {
          return;
        }
        // only update in GUI mode, as in non-GUI mode, we'd create infinite recursive calls when instantiating TOPPASScene in INIUpdater
#ifdef OPENMS_WINDOWSPLATFORM
        String extra_quotes = "\""; // note: double quoting required for Windows, as outer quotes are required by cmd.exe (arghh)...
#else
        String extra_quotes = "";
#endif

        String cmd = extra_quotes + "\"" + File::getExecutablePath() + "INIUpdater\" -in \"" + file + "\" -i " + extra_quotes;
        std::cerr << cmd << "\n\n";
        if (std::system(cmd.c_str()))
        {
          QMessageBox::warning(0, tr("INIUpdater failed"), tr("Updating using the INIUpdater tool failed. Please submit a bug report!\n"), QMessageBox::Ok);
          return;
        }
        // reload updated file
        load_param.load(file);
      }
    }
    else if (v_file > v_this_high)
    {
      if (this->gui_ && QMessageBox::warning(0, tr("TOPPAS file too new"), tr("The TOPPAS file you downloaded was created with a more recent version of TOPPAS. Shall we will try to open it?\n"
        "If this fails, update to the new TOPPAS version.\n"), QMessageBox::Yes, QMessageBox::No) == QMessageBox::No) return;
    }


    Param vertices_param = load_param.copy("vertices:",true);
    Param edges_param = load_param.copy("edges:",true);
    
    bool pre_1_9_toppas = true;
    if (load_param.exists("info:version")) pre_1_9_toppas = false; // using param names instead of indices for connecting edges

    if (load_param.exists("info:description"))
    {
      String text = String(load_param.getValue("info:description")).toQString();
      text.substitute("<![CDATA[", "");
      text.substitute("]]>", "");
      description_text_ = text.trim().toQString();
    }

    String current_type, current_id;
    TOPPASVertex* current_vertex = 0;
    QVector<TOPPASVertex*> vertex_vector;
    vertex_vector.resize((Size)(int)load_param.getValue("info:num_vertices"));
    
    //load all vertices
    for (Param::ParamIterator it = vertices_param.begin(); it != vertices_param.end(); ++it)
    {
      StringList substrings;
      it.getName().split(':', substrings);
      if (substrings.back() == "toppas_type") // next node (all nodes have a "toppas_type")
      {
        current_vertex = 0;
      	current_type = (it->value).toString();
      	current_id = substrings[0];
     		Int index = current_id.toInt();
      
				if (current_type == "input file list")
				{
					StringList file_names = vertices_param.getValue(current_id + ":file_names");
					QStringList file_names_qt;
					for (StringList::const_iterator str_it = file_names.begin(); str_it != file_names.end(); ++str_it)
					{
						file_names_qt.push_back(QDir::cleanPath(str_it->toQString()));
					}
					TOPPASInputFileListVertex* iflv = new TOPPASInputFileListVertex(file_names_qt);
					current_vertex = iflv;
				}
				else if (current_type == "output file list")
				{
					TOPPASOutputFileListVertex* oflv = new TOPPASOutputFileListVertex();
					
					connectOutputVertexSignals(oflv);
					
					current_vertex = oflv;
				}
				else if (current_type == "tool")
				{
					String tool_name = vertices_param.getValue(current_id + ":tool_name");
					String tool_type = vertices_param.getValue(current_id + ":tool_type");
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
            rb = vertices_param.getValue(current_id + ":round_based");
					}
          TOPPASMergerVertex* mv = new TOPPASMergerVertex(rb == "true");

          connectMergerVertexSignals(mv);

					current_vertex = mv;
				}
				else
				{
					std::cerr << "Unknown vertex type '" << current_type << "'" << std::endl;
				}
				
				if (current_vertex)
				{
					float x = vertices_param.getValue(current_id + ":x_pos");
					float y = vertices_param.getValue(current_id + ":y_pos");
					
					current_vertex->setPos(QPointF(x,y));
					
          // vertex parameters:
          if (vertices_param.exists(current_id + ":recycle_output")) // only since TOPPAS 1.9, so does not need to exist
          {
            String recycle = vertices_param.getValue(current_id + ":recycle_output");
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
    
    //load all edges
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
      	
      	TOPPASEdge* edge = new TOPPASEdge();
      	edge->setSourceVertex(tv_1);
      	edge->setTargetVertex(tv_2);
      	tv_1->addOutEdge(edge);
      	tv_2->addInEdge(edge);
      	
      	connectEdgeSignals(edge);
      	
      	addEdge(edge);
				
        String source_out_param = (++it)->value;
        String target_in_param = (++it)->value;
        if (pre_1_9_toppas)
        { // just indices stored - no way we can check
          edge->setSourceOutParam(source_out_param.toInt());
          edge->setTargetInParam(target_in_param.toInt());
        }
        else
        {
          QVector<TOPPASToolVertex::IOInfo> files;
          Int src_index = -1;
          Int tgt_index = -1;
          TOPPASToolVertex* tv_src = qobject_cast<TOPPASToolVertex*>(tv_1);
          if (source_out_param != "__no_name__" && tv_src)
          {
            tv_src->getOutputParameters(files);
            // search for the name
            for (int i=0;i<files.size();++i)
            {
              if (files[i].param_name == source_out_param)
              {
                src_index = i;
                break;
              }
            }
            if (src_index==-1) logTOPPOutput(String("Could not find output parameter called '" + source_out_param + "'. Check edge!").toQString());
          }

          tv_src = qobject_cast<TOPPASToolVertex*>(tv_2);
          if (target_in_param != "__no_name__" && tv_src)
          {
            tv_src->getInputParameters(files);
            // search for the name
            for (int i=0;i<files.size();++i)
            {
              if (files[i].param_name == target_in_param)
              {
                tgt_index = i;
                break;
              }
            }
            if (tgt_index==-1) logTOPPOutput(String("Could not find input parameter called '" + target_in_param + "'. Check edge!").toQString());
          }

          edge->setSourceOutParam(src_index);
          edge->setTargetInParam(tgt_index);
        }
      }
    }
    if (pre_1_9_toppas)
    { // just indices stored - no way we can check
      logTOPPOutput(String("Your TOPPAS file was build with an old version of TOPPAS and is susceptible to errors when used with new versions of OpenMS. "
                           "Check every edge for correct input/output parameter names and store the workflow using the current version of TOPPAS "
                           "(e.g using the \"Save as ...\" functionality to make the workflow more robust to changes in future versions of TOPP tools!").toQString());
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
		file_name_ = file;
		
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
		QRectF new_bounding_rect = tmp_scene->itemsBoundingRect();
		QRectF our_bounding_rect = itemsBoundingRect();
		qreal x_offset, y_offset;
		if (pos == QPointF())
		{
			y_offset = our_bounding_rect.bottom() - new_bounding_rect.top() + 40.0;
			x_offset = our_bounding_rect.left() - new_bounding_rect.left();
		}
		else
		{
			x_offset = pos.x() - new_bounding_rect.left();
			y_offset = pos.y() - new_bounding_rect.top();
		}
		Map<TOPPASVertex*,TOPPASVertex*> vertex_map;

		for (VertexIterator it = tmp_scene->verticesBegin(); it != tmp_scene->verticesEnd(); ++it)
		{
			TOPPASVertex* v = *it;
			TOPPASVertex* new_v = 0;
			
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
			TOPPASEdge* new_e = new TOPPASEdge();
			TOPPASVertex* old_source = (*it)->getSourceVertex();
			TOPPASVertex* old_target = (*it)->getTargetVertex();
			TOPPASVertex* new_source = vertex_map[old_source];
			TOPPASVertex* new_target = vertex_map[old_target];
			new_e->setSourceVertex(new_source);
			new_e->setTargetVertex(new_target);
			new_e->setSourceOutParam((*it)->getSourceOutParam());
			new_e->setTargetInParam((*it)->getTargetInParam());
			new_source->addOutEdge(new_e);
			new_target->addInEdge(new_e);
			
			connectEdgeSignals(new_e);
			
			addEdge(new_e);
		}

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
		foreach (QGraphicsItem* item, all_items)
		{
			item->setSelected(false);
		}
		update(sceneRect());
	}
	
	void TOPPASScene::checkIfWeAreDone()
	{
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{ // check if all output nodes are done
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
			if (oflv && !oflv->isFinished())
			{
				return;
			}
		}
		
		running_ = false;
		emit entirePipelineFinished();
	}
	
	void TOPPASScene::pipelineErrorSlot(const QString /*msg*/)
	{
		running_ = false;
    error_occured_ = true;
    abortPipeline();
		emit pipelineExecutionFailed();
	}
	
	void TOPPASScene::writeToLogFile_(const QString& text)
	{
		QFile logfile(out_dir_+QDir::separator()+"TOPPAS.log");
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
			std::cout	<< std::endl << text << std::endl;
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
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " started. Processing ...";
			
			if (!gui_)
			{
				std::cout	<< std::endl << text << std::endl;
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
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " finished!";
			
			if (!gui_)
			{
				std::cout	<< std::endl << text << std::endl;
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
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " failed!";
			
			if (!gui_)
			{
				std::cout	<< std::endl << text << std::endl;
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
			if (type != "")
			{
				text += " ("+type+")";
			}
			text += " crashed!";
			
			if (!gui_)
			{
				std::cout	<< std::endl << text << std::endl;
			}
			
			writeToLogFile_(text.toQString());
		}
	}
	
	void TOPPASScene::logOutputFileWritten(const String& file)
	{
		String text = "Output file '" + file + "' written.";
		
		if (!gui_)
		{
			std::cout	<< std::endl << text << std::endl;
		}
			
		writeToLogFile_(text.toQString());
	}
	
	void TOPPASScene::topoSort()
	{
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			(*it)->setTopoSortMarked(false);
		}
		
		bool topo_sort_finished = false;
		UInt topo_counter = 1;
		while (!topo_sort_finished)
		{
			bool some_vertex_not_finished = false;
			for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
			{
				if ((*it)->isTopoSortMarked())
				{
					continue;
				}
				some_vertex_not_finished = true;
				bool has_predecessors = false;
				for (TOPPASVertex::EdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
				{
					TOPPASVertex* v = (*e_it)->getSourceVertex();
					if (!(v->isTopoSortMarked()))
					{
						has_predecessors = true;
						break;
					}
				}
				if (!has_predecessors)
				{
					//update name of input node
					TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
					if (iflv)
					{
						//check if key was modified by user. if yes, don't update it
						QString old_topo_nr = QString::number((*it)->getTopoNr());
						if (old_topo_nr == iflv->getKey() || iflv->getKey() == "")
						{
							iflv->setKey(QString::number(topo_counter));
						}
					}
					
					(*it)->setTopoNr(topo_counter);
					(*it)->setTopoSortMarked(true);
					
					++topo_counter;
				}
			}
			if (!some_vertex_not_finished)
			{
				topo_sort_finished = true;
			}
		}
		
		update(sceneRect());
	}

	const QString& TOPPASScene::getOutDir()
	{
		return out_dir_;
	}
	
  const QString& TOPPASScene::getTempDir()
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
			for (TOPPASVertex::EdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
			{
				(*e_it)->prepareResize();
			}
			for (TOPPASVertex::EdgeIterator e_it = (*it)->outEdgesBegin(); e_it != (*it)->outEdgesEnd(); ++e_it)
			{
				(*e_it)->prepareResize();
			}
			
			(*it)->moveBy(dx,dy);
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

			if (x_int - prev_grid_x > (grid_step/2))
			{
				new_x += grid_step;
			}
			if (y_int - prev_grid_y > (grid_step/2))
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
			QString name = file_name_ == "" ? "Untitled" : File::basename(file_name_).toQString();
			QMessageBox::StandardButton ret;
			ret = QMessageBox::warning(views().first(), "Save changes?",
						"'"+name+"' has been modified.\n\nDo you want to save your changes?",
						QMessageBox::Save | QMessageBox::Discard | QMessageBox::Cancel);
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
	
	bool TOPPASScene::wasChanged()
	{
		return changed_;
	}
	
	bool TOPPASScene::isPipelineRunning()
	{
		return running_;
	}
	
	void TOPPASScene::abortPipeline()
	{
		emit terminateCurrentPipeline();
		resetProcessesQueue();
		running_ = false;
	}
	
	void TOPPASScene::resetProcessesQueue()
	{
		topp_processes_queue_.clear();
	}
	
	void TOPPASScene::setPipelineRunning(bool b)
	{
		running_ = b;
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
				TOPPASOutputFilesDialog tofd(out_dir_);
				if (tofd.exec())
				{
          setOutDir( tofd.getDirectory() );
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
		QGraphicsItem* clicked_item = itemAt(scene_pos);
		QMenu menu;

		if (clicked_item == 0)
		{
			QAction* new_action = menu.addAction("Paste");
			emit requestClipboardContent();
			if (clipboard_ == 0)
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
			bool found_output = false;
			bool found_merger = false;
			bool found_edge = false;
			bool disable_resume = false;
			//bool disable_toppview = true;

			foreach (TOPPASEdge* edge, edges_)
			{
				if (edge->isSelected())
				{
					found_edge = true;
					break;
				}
			}

			foreach (TOPPASVertex* tv, vertices_)
			{
				if (!tv->isSelected())
				{
					continue;
				}

				if (qobject_cast<TOPPASToolVertex*>(tv))
				{
					found_tool = true;
					// all predecessor nodes finished successfully? if not, disable resuming
					for (EdgeIterator it = tv->inEdgesBegin(); it != tv->inEdgesEnd(); ++it)
					{
						TOPPASToolVertex* pred_ttv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
						if (pred_ttv && (pred_ttv->getProgressColor() != Qt::green || !pred_ttv->isFinished()))
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
				if (qobject_cast<TOPPASOutputFileListVertex*>(tv))
				{
					found_output = true;
					continue;
				}
				if (qobject_cast<TOPPASMergerVertex*>(tv))
				{
					found_merger = true;
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
        action.insert("Toggle breakpoint");
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
				action.insert("Open files in TOPPView");
				action.insert("Open containing folder");
			}

			if (found_edge)
			{
				action.insert("Edit I/O mapping");
			}

      if (found_input || found_tool || found_merger)
      {
        action.insert("Toggle recycling mode");
      }

 			QList< QSet<QString> > all_actions;
			all_actions.push_back(action);

			QSet<QString> supported_actions_set = all_actions.first();
			foreach (const QSet<QString>& action_set, all_actions)
			{
				supported_actions_set.intersect(action_set);
			}
			QList<QString> supported_actions = supported_actions_set.toList();
			supported_actions << "Copy" << "Cut" << "Remove";
			foreach (const QString& supported_action, supported_actions)
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
			
			foreach (QGraphicsItem* gi, selectedItems())
			{

        if (text == "Toggle recycling mode")
        {
          TOPPASVertex* tv = dynamic_cast<TOPPASVertex*>(gi);
          if (tv)
          {
            tv->invertRecylingMode();
            tv->update(tv->boundingRect());
          }
          continue;
        }

				TOPPASEdge* edge = dynamic_cast<TOPPASEdge*>(gi);
				if (edge)
				{
					if (text == "Edit I/O mapping")
					{
						edge->showIOMappingDialog();
					}
					
					continue;
				}
				
				TOPPASToolVertex* ttv = dynamic_cast<TOPPASToolVertex*>(gi);
				if (ttv)
				{
					if (text == "Edit parameters")
					{
						ttv->editParam();
					}
					else if (text == "Resume")
					{
						if (askForOutputDir(false))
						{
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
				
				TOPPASInputFileListVertex* ifv = dynamic_cast<TOPPASInputFileListVertex*>(gi);
				if (ifv)
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
				
				TOPPASOutputFileListVertex* ofv = dynamic_cast<TOPPASOutputFileListVertex*>(gi);
				if (ofv)
				{
					if (text == "Open files in TOPPView")
					{
            QStringList out_files = ofv->getFileNames();
            emit openInTOPPView(out_files);
					}
					else if (text == "Open containing folder")
					{
						ofv->openContainingFolder();
					}

					continue;
				}
			}
		}
		
		event->accept();
	}	

	void TOPPASScene::enqueueProcess(QProcess* p, const QString& command, const QStringList& args)
	{
		topp_processes_queue_ << TOPPProcess(p, command, args);
	}
	
	void TOPPASScene::runNextProcess()
	{
    static bool used = false;

    if (used) return;
    
    used = true;

    int allowed_threads = 1; // change as desired

		while (!topp_processes_queue_.empty() && threads_active_ < allowed_threads)
		{
      ++threads_active_; // will be decreased, once the tool finishes
			TOPPProcess tp = topp_processes_queue_.first();
      topp_processes_queue_.pop_front();
      FakeProcess* p = qobject_cast<FakeProcess*>(tp.proc);
		  if (p)
        p->start(tp.command, tp.args);
      else
        tp.proc->start(tp.command, tp.args);
		}

    used = false;

	}
	
	bool TOPPASScene::sanityCheck()
	{
		QStringList strange_vertices;
		
		// ----- are there any input nodes and are files specified? ----

    /// check if we have any input nodes
		QVector<TOPPASInputFileListVertex*> input_nodes;
		foreach (TOPPASVertex* tv, vertices_)
		{
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
			if (iflv)
			{
				input_nodes.push_back(iflv);
			}
		}
		if (input_nodes.empty())
		{
			if (gui_)
			{
				QMessageBox::warning(0,"No input files","The pipeline does not contain any input file nodes!");
			}
			else
			{
				std::cerr << "The pipeline does not contain any input file nodes!" << std::endl;
			}
			return false;
		}

    /// warn about empty input nodes
		foreach (TOPPASInputFileListVertex* iflv, input_nodes)
		{
      if (iflv->getFileNames().empty())
			{
				strange_vertices.push_back(QString::number(iflv->getTopoNr()));
			}
		}
		if (!strange_vertices.empty())
		{
			if (gui_)
			{
				QMessageBox::warning(views().first(), "Empty input file nodes",
																QString("Node")
																+(strange_vertices.size()>1 ? "s " : " ")
																+strange_vertices.join(", ")
																+(strange_vertices.size()>1 ? " have " : " has ")
																+" an empty input file list!");
			}
			else
			{
				std::cerr << "Pipeline contains input file nodes without specified files!" << std::endl;
			}
			return false;
		}

    /// check if input files exist
		strange_vertices.clear();
		foreach (TOPPASInputFileListVertex* iflv, input_nodes)
		{
      if (!iflv->fileNamesValid())
      {
				strange_vertices.push_back(QString::number(iflv->getTopoNr()));
      }
		}
		if (!strange_vertices.empty())
		{
			if (gui_)
			{
				QMessageBox::warning(views().first(), "Input file names wrong",
																QString("Node")
																+(strange_vertices.size()>1 ? "s " : " ")
																+strange_vertices.join(", ")
																+(strange_vertices.size()>1 ? " have " : " has ")
																+" invalid (non-existing) input files!");
			}
			else
			{
				std::cerr << "Pipeline contains input file nodes with invalid (non-existing) input files!" << std::endl;
			}
			return false;
		}
		
		// ----- are there nodes without parents (besides input nodes)? -----
		strange_vertices.clear();
		foreach (TOPPASVertex* tv, vertices_)
		{
			if (qobject_cast<TOPPASInputFileListVertex*>(tv))
			{ // input nodes don't need a parent
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
			if (gui_)
			{
				QMessageBox::StandardButton ret;
				ret = QMessageBox::warning(views().first(), "Nodes without incoming edges",
																	 QString("Node")
																	 +(strange_vertices.size()>1 ? "s " : " ")
																	 +strange_vertices.join(", ")
																	 +" will never be reached.\n\nDo you still want to run the pipeline?",
																	 QMessageBox::Yes | QMessageBox::No);
				if (ret == QMessageBox::No)
				{
					return false;
				}
			}
			//else
			//{
			//	assume the pipeline was tested in the gui, continue
			//}
		}
		
		// ----- are there nodes without children (besides output nodes)? -----
		strange_vertices.clear();
		foreach (TOPPASVertex* tv, vertices_)
		{
			if (qobject_cast<TOPPASOutputFileListVertex*>(tv))
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
			if (gui_)
			{
				QMessageBox::StandardButton ret;
				ret = QMessageBox::warning(views().first(), "Nodes without outgoing edges",
																	 QString("Node")
																	 +(strange_vertices.size()>1 ? "s " : " ")
																	 +strange_vertices.join(", ")
																	 +(strange_vertices.size()>1 ? " have " : " has ")
																	 +"no outgoing edges.\n\nDo you still want to run the pipeline?",
																	 QMessageBox::Yes | QMessageBox::No);
				if (ret == QMessageBox::No)
				{
					return false;
				}
			}
			//else
			//{
			//	assume the pipeline was tested in the gui, continue
			//}
		}
				
		return true;
	}
	
	void TOPPASScene::connectVertexSignals(TOPPASVertex* tv)
	{
		connect(tv,SIGNAL(clicked()),this,SLOT(itemClicked()));
		connect(tv,SIGNAL(released()),this,SLOT(itemReleased()));
		connect(tv,SIGNAL(hoveringEdgePosChanged(const QPointF&)),this,SLOT(updateHoveringEdgePos(const QPointF&)));
		connect(tv,SIGNAL(newHoveringEdge(const QPointF&)),this,SLOT(addHoveringEdge(const QPointF&)));
		connect(tv,SIGNAL(finishHoveringEdge()),this,SLOT(finishHoveringEdge()));
		connect(tv,SIGNAL(itemDragged(qreal,qreal)),this,SLOT(moveSelectedItems(qreal,qreal)));
	}
	
	void TOPPASScene::connectToolVertexSignals(TOPPASToolVertex* ttv)
	{
		connect(ttv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(logTOPPOutput(const QString&)));
		connect(ttv, SIGNAL(toolStarted()), this, SLOT(logToolStarted()));
		connect(ttv, SIGNAL(toolFinished()), this, SLOT(logToolFinished()));
		connect(ttv, SIGNAL(toolFailed()), this, SLOT(logToolFailed()));
		connect(ttv, SIGNAL(toolCrashed()), this, SLOT(logToolCrashed()));
		
		connect(ttv,SIGNAL(toolFailed(const QString&)),this,SLOT(pipelineErrorSlot(QString)));
		connect(ttv,SIGNAL(toolCrashed()),this,SLOT(pipelineErrorSlot()));

    connect(ttv, SIGNAL(somethingHasChanged()), this, SLOT(abortPipeline()));
	}


	void TOPPASScene::connectMergerVertexSignals(TOPPASMergerVertex* tmv)
  {
    connect(tmv,SIGNAL(mergeFailed(QString)),this,SLOT(pipelineErrorSlot(QString)));
    connect(tmv, SIGNAL(somethingHasChanged()), this, SLOT(abortPipeline()));
	}

	void TOPPASScene::connectOutputVertexSignals(TOPPASOutputFileListVertex* oflv)
	{
		connect(oflv, SIGNAL(iAmDone()), this, SLOT(checkIfWeAreDone()));
		connect(oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(logOutputFileWritten(const String&)));
	}
	
	void TOPPASScene::connectEdgeSignals(TOPPASEdge* e)
	{
		TOPPASVertex* source = e->getSourceVertex();
		TOPPASVertex* target = e->getTargetVertex();
		connect(source, SIGNAL(somethingHasChanged()), e, SLOT(sourceHasChanged()));
		connect(e, SIGNAL(somethingHasChanged()), target, SLOT(inEdgeHasChanged()));
    connect(e, SIGNAL(somethingHasChanged()), this, SLOT(abortPipeline()));
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
				foreach (const TOPPASResource& res, resource_list)
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
						QMessageBox::warning(0,"Non-unique input node names","Some of the input nodes have the same names. Cannot create resource file.");
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
				foreach (const QString& file, files)
				{
					resource_list << TOPPASResource(file);
				}
				resources.add(key,resource_list);
			}
		}
	}
	
	bool TOPPASScene::refreshParameters()
	{
		bool change = false;
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(*it);
			if (ttv && ttv->refreshParameters())
			{
				change = true;
			}
		}
		
		return change;
	}
  
  bool TOPPASScene::isDryRun() const
  {
    return dry_run_;
  }

  void TOPPASScene::quitWithError()
  {
    exit(1);
  }
	
} //namespace OpenMS

