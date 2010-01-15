// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASVertex.h>
#include <OpenMS/VISUAL/TOPPASWidget.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASOutputFilesDialog.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtCore/QDir>
#include <QtCore/QFileInfo>
#include <QtGui/QMessageBox>
#include <QtCore/QSet>

namespace OpenMS
{
	
	TOPPASScene::TOPPASScene(QObject* parent, const String& tmp_path, bool gui)
		:	QGraphicsScene(parent),
			action_mode_(AM_NEW_EDGE),
			vertices_(),
			edges_(),
			hover_edge_(0),
			potential_target_(0),
			file_name_(),
			tmp_path_(tmp_path),
			gui_(gui),
			out_dir_(QDir::currentPath()),
			changed_(false),
			running_(false),
			user_specified_out_dir_(false)
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
		
		// unselect all items except for the one under the cursor, but only if no multiple selection
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
					hover_edge_->setColor(Qt::green);
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
		//find backedges via DFS
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

	void TOPPASScene::runPipeline()
	{
		//reset all nodes
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			(*it)->reset(true);
		}
		update(sceneRect());
		
		//check if pipeline OK
		if (!sanityCheck())
		{
			return;
		}
		
		//ask for output directory
		if (!askForOutputDir(true))
		{
			return;
		}
		
		//reset processes
		topp_processes_queue_.clear();
		
		// start at input nodes
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
			if (iflv)
			{
				running_ = true;
				iflv->startPipeline();
			}
		}
	}
	
	void TOPPASScene::store(const String& file)
	{
		Param save_param;
		
		save_param.setValue("info:num_vertices", DataValue(vertices_.size()));
		save_param.setValue("info:num_edges", DataValue(edges_.size()));
		
		// store all vertices (together with all parameters)
		UInt counter = 0;
		foreach (TOPPASVertex* tv, vertices_)
		{
			String id(counter);
			tv->setID(counter);
			counter++;
			
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
			if (iflv)
			{
				// store file names relative to toppas file
				QDir save_dir(File::path(file).toQString());
				const QStringList& files_qt = iflv->getFilenames();
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
				save_param.setValue("vertices:"+id+":list_mode", DataValue("false")); // obsolete, but keep it for compatibility with older versions..
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
		counter = 0;
		foreach (TOPPASEdge* te, edges_)
		{
			if (!(te->getSourceVertex() && te->getTargetVertex()))
			{
				continue;
			}
			
			save_param.setValue("edges:"+String(counter)+":source/target:", DataValue(String(te->getSourceVertex()->getID()) + "/" + String(te->getTargetVertex()->getID())));
			save_param.setValue("edges:"+String(counter)+":source_out_param:", DataValue(te->getSourceOutParam()));
			save_param.setValue("edges:"+String(counter)+":target_in_param:", DataValue(te->getTargetInParam()));
			
			counter++;
		}
		
		//save file
		save_param.store(file);
		changed_ = false;
		file_name_ = file;
	}
	
	void TOPPASScene::load(const String& file)
  {
    Param load_param;
    load_param.load(file);
    Param vertices_param = load_param.copy("vertices:",true);
    Param edges_param = load_param.copy("edges:",true);
    
    String current_type, current_id;
    TOPPASVertex* current_vertex = 0;
    QVector<TOPPASVertex*> vertex_vector;
    vertex_vector.resize((Size)(int)load_param.getValue("info:num_vertices"));
    
    //load all vertices
    for (Param::ParamIterator it = vertices_param.begin(); it != vertices_param.end(); ++it)
    {
      StringList substrings;
      it.getName().split(':', substrings);
      if (substrings.back() == "toppas_type") // next node (all nodes begin with "toppas_type")
      {
      	current_type = (it->value).toString();
      	current_id = substrings[0];
     		Int index = current_id.toInt();
      
				if (current_type == "input file list")
				{
					StringList file_names = vertices_param.getValue(current_id + ":file_names");
					// make file names absolute again
					QDir load_dir(File::path(file).toQString());
					QStringList file_names_qt;
					for (StringList::const_iterator str_it = file_names.begin(); str_it != file_names.end(); ++str_it)
					{
						file_names_qt.push_back(QDir::cleanPath(load_dir.absoluteFilePath(str_it->toQString())));
					}
					TOPPASInputFileListVertex* iflv = new TOPPASInputFileListVertex(file_names_qt);
					current_vertex = iflv;
				}
				else if (current_type == "output file list")
				{
					TOPPASOutputFileListVertex* oflv = new TOPPASOutputFileListVertex();
					connect (oflv, SIGNAL(iAmDone()), this, SLOT(checkIfWeAreDone()));
					if (!gui_)
					{
						connect (oflv, SIGNAL(outputFileWritten(const String&)), this, SLOT(noGuiOutputFileWritten(const String&)));
					}
					current_vertex = oflv;
				}
				else if (current_type == "tool")
				{
					String tool_name = vertices_param.getValue(current_id + ":tool_name");
					String tool_type = vertices_param.getValue(current_id + ":tool_type");
					Param param_param = vertices_param.copy(current_id + ":parameters:", true);
					TOPPASToolVertex* tv = new TOPPASToolVertex(tool_name, tool_type, tmp_path_);
					tv->setParam(param_param);
					connect(tv,SIGNAL(toolStarted()),this,SLOT(setPipelineRunning()));
					connect(tv,SIGNAL(toolFailed()),this,SLOT(pipelineErrorSlot()));
					connect(tv,SIGNAL(toolCrashed()),this,SLOT(pipelineErrorSlot()));
					if (!gui_)
					{
						connect (tv, SIGNAL(toppOutputReady(const QString&)), this, SLOT(noGuiTOPPOutput(const QString&)));
						connect (tv, SIGNAL(toolStarted()), this, SLOT(noGuiToolStarted()));
						connect (tv, SIGNAL(toolFinished()), this, SLOT(noGuiToolFinished()));
						connect (tv, SIGNAL(toolFailed()), this, SLOT(noGuiToolFailed()));
						connect (tv, SIGNAL(toolCrashed()), this, SLOT(noGuiToolCrashed()));
					}
					
					current_vertex = tv;
				}
				else if (current_type == "merger")
				{
					TOPPASMergerVertex* mv = new TOPPASMergerVertex();
					if (vertices_param.exists(current_id + ":round_based"))
					{
						String rb = vertices_param.getValue(current_id + ":round_based");
						mv->setRoundBasedMode(rb == "true" ? true : false);
					}
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
					current_vertex->setID((UInt)(current_id.toInt()));
					
					addVertex(current_vertex);
					
					connect(current_vertex,SIGNAL(clicked()),this,SLOT(itemClicked()));
					connect(current_vertex,SIGNAL(released()),this,SLOT(itemReleased()));
					connect(current_vertex,SIGNAL(hoveringEdgePosChanged(const QPointF&)),this,SLOT(updateHoveringEdgePos(const QPointF&)));
					connect(current_vertex,SIGNAL(newHoveringEdge(const QPointF&)),this,SLOT(addHoveringEdge(const QPointF&)));
					connect(current_vertex,SIGNAL(finishHoveringEdge()),this,SLOT(finishHoveringEdge()));
					connect(current_vertex,SIGNAL(itemDragged(qreal,qreal)),this,SLOT(moveSelectedItems(qreal,qreal)));
					
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
      	addEdge(edge);
      	connect (tv_1, SIGNAL(somethingHasChanged()), edge, SLOT(sourceHasChanged()));
				connect (edge, SIGNAL(somethingHasChanged()), tv_2, SLOT(inEdgeHasChanged()));
				
      	int source_out_param = (++it)->value;
      	int target_in_param = (++it)->value;
      	edge->setSourceOutParam(source_out_param);
      	edge->setTargetInParam(target_in_param);
      }
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
		
		file_name_ = file;
		
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
		{
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
			if (oflv && !oflv->isFinished())
			{
				return;
			}
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(*it);
			if (mv && !mv->mergeComplete())
			{
				return;
			}
		}
		
		running_ = false;
		emit entirePipelineFinished();
	}
	
	void TOPPASScene::pipelineErrorSlot()
	{
		running_ = false;
		emit pipelineExecutionFailed();
	}
	
	void TOPPASScene::noGuiTOPPOutput(const QString& out)
	{
		TOPPASToolVertex* sender = qobject_cast<TOPPASToolVertex*>(QObject::sender());
		if (!sender)
		{
			return;
		}
		String tool = sender->getName();
		if (sender->getType() != "")
		{
			tool += " ("+sender->getType()+")";
		}
		std::cout	<< std::endl
							<< tool
							<< std::endl
							<< String(out)
							<< std::endl;
	}
	
	void TOPPASScene::noGuiToolStarted()
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
			
			std::cout	<< std::endl
								<< text
								<< std::endl;
		}
	}
	
	void TOPPASScene::noGuiToolFinished()
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
			
			std::cout	<< std::endl
								<< text
								<< std::endl;
		}
	}
	
	void TOPPASScene::noGuiToolFailed()
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
			
			std::cout	<< std::endl
								<< text
								<< std::endl;
		}
	}
	
	void TOPPASScene::noGuiToolCrashed()
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
			
			std::cout	<< std::endl
								<< text
								<< std::endl;
		}
	}
	
	void TOPPASScene::noGuiOutputFileWritten(const String& file)
	{
		String text = "Output file '"+file+"' written.";
		std::cout	<< std::endl
							<< text
							<< std::endl;
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
				// ignore input vertices (need no tmp directory with number)
				TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(*it);
				if (iflv || (*it)->isTopoSortMarked())
				{
					continue;
				}
				some_vertex_not_finished = true;
				bool has_predecessors = false;
				for (TOPPASVertex::EdgeIterator e_it = (*it)->inEdgesBegin(); e_it != (*it)->inEdgesEnd(); ++e_it)
				{
					TOPPASVertex* v = (*e_it)->getSourceVertex();
					TOPPASInputFileListVertex* v_iflv = qobject_cast<TOPPASInputFileListVertex*>(v);
					if (!v_iflv && !(v->isTopoSortMarked()))
					{
						has_predecessors = true;
						break;
					}
				}
				if (!has_predecessors)
				{
					(*it)->setTopoSortMarked(true);
					(*it)->setTopoNr(topo_counter++);
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
	
	void TOPPASScene::setOutDir(const QString& dir)
	{
		out_dir_ = dir;
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
		
		changed_ = true;
	}

	void TOPPASScene::snapToGrid()
	{
		int grid_step = 20;

		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
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
		changed_ = b;
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
	
	bool TOPPASScene::askForOutputDir(bool always_ask)
	{
		if (gui_)
		{
			if (always_ask || !user_specified_out_dir_)
			{
				TOPPASOutputFilesDialog tofd(out_dir_);
				if (tofd.exec())
				{
					out_dir_ = tofd.getDirectory();
					user_specified_out_dir_ = true;
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
		
		if (clicked_item == 0)
		{
			return;
		}
		
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
		
		QMenu menu;
		QList<QSet<QString> > all_actions;
		
		if (found_tool)
		{
			QSet<QString> tool_actions;
			tool_actions.insert("Edit parameters");
			tool_actions.insert("Resume");
			tool_actions.insert("Open files in TOPPView");
			tool_actions.insert("Remove");
			all_actions.push_back(tool_actions);
		}
		
		if (found_input)
		{
			QSet<QString> input_actions;
			input_actions.insert("Change files");
			input_actions.insert("Open files in TOPPView");
			input_actions.insert("Remove");
			all_actions.push_back(input_actions);
		}
		
		if (found_output)
		{
			QSet<QString> output_actions;
			output_actions.insert("Open files in TOPPView");
			output_actions.insert("Remove");
			all_actions.push_back(output_actions);
		}
		
		if (found_edge)
		{
			QSet<QString> edge_actions;
			edge_actions.insert("Edit I/O mapping");
			edge_actions.insert("Remove");
			all_actions.push_back(edge_actions);
		}
		
		if (found_merger)
		{
			QSet<QString> merger_actions;
			merger_actions.insert("Remove");
			merger_actions.insert("Change mode");
			all_actions.push_back(merger_actions);
		}
		
		QSet<QString> supported_actions_set = all_actions.first();
		foreach (const QSet<QString>& action_set, all_actions)
		{
			supported_actions_set.intersect(action_set);
		}
		QList<QString> supported_actions = supported_actions_set.toList();
		
		foreach (const QString& supported_action, supported_actions)
		{
			QAction* new_action = menu.addAction(supported_action);
			if (supported_action == "Resume" && disable_resume)
			{
				new_action->setEnabled(false);
			}
		}
		
		// ------ execute action on all selected items ------
		
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
			
			foreach (QGraphicsItem* gi, selectedItems())
			{
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
							ttv->runToolIfInputReady();
						}
					}
					else if (text == "Open files in TOPPView")
					{
						ttv->openInTOPPView();
					}
					
					continue;
				}
				
				TOPPASInputFileListVertex* ifv = dynamic_cast<TOPPASInputFileListVertex*>(gi);
				if (ifv)
				{
					if (text == "Open files in TOPPView")
					{
						ifv->openInTOPPView();
					}
					else if (text == "Change files")
					{
						ifv->showFilesDialog();
					}
					
					continue;
				}
				
				TOPPASOutputFileListVertex* ofv = dynamic_cast<TOPPASOutputFileListVertex*>(gi);
				if (ofv)
				{
					if (text == "Open files in TOPPView")
					{
						ofv->openInTOPPView();
					}
					
					continue;
				}
				
				TOPPASMergerVertex* mv = dynamic_cast<TOPPASMergerVertex*>(gi);
				if (mv)
				{
					if (text == "Change mode")
					{
						mv->setRoundBasedMode(!mv->roundBasedMode());
						mv->update(mv->boundingRect());
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
		
		// run first process
		if (topp_processes_queue_.size() == 1)
		{
			const TOPPProcess& tp = topp_processes_queue_.first();
			tp.proc->start(tp.command, tp.args);
		}
	}
	
	void TOPPASScene::runNextProcess()
	{
		if (topp_processes_queue_.empty())
		{
			return;
		}
		
		topp_processes_queue_.removeFirst();
		if (!topp_processes_queue_.empty())
		{
			const TOPPProcess& tp = topp_processes_queue_.first();
			tp.proc->start(tp.command, tp.args);
		}
	}
	
	bool TOPPASScene::sanityCheck()
	{
		QStringList strange_vertices;
		
		// ----- are there any input nodes and are files specified? ----
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
		foreach (TOPPASInputFileListVertex* iflv, input_nodes)
		{
			if (iflv->getFilenames().empty())
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
		
		// ----- are there nodes without parents (besides input nodes)? -----
		strange_vertices.clear();
		foreach (TOPPASVertex* tv, vertices_)
		{
			if (qobject_cast<TOPPASInputFileListVertex*>(tv))
			{
				continue;
			}
			if (tv->inEdgesBegin() == tv->inEdgesEnd())
			{
				strange_vertices << QString::number(tv->getTopoNr());
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
		
		// ----- are there mergers with unequal input list lengths (per merge round and over entire run)? -----
		QStringList unequal_per_round;
		QStringList unequal_over_entire_run;
		
		foreach (TOPPASVertex* tv, vertices_)
		{
			if (qobject_cast<TOPPASInputFileListVertex*>(tv))
			{
				tv->checkListLengths(unequal_per_round, unequal_over_entire_run);
			}
		}
		
		if (!unequal_per_round.empty() || !unequal_over_entire_run.empty())
		{
			if (gui_)
			{
				QString message("");
				if (!unequal_per_round.empty())
				{
					message = QString("Node")
										+(unequal_per_round.size()>1 ? "s " : " ")
										+unequal_per_round.join(", ")
										+(unequal_per_round.size()>1 ? " have " : " has ")
										+"unequal input list lengths. Some files will not be processed.\n\n";
				}
				foreach(const QString& str, unequal_per_round)
				{
					unequal_over_entire_run.removeAll(str);
				}
				if (!unequal_over_entire_run.empty())
				{
					message += QString("Merger")
										+(unequal_over_entire_run.size()>1 ? "s " : " ")
										+unequal_over_entire_run.join(", ")+":\n"
										+"The overall number of files to be merged is not the same "
										+"for all incoming edges. This either means that some files "
										+"will not be merged or that one and the same file will be "
										+"merged several times.\n\n";
				}
				message += "Do you still want to continue?";
				
				QMessageBox::StandardButton ret;
				ret = QMessageBox::warning(views().first(), "Unequal input list lengths",
																	 message,
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
	
} //namespace OpenMS

