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
#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>

namespace OpenMS
{
	
	TOPPASScene::TOPPASScene(QObject* parent, const String& tmp_path)
		:	QGraphicsScene(parent),
			action_mode_(AM_NEW_EDGE),
			vertices_(),
			edges_(),
			hover_edge_(0),
			potential_target_(0),
			file_name_(),
			tmp_path_(tmp_path)
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
			vertex->setSelected(true);
		}
		foreach (TOPPASEdge* edge, edges_)
		{
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
		TOPPASVertex* sender = qobject_cast<TOPPASVertex*>(QObject::sender());
		if (!sender)
		{
			return;
		}
		bool was_selected = sender->isSelected();
		foreach (QGraphicsItem* item, items())
		{
			item->setSelected(false);
		}
		sender->setSelected(was_selected);
	}
	
	void TOPPASScene::itemDoubleClicked()
	{
		
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
		
		if (target && 
				target != hover_edge_->getSourceVertex() &&
				isEdgeAllowed_(hover_edge_->getSourceVertex(), target))
		{
			hover_edge_->setTargetVertex(target);
			TOPPASVertex* source = hover_edge_->getSourceVertex();
			source->addOutEdge(hover_edge_);
			target->addInEdge(hover_edge_);
			hover_edge_->determineEdgeType();
			hover_edge_->setColor(QColor(255,165,0));
			connect (source, SIGNAL(somethingHasChanged()), hover_edge_, SLOT(sourceHasChanged()));
			connect (hover_edge_, SIGNAL(somethingHasChanged()), target, SLOT(inEdgeHasChanged()));
			
			TOPPASIOMappingDialog dialog(hover_edge_);
			dialog.exec();
			
			hover_edge_->emitChanged();
		}
		else
		{
			edges_.removeAll(hover_edge_);
			removeItem(hover_edge_);
			delete hover_edge_;
			hover_edge_ = 0;
		}
		
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
		
		updateEdgeColors();
	}
	
	bool TOPPASScene::isEdgeAllowed_(TOPPASVertex* u, TOPPASVertex* v)
	{
		if (u == 0 ||
				v == 0 ||
				u == v ||
				// edges leading to input files make no sense:
				qobject_cast<TOPPASInputFileVertex*>(v) ||
				qobject_cast<TOPPASInputFileListVertex*>(v) ||
				// neither do edges coming from output files:
				qobject_cast<TOPPASOutputFileVertex*>(u) ||
				qobject_cast<TOPPASOutputFileListVertex*>(u) ||
				// nor edges from input to output without a tool in between:
				((qobject_cast<TOPPASInputFileVertex*>(u) || qobject_cast<TOPPASInputFileListVertex*>(u))
					&& (qobject_cast<TOPPASOutputFileVertex*>(v) || qobject_cast<TOPPASOutputFileListVertex*>(v))) ||
				// nor multiple incoming edges for a single output file/list node
				((qobject_cast<TOPPASOutputFileVertex*>(v) || qobject_cast<TOPPASOutputFileListVertex*>(v))
					&& v->inEdgesBegin() != v->inEdgesEnd()))
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
		update();
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
		// make sure all output file names are updated, first:
		updateOutputFileNames();
		
		// unset the finished flag for all TOPP tool nodes
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(*it);
			if (tv)
			{
				tv->setFinished(false);
				tv->setStartedHere(false);
			}
		}
		
		// start recursive execution at every output node
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			TOPPASOutputFileVertex* ofv = qobject_cast<TOPPASOutputFileVertex*>(*it);
			if (ofv)
			{
				ofv->startComputation();
				continue;
			}
			
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
			if (oflv)
			{
				oflv->startComputation();
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
			
			TOPPASInputFileVertex* ifv = qobject_cast<TOPPASInputFileVertex*>(tv);
			if (ifv)
			{
				save_param.setValue("vertices:"+id+":type", DataValue("input file"));
				save_param.setValue("vertices:"+id+":file_name", DataValue(String(ifv->getFilename())));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASInputFileListVertex* iflv = qobject_cast<TOPPASInputFileListVertex*>(tv);
			if (iflv)
			{
				const QStringList& files_qt = iflv->getFilenames();
				StringList files;
				foreach (QString file_qt, files_qt)
				{
					files.push_back(String(file_qt));
				}
				save_param.setValue("vertices:"+id+":type", DataValue("input file list"));
				save_param.setValue("vertices:"+id+":file_names", DataValue(files));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASOutputFileVertex* ofv = qobject_cast<TOPPASOutputFileVertex*>(tv);
			if (ofv)
			{
				save_param.setValue("vertices:"+id+":type", DataValue("output file"));
				save_param.setValue("vertices:"+id+":file_name", DataValue(String(ofv->getFilename())));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
			if (oflv)
			{
				const QStringList& files_qt = oflv->getFilenames();
				StringList files;
				foreach (QString file_qt, files_qt)
				{
					files.push_back(String(file_qt));
				}
				save_param.setValue("vertices:"+id+":type", DataValue("output file list"));
				save_param.setValue("vertices:"+id+":file_names", DataValue(files));
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
			if (ttv)
			{
				save_param.setValue("vertices:"+id+":type", DataValue("tool"));
				save_param.setValue("vertices:"+id+":tool_name", DataValue(ttv->getName()));
				save_param.setValue("vertices:"+id+":tool_type", DataValue(ttv->getType()));
				save_param.insert("vertices:"+id+":parameters:", ttv->getParam());
				save_param.setValue("vertices:"+id+":x_pos", DataValue(tv->x()));
				save_param.setValue("vertices:"+id+":y_pos", DataValue(tv->y()));
				continue;
			}
			// NO CODE HERE BECAUSE OF THE "continue"s
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
      if (substrings.back() == "type") // next node (all nodes begin with "type")
      {
      	current_type = (it->value).toString();
      	current_id = substrings[0];
     		Int index = current_id.toInt();
      
				if (current_type == "input file")
				{
					QString file_name = vertices_param.getValue(current_id + ":file_name").toQString();
					TOPPASInputFileVertex* ifv = new TOPPASInputFileVertex(file_name);
					current_vertex = ifv;
				}
				else if (current_type == "input file list")
				{
					StringList file_names = vertices_param.getValue(current_id + ":file_names");
					QStringList file_names_qt;
					for (StringList::const_iterator str_it = file_names.begin(); str_it != file_names.end(); ++str_it)
					{
						file_names_qt.push_back(str_it->toQString());
					}
					TOPPASInputFileListVertex* iflv = new TOPPASInputFileListVertex(file_names_qt);
					current_vertex = iflv;
				}
				else if (current_type == "output file")
				{
					QString file_name = vertices_param.getValue(current_id + ":file_name").toQString();
					TOPPASOutputFileVertex* ofv = new TOPPASOutputFileVertex(file_name);
					current_vertex = ofv;
				}
				else if (current_type == "output file list")
				{
					StringList file_names = vertices_param.getValue(current_id + ":file_names");
					QStringList file_names_qt;
					for (StringList::const_iterator str_it = file_names.begin(); str_it != file_names.end(); ++str_it)
					{
						file_names_qt.push_back(str_it->toQString());
					}
					TOPPASOutputFileListVertex* oflv = new TOPPASOutputFileListVertex(file_names_qt);
					current_vertex = oflv;
				}
				else if (current_type == "tool")
				{
					String tool_name = vertices_param.getValue(current_id + ":tool_name");
					String tool_type = vertices_param.getValue(current_id + ":tool_type");
					Param param_param = vertices_param.copy(current_id + ":parameters:", true);
					TOPPASToolVertex* tv = new TOPPASToolVertex(tool_name, tool_type, tmp_path_);
					tv->setParam(param_param);
					current_vertex = tv;
				}
				else
				{
					std::cerr << "This should not have happened." << std::endl;
				}
				
				if (current_vertex)
				{
					float x = vertices_param.getValue(current_id + ":x_pos");
					float y = vertices_param.getValue(current_id + ":y_pos");
				
					current_vertex->setPos(QPointF(x,y));
					current_vertex->setID((UInt)(current_id.toInt()));
					
					addVertex(current_vertex);
					
					connect(current_vertex,SIGNAL(clicked()),this,SLOT(itemClicked()));
					connect(current_vertex,SIGNAL(doubleClicked()),this,SLOT(itemDoubleClicked()));
					connect(current_vertex,SIGNAL(hoveringEdgePosChanged(const QPointF&)),this,SLOT(updateHoveringEdgePos(const QPointF&)));
					connect(current_vertex,SIGNAL(newHoveringEdge(const QPointF&)),this,SLOT(addHoveringEdge(const QPointF&)));
					connect(current_vertex,SIGNAL(finishHoveringEdge()),this,SLOT(finishHoveringEdge()));
					
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
					std::cerr << "This should not have happened." << std::endl;
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
      	edge->determineEdgeType();	
      	addEdge(edge);
      	connect (tv_1, SIGNAL(somethingHasChanged()), edge, SLOT(sourceHasChanged()));
				connect (edge, SIGNAL(somethingHasChanged()), tv_2, SLOT(inEdgeHasChanged()));
				
      	int source_out_param = (++it)->value;
      	int target_in_param = (++it)->value;
      	edge->setSourceOutParam(source_out_param);
      	edge->setTargetInParam(target_in_param);
      }
    }
    
    // update status of output file list vertices (depends on input)
    foreach (TOPPASVertex* v, vertices_)
    {
    	TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(v);
    	if (oflv)
    	{
    		oflv->updateStatus();
    	}
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
	
	void TOPPASScene::updateOutputFileNames()
	{
		for (VertexIterator it = verticesBegin(); it != verticesEnd(); ++it)
		{
			TOPPASOutputFileVertex* ofv = qobject_cast<TOPPASOutputFileVertex*>(*it);
			if (ofv && ofv->inEdgesBegin() != ofv->inEdgesEnd())
			{
				TOPPASEdge* in_edge = *(ofv->inEdgesBegin());
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(in_edge->getSourceVertex());
				tv->updateOutputFileNames();
				
				continue;
			}
			
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(*it);
			if (oflv && oflv->inEdgesBegin() != oflv->inEdgesEnd())
			{
				TOPPASEdge* in_edge = *(oflv->inEdgesBegin());
				TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>(in_edge->getSourceVertex());
				tv->updateOutputFileNames();
			}
		}
	}
	
} //namespace OpenMS

