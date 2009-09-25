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

#include <OpenMS/VISUAL/TOPPASEdge.h>
#include <OpenMS/VISUAL/TOPPASScene.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/TOPPASMergerVertex.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>
#include <QtGui/QMessageBox>
#include <QtGui/QMenu>


namespace OpenMS
{	
	
	TOPPASEdge::TOPPASEdge()
		:	QObject(),
			QGraphicsItem(),
			from_(0),
			to_(0),
			hover_pos_(),
			color_(),
			source_out_param_(-1),
			target_in_param_(-1)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge::TOPPASEdge(TOPPASVertex* from, const QPointF& hover_pos)
		:	QObject(),
			QGraphicsItem(),
			from_(from),
			to_(0),
			hover_pos_(hover_pos),
			color_(),
			source_out_param_(-1),
			target_in_param_(-1)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge::TOPPASEdge(const TOPPASEdge& rhs)
		:	QObject(),
			QGraphicsItem(),
			from_(rhs.from_),
			to_(rhs.to_),
			hover_pos_(rhs.hover_pos_),
			color_(rhs.color_),
			source_out_param_(rhs.source_out_param_),
			target_in_param_(rhs.target_in_param_)
	{
		setFlag(QGraphicsItem::ItemIsSelectable, true);
	}
	
	TOPPASEdge& TOPPASEdge::operator= (const TOPPASEdge& rhs)
	{
		from_ = rhs.from_;
		to_ = rhs.to_;
		hover_pos_ = rhs.hover_pos_;
		color_ = rhs.color_;
		source_out_param_ = rhs.source_out_param_;
		target_in_param_ = rhs.target_in_param_;
		
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		
		return *this;
	}
	
	TOPPASEdge::~TOPPASEdge()
	{
		// notify our childs that we are dying ;)
		emit somethingHasChanged();
		
		if (from_)
		{
			from_->removeOutEdge(this);
			disconnect (from_, SIGNAL(somethingHasChanged()), this, SLOT(sourceHasChanged()));
		}
		if (to_)
		{
			to_->removeInEdge(this);
			disconnect (this, SIGNAL(somethingHasChanged()), to_, SLOT(inEdgeHasChanged()));
		}
	}
	
	QRectF TOPPASEdge::boundingRect() const
	{
		qreal min_x = startPos().x() < endPos().x() ? startPos().x() : endPos().x();
		qreal min_y = startPos().y() < endPos().y() ? startPos().y() : endPos().y();
		qreal max_x = startPos().x() > endPos().x() ? startPos().x() : endPos().x();
		qreal max_y = startPos().y() > endPos().y() ? startPos().y() : endPos().y();
		
		return QRectF(QPointF(min_x-11.0,min_y-11.0), QPointF(max_x+11.0,max_y+11.0));
	}
	
	QPainterPath TOPPASEdge::shape () const
	{
		QPainterPath shape_1;				
		shape_1.moveTo(startPos() + QPointF(-10,-10));
		shape_1.lineTo(endPos() + QPointF(-10,-10));
		shape_1.lineTo(endPos() + QPointF(10,10));
		shape_1.lineTo(startPos() + QPointF(10,10));
		shape_1.closeSubpath();
		
		QPainterPath shape_2;
		shape_2.moveTo(startPos() + QPointF(-10,10));
		shape_2.lineTo(endPos() + QPointF(-10,10));
		shape_2.lineTo(endPos() + QPointF(10,-10));
		shape_2.lineTo(startPos() + QPointF(10,-10));
		shape_2.closeSubpath();
				
		return shape_1.united(shape_2);
	}
	
	void TOPPASEdge::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		painter->setBrush(Qt::white);

		QPen pen(color_);
		if (isSelected())
		{
			pen.setWidth(2);
		}
		painter->setPen(pen);
		
		painter->drawLine(startPos(),endPos());
		
		// draw arrow head
		QPointF delta = endPos() - startPos();
		qreal angle;
		if (delta.x() == 0.0)
		{
			angle = endPos().y() > startPos().y() ? 90.0 : 270.0;
		}
		else
		{
			angle = delta.y() / delta.x();
			angle = std::atan(angle);
			angle = (angle / 3.14159265) * 180.0;
			if (delta.x() < 0.0)
			{
				angle += 180;
			}
		}
		
		painter->save();
		painter->translate(endPos());
		painter->rotate(angle);
		QPainterPath path;
		path.moveTo(QPointF(0,0));
		path.lineTo(QPointF(-10,4));
		path.lineTo(QPointF(-10,-4));
		path.closeSubpath();
		painter->drawPath(path);
		painter->restore();
	}
	
	QPointF TOPPASEdge::startPos() const
	{
		if (from_)
		{
			return mapFromScene(from_->scenePos());
		}
		else
		{
			return QPointF();
		}
	}
	
	QPointF TOPPASEdge::endPos() const
	{
		QPointF position;
		
		if (!to_)
		{
			// we do not have a target vertex yet
			position = hover_pos_;
		}
		else
		{
			// we have a target node --> line should end at its border
			
			QList<QPointF> point_list;
			
			QPointF target_pos = mapFromScene(to_->scenePos());
			QRectF target_boundings = mapFromItem(to_, to_->shape()).boundingRect();
			QPointF delta = target_pos - startPos();
			qreal slope;
			if (delta.x() == 0)
			{
				slope = std::numeric_limits<double>::infinity();
			}
			else
			{
				slope = delta.y() / delta.x();
			}
			
			qreal y_1 = startPos().y() + slope * (target_boundings.left() - startPos().x());
			qreal y_2 = startPos().y() + slope * (target_boundings.right() - startPos().x());
			
			slope = 1.0 / slope;
			
			qreal x_3 = startPos().x() + slope * (target_boundings.top() - startPos().y());
			qreal x_4 = startPos().x() + slope * (target_boundings.bottom() - startPos().y());
			
			if (y_1 <= target_boundings.bottom() && y_1 >= target_boundings.top())
			{
				point_list.push_back(QPointF(target_boundings.left(), y_1));
			}
			if (y_2 <= target_boundings.bottom() && y_2 >= target_boundings.top())
			{
				point_list.push_back(QPointF(target_boundings.right(), y_2));
			}
			if (x_3 <= target_boundings.right() && x_3 >= target_boundings.left())
			{
				point_list.push_back(QPointF(x_3, target_boundings.top()));
			}
			if (x_4 <= target_boundings.right() && x_4 >= target_boundings.left())
			{
				point_list.push_back(QPointF(x_4, target_boundings.bottom()));
			}
			
			position = nearestPoint_(startPos(), point_list);
		}
		
		return position;
	}
	
	void TOPPASEdge::setHoverPos(const QPointF& pos)
	{
		prepareResize();
		hover_pos_ = pos;
		update();
	}
	
	void TOPPASEdge::setTargetVertex(TOPPASVertex* tv)
	{
		to_ = tv;
	}
	
	void TOPPASEdge::setSourceVertex(TOPPASVertex* tv)
	{
		from_ = tv;
	}
	
	TOPPASVertex* TOPPASEdge::getSourceVertex()
	{
		return from_;
	}
	
	TOPPASVertex* TOPPASEdge::getTargetVertex()
	{
		return to_;
	}
	
	void TOPPASEdge::prepareResize()
	{
		prepareGeometryChange();
	}
	
	QPointF TOPPASEdge::nearestPoint_(const QPointF& origin, const QList<QPointF>& list) const
	{
		if (list.empty())
		{
			return QPointF();
		}
		QPointF nearest = list.first();
		qreal min_distance = std::numeric_limits<double>::max();
		
		for (QList<QPointF>::const_iterator it = list.begin(); it != list.end(); ++it)
		{
			qreal sqr_distance = (it->x() - origin.x()) * (it->x() - origin.x()) +
												(it->y() - origin.y()) * (it->y() - origin.y());
			if (sqr_distance < min_distance)
			{
				min_distance = sqr_distance;
				nearest = *it;
			}
		}
		
		return nearest;
	}
	
	void TOPPASEdge::setColor(const QColor& color)
	{
		color_ = color;
	}

	void TOPPASEdge::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
		TOPPASIOMappingDialog dialog(this);
		if (dialog.exec())
		{
			emit somethingHasChanged();
		}
	}

	TOPPASEdge::EdgeStatus TOPPASEdge::getToolToolStatus_(TOPPASToolVertex* source_tool, int source_param_index, TOPPASToolVertex* target_tool, int target_param_index)
	{
		QVector<TOPPASToolVertex::IOInfo> source_output_files;
		QVector<TOPPASToolVertex::IOInfo> target_input_files;
		source_tool->getOutputParameters(source_output_files);
		const TOPPASToolVertex::IOInfo& source_param = source_output_files[source_param_index];
		StringList source_param_types = source_param.valid_types;
		target_tool->getInputParameters(target_input_files);
		const TOPPASToolVertex::IOInfo& target_param = target_input_files[target_param_index];
		StringList target_param_types = target_param.valid_types;
		
		if (source_param_types.size() == 0 || target_param_types.size() == 0)
		{
			// no type specified --> allow edge
			return ES_VALID;
		}
		else
		{
			// check file type compatibility
			bool found_match = false;
			for (StringList::iterator s_it = source_param_types.begin(); s_it != source_param_types.end(); ++s_it)
			{
				String ext_1 = *s_it;
				ext_1.toLower();
				found_match = false;
				for (StringList::iterator t_it = target_param_types.begin(); t_it != target_param_types.end(); ++t_it)
				{
					String ext_2 = *t_it;
					ext_2.toLower();
					if (ext_1 == ext_2)
					{
						found_match = true;
						break;
					}
				}
				if (found_match)
				{
					break;
				}
			}

			if (!found_match)
			{
				return ES_FILE_EXT_MISMATCH;
			}

			return ES_VALID;
		}
	}

	TOPPASEdge::EdgeStatus TOPPASEdge::getListToolStatus_(TOPPASInputFileListVertex* source_input_list, TOPPASToolVertex* target_tool, int target_param_index)
	{
		const QStringList& file_names = source_input_list->getFilenames();
		if (file_names.empty())
		{
			// file names are not specified yet
			return ES_NOT_READY_YET;
		}

		QVector<TOPPASToolVertex::IOInfo> target_input_files;
		target_tool->getInputParameters(target_input_files);
		const TOPPASToolVertex::IOInfo& target_param = target_input_files[target_param_index];
		StringList target_param_types = target_param.valid_types;

		if (target_param_types.empty())
		{
			// no file types specified --> allow
			return ES_VALID;
		}

		// check file type compatibility
		bool type_mismatch = false;
		foreach (const QString& q_file_name, file_names)
		{
			type_mismatch = true;
			const String& file_name = String(q_file_name);
			String::SizeType extension_start_index = file_name.rfind(".");
			if (extension_start_index != String::npos)
			{
				String extension = file_name.substr(extension_start_index+1);
				extension.toLower();
				for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
				{
					String other_ext = *it;
					other_ext.toLower();
					if (extension == other_ext)
					{
						type_mismatch = false;
						break;
					}
				}
			}

			if (type_mismatch)
			{
				return ES_FILE_EXT_MISMATCH;
			}
		}

		// all file types ok
		return ES_VALID;
	}
		
	TOPPASEdge::EdgeStatus TOPPASEdge::getEdgeStatus()
	{
		TOPPASVertex* source = getSourceVertex();
		TOPPASVertex* target = getTargetVertex();
		TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
		TOPPASMergerVertex* target_merger = qobject_cast<TOPPASMergerVertex*>(target);
		TOPPASInputFileListVertex* source_input_list = qobject_cast<TOPPASInputFileListVertex*>(source);
		TOPPASOutputFileListVertex* target_output_list = qobject_cast<TOPPASOutputFileListVertex*>(target);
		TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
		TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);
		
		if (target_output_list)
		// edges to output vertices are always valid (if finishHoveringEdge_() allowed to construct them in the first place)
		{
			return ES_VALID;
		}

		if (source_tool && source_out_param_ < 0)
		{
			return ES_NO_SOURCE_PARAM;
		}

		if (target_tool && target_in_param_ < 0)
		{
			return ES_NO_TARGET_PARAM;
		}
		
		if (source_tool && target_tool)
		{
			return getToolToolStatus_(source_tool, source_out_param_, target_tool, target_in_param_);
		}
		
		if (source_input_list && target_tool)
		{
			return getListToolStatus_(source_input_list, target_tool, target_in_param_);
		}

		if (source_merger)
		{
			//check compatibility of source with all target nodes of merger
			for (TOPPASVertex::EdgeIterator e_it = source_merger->inEdgesBegin(); e_it != source_merger->inEdgesEnd(); ++e_it)
			{
				TOPPASEdge* merger_in_edge = *e_it;
				TOPPASToolVertex* merger_in_tool = qobject_cast<TOPPASToolVertex*>(merger_in_edge->getSourceVertex());
				TOPPASInputFileListVertex* merger_in_list = qobject_cast<TOPPASInputFileListVertex*>(merger_in_edge->getSourceVertex());
				
				if (merger_in_tool && target_tool)
				{
					EdgeStatus es = getToolToolStatus_(merger_in_tool, merger_in_edge->getSourceOutParam(), target_tool, target_in_param_);
					if (es != ES_VALID)
					{
						return es;
					}
				}
				else if (merger_in_list)
				{
					if (target_output_list)
					{
						// [input] -> [merger] -> [output]) makes no sense
						return ES_MERGER_WITHOUT_TOOL;
					}
					else if (target_tool)
					{
						EdgeStatus es = getListToolStatus_(merger_in_list, target_tool, target_in_param_);
						if (es != ES_VALID)
						{
							return es;
						}
					}
				}
			}
			// no incompatible merger target found
			return ES_VALID;
		}

		if (target_merger)
		{
			//check compatibility of source with all target nodes of merger
			for (TOPPASVertex::EdgeIterator e_it = target_merger->outEdgesBegin(); e_it != target_merger->outEdgesEnd(); ++e_it)
			{
				TOPPASEdge* merger_out_edge = *e_it;
				TOPPASToolVertex* merger_out_tool = qobject_cast<TOPPASToolVertex*>(merger_out_edge->getTargetVertex());
				
				if (source_tool)
				{
					if (merger_out_tool)
					{
						EdgeStatus es = getToolToolStatus_(source_tool, source_out_param_, merger_out_tool, merger_out_edge->getTargetInParam());
						if (es != ES_VALID)
						{
							return es;
						}
					}
				}
				else if (source_input_list)
				{
					if (merger_out_tool)
					{
						EdgeStatus es = getListToolStatus_(source_input_list, merger_out_tool, merger_out_edge->getTargetInParam());
						if (es != ES_VALID)
						{
							return es;
						}
					}	
					else
					// [input] -> [merger] -> [output]) makes no sense
					{
						return ES_MERGER_WITHOUT_TOOL;
					}
				}
			}
			// no incompatible merger target found
			return ES_VALID;
		}

		return ES_UNKNOWN;
	}
	
	void TOPPASEdge::setSourceOutParam(int out)
	{
		source_out_param_ = out;
	}
	
	int TOPPASEdge::getSourceOutParam()
	{
		return source_out_param_;
	}
	
	void TOPPASEdge::setTargetInParam(int in)
	{
		target_in_param_ = in;
	}
	
	int TOPPASEdge::getTargetInParam()
	{
		return target_in_param_;
	}
	
	void TOPPASEdge::updateColor()
	{
		EdgeStatus es = getEdgeStatus();
		if (es == ES_VALID)
		{
			setColor(Qt::green);
		}
		else if (es == ES_NOT_READY_YET)
		{
			setColor(QColor(255,165,0));
		}
		else
		{
			setColor(Qt::red);
		}
		update(boundingRect());
	}
	
	void TOPPASEdge::sourceHasChanged()
	{
		emit somethingHasChanged();
	}
	
	void TOPPASEdge::emitChanged()
	{
		emit somethingHasChanged();
	}
	
	void TOPPASEdge::contextMenuEvent(QGraphicsSceneContextMenuEvent* event)
	{
		TOPPASScene* ts = qobject_cast<TOPPASScene*>(scene());
		ts->unselectAll();
		setSelected(true);
		
		QMenu menu;
		menu.addAction("Edit I/O mapping");
		menu.addAction("Remove");
		
		QAction* selected_action = menu.exec(event->screenPos());
		if (selected_action)
		{
			QString text = selected_action->text();
			if (text == "Edit I/O mapping")
			{
				TOPPASIOMappingDialog dialog(this);
				if (dialog.exec())
				{
					emit somethingHasChanged();
				}
			}
			else if (text == "Remove")
			{
				ts->removeSelected();
			}
			event->accept();
		}
		else
		{
			event->ignore();	
		}
	}
	
} //namespace


