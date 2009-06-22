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
#include <OpenMS/VISUAL/TOPPASInputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/DIALOGS/TOPPASIOMappingDialog.h>

#include <QtGui/QPainter>
#include <QtGui/QPainterPath>

namespace OpenMS
{	
	
	TOPPASEdge::TOPPASEdge()
		:	QObject(),
			QGraphicsItem(),
			from_(0),
			to_(0),
			hover_pos_(),
			color_(),
			edge_type_(ET_INVALID),
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
			edge_type_(ET_INVALID),
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
			edge_type_(rhs.edge_type_),
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
		edge_type_ = rhs.edge_type_;
		source_out_param_ = rhs.source_out_param_;
		target_in_param_ = rhs.target_in_param_;
		
		setFlag(QGraphicsItem::ItemIsSelectable, true);
		
		return *this;
	}
	
	TOPPASEdge::~TOPPASEdge()
	{
		if (from_)
		{
			from_->removeOutEdge(this);
		}
		if (to_)
		{
			to_->removeInEdge(this);
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
		dialog.exec();
	}
	
	void TOPPASEdge::determineEdgeType()
	{
		bool source_vertex_is_a_tool = false;
		bool source_vertex_is_a_list = false;
		TOPPASVertex* source = getSourceVertex();
		TOPPASVertex* target = getTargetVertex();
		
		if (source == 0 || target == 0)
		{
			return;
		}
		
		if (qobject_cast<TOPPASToolVertex*>(source))
		{
			source_vertex_is_a_tool = true;
		}
		else if (qobject_cast<TOPPASInputFileListVertex*>(source))
		{
			source_vertex_is_a_list = true;
			// fill that
		}
		if (source_vertex_is_a_tool)
		{
			if (qobject_cast<TOPPASToolVertex*>(target))
			{
				edge_type_ = ET_TOOL_TO_TOOL;
			}
			else if (qobject_cast<TOPPASOutputFileListVertex*>(target))
			{
				edge_type_ = ET_TOOL_TO_LIST;
				// here too
			}
			else
			{
				edge_type_ = ET_TOOL_TO_FILE;
				// and so on
			}
		}
		else if (qobject_cast<TOPPASToolVertex*>(target))
		{			
			if (source_vertex_is_a_list)
			{
				edge_type_ = ET_LIST_TO_TOOL;
			}
			else
			{
				edge_type_ = ET_FILE_TO_TOOL;
			}
		}
		else
		{
			edge_type_ = ET_INVALID;
		}
	}
	
	TOPPASEdge::EdgeType TOPPASEdge::getEdgeType()
	{
		return edge_type_;
	}
	
	bool TOPPASEdge::getValidity()
	{
		TOPPASVertex* source = getSourceVertex();
		TOPPASVertex* target = getTargetVertex();
		QVector<TOPPASToolVertex::IOInfo> source_output_files;
		QVector<TOPPASToolVertex::IOInfo> target_input_files;
		StringList source_param_types;
		StringList target_param_types;
		bool source_param_has_list_type = false;
		bool target_param_has_list_type = false;
		bool valid = false;
		
		TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
		if (source_tool && source_out_param_ >= 0)
		{
			source_tool->getOutputFiles(source_output_files);
			TOPPASToolVertex::IOInfo& source_param = source_output_files[(Size)(source_out_param_)];
			source_param_types = source_param.valid_types;
			source_param_has_list_type = source_param.type == TOPPASToolVertex::IOInfo::IOT_LIST;
		}
		TOPPASToolVertex* target_tool = qobject_cast<TOPPASToolVertex*>(target);
		if (target_tool && target_in_param_ >= 0)
		{
			target_tool->getInputFiles(target_input_files);
			TOPPASToolVertex::IOInfo& target_param = target_input_files[(Size)(target_in_param_)];
			target_param_types = target_param.valid_types;
			target_param_has_list_type = target_param.type == TOPPASToolVertex::IOInfo::IOT_LIST;
		}
		if (edge_type_ == ET_FILE_TO_TOOL)
		{
			if (target_param_has_list_type)
			{
				return false;
			}
			else if (target_param_types.empty())
			{
				// no restrictions specified
				valid = true;
			}
			else
			{
				const String& file_name = String(qobject_cast<TOPPASInputFileVertex*>(source)->getFilename());
				String::SizeType extension_start_index = file_name.rfind(".");
				if (extension_start_index != String::npos)
				{
					const String& extension = file_name.substr(extension_start_index+1);
					for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
					{
						if (*it == extension)
						{
							valid = true;
							break;
						}
					}
				}
			}
		}
		else if (edge_type_ == ET_LIST_TO_TOOL)
		{
			if (!target_param_has_list_type)
			{
				return false;
			}
			else if (target_param_types.empty())
			{
				// no restrictions specified
				valid = true;
			}
			else
			{
				const QStringList& file_names = qobject_cast<TOPPASInputFileListVertex*>(source)->getFilenames();
				foreach (const QString& q_file_name, file_names)
				{
					const String& file_name = String(q_file_name);
					String::SizeType extension_start_index = file_name.rfind(".");
					if (extension_start_index != String::npos)
					{
						const String& extension = file_name.substr(extension_start_index+1);
						for (StringList::iterator it = target_param_types.begin(); it != target_param_types.end(); ++it)
						{
							if (*it == extension)
							{
								valid = true;
								break;
							}
						}
					}
				}
			}
		}
		else if (edge_type_ == ET_TOOL_TO_FILE)
		{
			if (source_param_has_list_type)
			{
				return false;
			}
			
			valid = true;
		}
		else if (edge_type_ == ET_TOOL_TO_LIST)
		{
			if (!source_param_has_list_type)
			{
				return false;
			}
			
			valid = true;
		}
		else if (edge_type_ == ET_TOOL_TO_TOOL)
		{
			// check
			valid = true;
		}
		else
		{
			// should not happen
			return false;
		}
		
		return valid;
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
	
} //namespace


