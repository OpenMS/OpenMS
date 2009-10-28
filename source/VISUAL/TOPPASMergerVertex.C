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

#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

namespace OpenMS
{
	TOPPASMergerVertex::TOPPASMergerVertex()
		:	TOPPASVertex(),
			round_based_mode_(true),
			merge_complete_(false),
			merge_counter_(0)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASMergerVertex::TOPPASMergerVertex(const TOPPASMergerVertex& rhs)
		:	TOPPASVertex(rhs),
			round_based_mode_(rhs.round_based_mode_),
			merge_complete_(rhs.merge_complete_),
			merge_counter_(rhs.merge_counter_)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASMergerVertex::~TOPPASMergerVertex()
	{
	
	}
	
	TOPPASMergerVertex& TOPPASMergerVertex::operator= (const TOPPASMergerVertex& rhs)
	{
		TOPPASVertex::operator=(rhs);
		round_based_mode_ = rhs.round_based_mode_;
		merge_complete_ = rhs.merge_complete_;
		merge_counter_ = rhs.merge_counter_;
		
		return *this;
	}
	
	void TOPPASMergerVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
	}
	
	QStringList TOPPASMergerVertex::getCurrentOutputList()
	{
		QStringList out_files;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* source = (*it)->getSourceVertex();
			int param_index = (*it)->getSourceOutParam();
			
			TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
			if (source_tool)
			{
				if (round_based_mode_)
				{
					const QVector<QStringList>& output_files = source_tool->getCurrentOutputFileNames();
					const QStringList& file_names = output_files[param_index];
					out_files << file_names[merge_counter_];
				}
				else
				{
					const QVector<QStringList>& output_files = source_tool->getAllWrittenOutputFileNames();
					const QStringList& file_names = output_files[param_index];
					out_files << file_names;
				}
				
				continue;
			}
			TOPPASInputFileListVertex* source_list = qobject_cast<TOPPASInputFileListVertex*>(source);
			if (source_list)
			{
				if (round_based_mode_)
				{
					out_files << source_list->getFilenames()[merge_counter_];
				}
				else
				{
					out_files << source_list->getFilenames();
				}
				continue;
			}
			TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
			if (source_merger)
			{
				if (round_based_mode_)
				{
					out_files << source_merger->getCurrentOutputList()[merge_counter_];
				}
				else
				{
					out_files << source_merger->getAllCollectedFiles_();
				}
				continue;
			}
		}

		return out_files;
	}
	
	QStringList TOPPASMergerVertex::getAllCollectedFiles_()
	{
		QStringList out_files;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* source = (*it)->getSourceVertex();
			
			TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
			if (source_tool)
			{
				const QVector<QStringList>& output_files = source_tool->getAllWrittenOutputFileNames();
				int param_index = (*it)->getSourceOutParam();
				const QStringList& file_names = output_files[param_index];
				out_files << file_names;
				continue;
			}
			TOPPASInputFileListVertex* source_list = qobject_cast<TOPPASInputFileListVertex*>(source);
			if (source_list)
			{
				out_files << source_list->getFilenames();
				continue;
			}
			TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
			if (source_merger)
			{
				out_files << source_merger->getAllCollectedFiles_();
				continue;
			}
		}
		
		return out_files;
	}

	bool TOPPASMergerVertex::allInputsReady()
	{
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
		  TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
		  if (tv && !tv->isFinished())
		  {
		    // some tool that we depend on has not finished execution yet --> do not start yet
		    return false;
		  }
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
			if (mv && !mv->allInputsReady())
			{
				// some input of a merger that we depend on not finished yet
				return false;
			}
		}
		return true;
	}

	void TOPPASMergerVertex::forwardPipelineExecution()
	{
		if (!allInputsReady())
		{
			return;
		}
		
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			TOPPASToolVertex* ttv = qobject_cast<TOPPASToolVertex*>(tv);
			if (ttv)
			{
				ttv->runToolIfInputReady();
				continue;
			}
			TOPPASOutputFileListVertex* oflv = qobject_cast<TOPPASOutputFileListVertex*>(tv);
			if (oflv)
			{
				oflv->finish();
				continue;
			}
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>(tv);
			if (mv)
			{
				mv->forwardPipelineExecution();
				continue;
			}
		}
		
		++merge_counter_;
	}

	void TOPPASMergerVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		QPen pen(pen_color_, 1, Qt::SolidLine, Qt::FlatCap, Qt::MiterJoin);
		if (isSelected())
		{
			pen.setWidth(2);
			painter->setBrush(brush_color_.darker(130));
			pen.setColor(Qt::darkBlue);
		}
		else
		{
			painter->setBrush(brush_color_);
		}
		painter->setPen(pen);
		
		QPainterPath path;
		path.addRoundRect(-40.0, -40.0, 80.0, 80.0, 20, 20);
 		painter->drawPath(path);
 		
 		pen.setColor(pen_color_);
 		painter->setPen(pen);
		QString text = round_based_mode_ ? "Merge" : "Merge all";
		QRectF text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
		painter->drawText(-(int)(text_boundings.width()/2.0), (int)(text_boundings.height()/4.0), text);
	}
	
	QRectF TOPPASMergerVertex::boundingRect() const
	{
		return QRectF(-41,-41,82,82);
	}
	
	QPainterPath TOPPASMergerVertex::shape () const
	{
		QPainterPath shape;
		shape.addRoundRect(-41.0, -41.0, 82.0, 82.0, 20, 20);
		return shape;
	}
	
	bool TOPPASMergerVertex::mergeComplete()
	{
		return merge_complete_;
	}
	
	void TOPPASMergerVertex::setMergeComplete(bool b)
	{
		merge_complete_ = b;
	}
	
	bool TOPPASMergerVertex::roundBasedMode()
	{
		return round_based_mode_;
	}
	
	void TOPPASMergerVertex::setRoundBasedMode(bool b)
	{
		round_based_mode_ = b;
	}
	
	int TOPPASMergerVertex::minInputListLength_()
	{
		int min_length = std::numeric_limits<int>::max();
		
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* source = (*it)->getSourceVertex();
			TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
			if (source_tool)
			{
				const QVector<QStringList>& output_files = source_tool->getCurrentOutputFileNames();
				int param_index = (*it)->getSourceOutParam();
				const QStringList& file_names = output_files[param_index];
				
				if (file_names.size() < min_length)
				{
					min_length = file_names.size();
				}
				continue;
			}
			TOPPASInputFileListVertex* source_list = qobject_cast<TOPPASInputFileListVertex*>(source);
			if (source_list)
			{
				if (source_list->getFilenames().size() < min_length)
				{
					min_length = source_list->getFilenames().size();
				}
				continue;
			}
			TOPPASMergerVertex* source_merger = qobject_cast<TOPPASMergerVertex*>(source);
			if (source_merger)
			{
				if (source_merger->getCurrentOutputList().size() < min_length)
				{
					min_length = source_merger->getCurrentOutputList().size();
				}
				continue;
			}
		}
		
		return min_length;
	}
	
	int TOPPASMergerVertex::numIterations_()
	{
		if (round_based_mode_)
		{
			return minInputListLength_();
		}
		else
		{
			return 1;
		}
	}
	
	void TOPPASMergerVertex::propagateDownwardsMergeComplete()
	{
		/*	We do not have to forward this (TODO: or do we?) (unlike in the base class implementation).
				Proceed in pipeline if this merger is a "collecting" merger (instead of round based)	*/
		
		
	}
	
	void TOPPASMergerVertex::propagateUpwardsSubtreeFinished()
	{
		/*	Proceed with next merging round, if this merger is still running and
				if all children have finished already. If merge is complete, propagate
				this upwards. */
		
		// all children's subtrees finished?
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* target = (*it)->getTargetVertex();
			if (!target->isSubtreeFinished())
			{
				return;
			}
		}
		
		if (merge_counter_ < numIterations_())
		{
			resetSubtree(false);
			forwardPipelineExecution();
		}
		else
		{
			merge_complete_ = true;
			//actually propagate upwards
			TOPPASVertex::propagateUpwardsSubtreeFinished();
		}
	}
	
	void TOPPASMergerVertex::reset(bool /*reset_all_files*/)
	{
		TOPPASVertex::reset();
		merge_complete_ = false;
		merge_counter_ = 0;
	}
	
	bool TOPPASMergerVertex::isSubtreeFinished()
	{
		if (!round_based_mode_)
		{
			/*	HACK: the (waiting) merger has to wait for all round-based mergers upstream
					to finish, which can only happen if it claims to have finished here	*/
			return true;
		}
		else
		{
			return TOPPASVertex::isSubtreeFinished();
		}
	}
}
