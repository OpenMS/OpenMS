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
			merge_counter_(0),
			notified_parents_counter_(0)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASMergerVertex::TOPPASMergerVertex(const TOPPASMergerVertex& rhs)
		:	TOPPASVertex(rhs),
			round_based_mode_(rhs.round_based_mode_),
			merge_counter_(rhs.merge_counter_),
			notified_parents_counter_(rhs.notified_parents_counter_)
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
		merge_counter_ = rhs.merge_counter_;
		notified_parents_counter_ = rhs.notified_parents_counter_;
		
		return *this;
	}
	
	void TOPPASMergerVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
	}
	
	QStringList TOPPASMergerVertex::getCurrentOutputList()
	{
		__DEBUG_BEGIN_METHOD__
		
		int tmp_merge_counter = merge_counter_;
		if (mergeComplete())
		{
			/*	Only relevant if nested mergers have incompatible input list lengths.
					If this merger has finished merging already (and so do all upstream
					mergers, if any) and this merger is still asked to return its output files
					(by another merger further downstream), then we just return the list
					from the last merging round, although the result is probably bogus. */
			--tmp_merge_counter;
		}
		
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
					out_files << file_names[tmp_merge_counter];
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
					out_files << source_list->getFilenames()[tmp_merge_counter];
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
					out_files << source_merger->getCurrentOutputList()[tmp_merge_counter];
				}
				else
				{
					out_files << source_merger->getAllCollectedFiles_();
				}
				continue;
			}
		}
		
		__DEBUG_END_METHOD__
		
		return out_files;
	}
	
	QStringList TOPPASMergerVertex::getAllCollectedFiles_()
	{
		__DEBUG_BEGIN_METHOD__
		
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
		
		__DEBUG_END_METHOD__
		return out_files;
	}

	bool TOPPASMergerVertex::allInputsReady()
	{
		__DEBUG_BEGIN_METHOD__
		
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
		  TOPPASToolVertex* tv = qobject_cast<TOPPASToolVertex*>((*it)->getSourceVertex());
		  if (tv && !tv->isFinished())
		  {
		    // some tool that we depend on has not finished execution yet --> do not start yet
				__DEBUG_END_METHOD__
		    return false;
		  }
			TOPPASMergerVertex* mv = qobject_cast<TOPPASMergerVertex*>((*it)->getSourceVertex());
			if (mv && !mv->allInputsReady())
			{
				// some input of a merger that we depend on not finished yet
				__DEBUG_END_METHOD__
				return false;
			}
		}
		__DEBUG_END_METHOD__
		return true;
	}
	
	void TOPPASMergerVertex::forwardPipelineExecution(bool start_merge_all)
	{
		__DEBUG_BEGIN_METHOD__
		
		if (!allInputsReady())
		{
			debugOut_("Not run because !allInputsReady()");
			
			__DEBUG_END_METHOD__
			return;
		}
		
		if (!round_based_mode_ && !start_merge_all)
		{
			debugOut_("Not run because !round_based_mode_ && !start_merge_all");
			
			__DEBUG_END_METHOD__
			return;
		}
		
		if (subtree_finished_ && notified_parents_counter_ != in_edges_.size()) // merger was run and has finished, but not all parents have been informed yet
		{
			debugOut_("Not run because !notified_parents_counter != in_edges_.size()");
			
			__DEBUG_END_METHOD__
			return;
		}
		
		debugOut_("About to run");
		
		if (subtree_finished_)
		{
			// this merger was run and has finished already and is now asked to run a second time --> reset subtree first
			debugOut_("Resetting subtree");
			
			TOPPASVertex::resetSubtreeUpToNextMergers(true);
		}
		
		debugOut_("Proceeding in pipeline...");
		
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
		
		debugOut_(String("All children nodes run! Incremented merge_counter_ to ")+merge_counter_);
		
		__DEBUG_END_METHOD__
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
		__DEBUG_BEGIN_METHOD__
		
		bool complete = merge_counter_ == numIterations_();
		debugOut_(String("Returning ")+(complete ? "true" : "false"));
		
		__DEBUG_END_METHOD__
		return complete;
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
		__DEBUG_BEGIN_METHOD__
		
		int min_length = std::numeric_limits<int>::max();
		
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* source = (*it)->getSourceVertex();
			
			debugOut_(String("Checking parent ")+source->getTopoNr());
			
			TOPPASToolVertex* source_tool = qobject_cast<TOPPASToolVertex*>(source);
			if (source_tool)
			{
				const QVector<QStringList>& output_files = source_tool->getCurrentOutputFileNames();
				if (output_files.empty())
				{
					// that node was reset
					debugOut_("Parent file list empty, returning -1");
					
					__DEBUG_END_METHOD__
					return -1;
				}
				
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
				// no check needed here, input files are always ready	
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
		
		debugOut_(String("Returning ")+min_length);
		
		__DEBUG_END_METHOD__
		
		return min_length;
	}
	
	int TOPPASMergerVertex::numIterations_()
	{
		__DEBUG_BEGIN_METHOD__
		
		int iterations = round_based_mode_ ? minInputListLength_() : 1;
		debugOut_(String("Returning ")+iterations);
		
		__DEBUG_END_METHOD__
		return iterations;
	}
	
	void TOPPASMergerVertex::checkIfSubtreeFinished()
	{
		__DEBUG_BEGIN_METHOD__
		
		// parents still ready? (might not be the case if they have been reset --> wait)
		if (!allInputsReady())
		{
			return;
		}

		// all children's subtrees finished?
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* target = (*it)->getTargetVertex();
			
			debugOut_(String("Checking if subtree of child ")+target->getTopoNr()+" finished...");
			
			if (!target->isSubtreeFinished())
			{
				debugOut_("... subtree NOT finished.");
				__DEBUG_END_METHOD__
				return;
			}
			else
			{
				debugOut_("... subtree finished.");
			}
		}
		
		/*	Proceed with next merging round, if this merger is still running and
				if all children have finished already. If merge is complete, propagate
				this upwards. */
		
		if (!mergeComplete())
		{
			debugOut_("Merge not complete yet, resetting subtree...");
			// proceed with next merge iteration
			TOPPASVertex::resetSubtreeUpToNextMergers(false);
			debugOut_("Subtree successfully reset");
			forwardPipelineExecution();
			debugOut_("forwardPipelineExecution() successfully called");
		}
		else
		{
			debugOut_("Merge complete - notifying parents");
			// merge complete --> propagate this upwards (for other round-based mergers)
			subtree_finished_ = true;
			for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
			{
				TOPPASVertex* source = (*it)->getSourceVertex();
				debugOut_(String("Notifying parent ")+source->getTopoNr()+" that merging round complete");
				++notified_parents_counter_;
				debugOut_(String("Increased notified_parents_counter_ to ")+notified_parents_counter_);
				source->checkIfSubtreeFinished();
			}
			debugOut_("Upstream notification successful");
		}
		
		__DEBUG_END_METHOD__
	}
	
	void TOPPASMergerVertex::reset(bool /*reset_all_files*/)
	{
		__DEBUG_BEGIN_METHOD__
		
		TOPPASVertex::reset(false);
		merge_counter_ = 0;
		notified_parents_counter_ = 0;
		
		__DEBUG_END_METHOD__
	}
	
	bool TOPPASMergerVertex::isSubtreeFinished()
	{
		if (!round_based_mode_)
		{
			/*	HACK: the waiting "merge all" merger has to wait for all round-based mergers upstream
					to finish, which can only happen if it claims to have finished here	*/
			return true;
		}
		else
		{
			return TOPPASVertex::isSubtreeFinished();
		}
	}
	
	void TOPPASMergerVertex::resetSubtreeUpToNextMergers(bool /*including_this_node*/)
	{
		__DEBUG_BEGIN_METHOD__
		debugOut_("Not doing anything");
		// this node _is_ a merger --> stop here, do not reset
		__DEBUG_END_METHOD__
	}
}
