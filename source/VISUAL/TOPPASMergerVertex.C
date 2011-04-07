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

#include <OpenMS/VISUAL/TOPPASMergerVertex.h>
#include <OpenMS/VISUAL/TOPPASInputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASOutputFileListVertex.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>
#include <OpenMS/VISUAL/TOPPASScene.h>

namespace OpenMS
{
	TOPPASMergerVertex::TOPPASMergerVertex()
		:	TOPPASVertex(),
			round_based_mode_(true)
	{
		pen_color_ = Qt::black;
		brush_color_ = Qt::lightGray;
	}
	
	TOPPASMergerVertex::TOPPASMergerVertex(const TOPPASMergerVertex& rhs)
		:	TOPPASVertex(rhs),
			round_based_mode_(rhs.round_based_mode_)
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

		return *this;
	}
	
	void TOPPASMergerVertex::mouseDoubleClickEvent(QGraphicsSceneMouseEvent* /*e*/)
	{
	  setRoundBasedMode(!round_based_mode_);
	}

	void TOPPASMergerVertex::paint(QPainter* painter, const QStyleOptionGraphicsItem* /*option*/, QWidget* /*widget*/)
	{
		__DEBUG_BEGIN_METHOD__
		
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
		
    if (round_total_ != -1) // draw round number
    {
      text = QString::number(round_counter_)+" / "+QString::number(round_total_);
  		text_boundings = painter->boundingRect(QRectF(0,0,0,0), Qt::AlignCenter, text);
	  	painter->drawText(-(int)(text_boundings.width()/2.0), 31, text);
    }

		//topo sort number
		qreal x_pos = -36.0;
		qreal y_pos = -23.0; 
		painter->drawText(x_pos, y_pos, QString::number(topo_nr_));
		
		
		// recycling status
    if (this->allow_output_recycling_)
    {
      painter->setPen(Qt::green);
      painter->drawChord(-7,-32, 14, 14, 0*16, 180*16);
    }

    __DEBUG_END_METHOD__
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
	
	bool TOPPASMergerVertex::roundBasedMode()
	{
		return round_based_mode_;
	}
	
	void TOPPASMergerVertex::setRoundBasedMode(bool b)
	{
		reset();
		round_based_mode_ = b;
		emit somethingHasChanged();
	}
	
	
	void TOPPASMergerVertex::markUnreachable()
	{
		//only mark as unreachable if all inputs are unreachable. otherwise the dead inputs will just be ignored.
		bool some_input_reachable_ = false;
		for (EdgeIterator it = inEdgesBegin(); it != inEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getSourceVertex();
			if (tv->isReachable())
			{
				some_input_reachable_ = true;
				break;
			}
		}
		if (!some_input_reachable_)
		{
			TOPPASVertex::markUnreachable();
		}
	}

  void TOPPASMergerVertex::run()
  {
    //check if everything ready
		if (!isUpstreamReady())	return;

    RoundPackages pkg;
    String error_msg("");
    bool success = buildRoundPackages(pkg, error_msg);
    if (!success)
    {
      std::cerr << "Could not retrieve input files from upstream nodes...\n";
      emit mergeFailed((String("Merger #") + this->getTopoNr() + " failed. " + error_msg).toQString());
      return;
    }

    /// update round status
    Size input_rounds = pkg.size();
    round_total_ = (round_based_mode_ ? (int) input_rounds : 1) ; // for round based: take number of rounds from previous tool(s) - should all be equal
    round_counter_ = 0;        // once round_counter_ reaches round_total_, we are done

    // clear output file list
		output_files_.clear();
    output_files_.resize(round_total_); // #rounds

    // Do the virtual merging (nothing more than reorganizing filenames)
    for (Size round = 0; round < input_rounds; ++round)
    {
      QStringList files;
      // warning: ite->first (i.e. target-in param could be -1,-2,... etc to cover all incoming edges (they all have -1 theoretically - see buildRoundPackages())
      for (RoundPackageConstIt ite = pkg[round].begin();
                               ite!= pkg[round].end();
                               ++ite)
			{
        files.append( ite->second.filenames ); // concat filenames from all incoming edges
      }
      Size round_index = (round_based_mode_ ? round : 0);
      output_files_[round_index][-1].filenames.append ( files ); // concat over all rounds (if required)
    }

    round_counter_ = round_total_;
    finished_ = true;
    
    // call all childs, proceed in pipeline
		for (EdgeIterator it = outEdgesBegin(); it != outEdgesEnd(); ++it)
		{
			TOPPASVertex* tv = (*it)->getTargetVertex();
			debugOut_(String("Starting child ") + tv->getTopoNr());
			tv->run();
		}

  }
}
