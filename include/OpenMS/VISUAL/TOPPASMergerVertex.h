// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#ifndef OPENMS_VISUAL_TOPPASMERGERVERTEX_H
#define OPENMS_VISUAL_TOPPASMERGERVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
	/**
		@brief A special vertex that allows to merge several inputs.
		
		A special vertex that allows to merge several inputs. Mergers have two modes: The normal,
		round-based merging mode and a "wait & merge all" mode. In round-based mode, a merger
		first takes the first files of each incoming file list and merges them into a list (which
		has as many elements as the merger has incoming edges). All tools this merger has outgoing
		edges to are called with this merged list as input files.
		As soon as they have all been processed, the second files of each incoming list are merged
		and the tools below the merger are called again, and so on.
		If several mergers are nested, mergers further upstream wait until mergers further downstream
		have performed all their merging rounds and then perform their own next merging round.
		
		In "wait & merge all" mode, the merger first waits for all upstream mergers to finish all
		their merging rounds and then merges all collected files from all merging rounds for all
		incoming edges into one single list and calls the next tool with this list of files as input.
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASMergerVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			/// Default constructor
			TOPPASMergerVertex();
			/// Copy constructor
			TOPPASMergerVertex(const TOPPASMergerVertex& rhs);
			/// Destructor
			virtual ~TOPPASMergerVertex();
			/// Assignment operator
			TOPPASMergerVertex& operator= (const TOPPASMergerVertex& rhs);
			/// Returns the current list of output files
			QStringList getCurrentOutputList();
			/// Forwards the pipeline execution downstream
			void forwardPipelineExecution();
			/// Determines whether all inputs are ready
			bool allInputsReady();
			/// Determines whether all merge rounds have been performed
			bool mergeComplete();
			/// Determines whether this merger is merging round based or merging all inputs into one list
			bool roundBasedMode();
			/// Sets whether this merger is merging round based or merging all inputs into one list
			void setRoundBasedMode(bool b);
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			// documented in base class
			virtual void checkIfSubtreeFinished();
			// documented in base class
			virtual bool areAllUpstreamMergersFinished();
			// documented in base class
			virtual void reset(bool reset_all_files = false);
			// documented in base class
			virtual bool isSubtreeFinished();
			// documented in base class
			virtual void checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run);
			// documented in base class
			virtual void markUnreachable();
			
		protected:

			/// Stores whether this merger is merging round based or merging all inputs into one list
			bool round_based_mode_;
			/// The counter for the merging process
			int merge_counter_;
			/// Is set to true while the parents are notified that this node has finished merging
			bool currently_notifying_parents_;
			/// The maximum length of all incoming lists
			int max_input_list_length_;
			/// Stores the last list of output files that was processed
			QStringList last_output_files_;
			// documented in base class
			using TOPPASVertex::reachable_;

			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
			/// Returns the number of iterations we have to perform
			int numIterations_();
			/// Returns the list of all written output files (possibly over several merging rounds) of the parents
			QStringList getAllCollectedFiles_();
			
	};
}

#endif
