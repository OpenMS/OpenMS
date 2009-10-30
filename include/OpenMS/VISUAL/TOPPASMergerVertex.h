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

#ifndef OPENMS_VISUAL_TOPPASMERGERVERTEX_H
#define OPENMS_VISUAL_TOPPASMERGERVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
	/**
		@brief A special vertex that allows to merge several inputs into a single output file list
	
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
			void forwardPipelineExecution(bool start_merge_all = false);
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
			virtual void checkIfAllUpstreamMergersFinished();
			// documented in base class
			virtual void checkIfSubtreeFinished();
			// documented in base class
			virtual void reset(bool reset_all_files = false, bool mergers_finished = true);
			// documented in base class
			virtual bool isSubtreeFinished();
			// documented in base class
			virtual bool isAllUpstreamMergersFinished();
			
		protected:

			/// Stores whether this merger is merging round based or merging all inputs into one list
			bool round_based_mode_;
			/// The counter for the merging process
			int merge_counter_;

			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
			/// Returns the minimum length of all incoming lists
			int minInputListLength_();
			/// Returns the number of iterations we have to perform
			int numIterations_();
			/// Returns the list of all written output files (during the entire pipeline execution) of the parents
			QStringList getAllCollectedFiles_();

			
	};
}

#endif
