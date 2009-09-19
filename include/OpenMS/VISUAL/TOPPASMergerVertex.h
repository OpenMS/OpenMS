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
			/// Returns the list of output files
			QStringList getOutputList();
			/// Starts the pipeline execution recursively	
			void runRecursively();
			/// Forwards the pipeline execution downstream
			void forwardPipelineExecution();
			/// Determines whether all inputs are ready
			bool allInputsReady();
			/// Updates all temporary output file names
			void updateOutputFileNames();
			/// Sets whether the currently running pipeline has already been started at this vertex
      void setStartedHere(bool b);
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			
		protected:

			// Stores whether the currently running pipeline has already been started at this vertex
			bool started_here_;

			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
			
	};
}

#endif
