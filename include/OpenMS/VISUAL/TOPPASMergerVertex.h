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
		has as many elements as the merger has incoming edges).
		
		In "wait & merge all" mode, the merger first waits for all upstream mergers to finish all
		their merging rounds and then merges all collected files from all merging rounds for all
		incoming edges into one single list and calls the next tool with this list of files as input.
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_GUI_DLLAPI TOPPASMergerVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			/// Default constructor
			TOPPASMergerVertex();
      /// Constructor
      TOPPASMergerVertex(bool round_based);
			/// Copy constructor
			TOPPASMergerVertex(const TOPPASMergerVertex& rhs);
			/// Destructor
			virtual ~TOPPASMergerVertex();
			/// Assignment operator
			TOPPASMergerVertex& operator= (const TOPPASMergerVertex& rhs);
      /// check if upstream nodes are finished and call downstream nodes
      virtual void run();
			/// Determines whether all inputs are ready (only a problem in mergers, when called from upstream)
			bool allInputsReady();
			/// Determines whether this merger is merging round based or merging all inputs into one list
			bool roundBasedMode();
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			// documented in base class
			virtual void markUnreachable();
			
      public slots:

      signals:
      /// Emitted when merging upstream data failed
      void mergeFailed(const QString message);

		protected:

			/// Stores whether this merger is merging round based or merging all inputs into one list
			bool round_based_mode_;

      ///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
		
	};
}

#endif
