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

#ifndef OPENMS_VISUAL_TOPPASOUTPUTFILELISTVERTEX_H
#define OPENMS_VISUAL_TOPPASOUTPUTFILELISTVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
	/**
		@brief A vertex representing an output file list
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASOutputFileListVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			
			/// Default constructor
			TOPPASOutputFileListVertex();
			/// Constructor
			TOPPASOutputFileListVertex(const QStringList& files);
			/// Copy constructor
			TOPPASOutputFileListVertex(const TOPPASOutputFileListVertex& rhs);
			/// Destructor
			virtual ~TOPPASOutputFileListVertex();
			/// Assignment operator
			TOPPASOutputFileListVertex& operator= (const TOPPASOutputFileListVertex& rhs);
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			/// Starts the workflow ending in this node
			void startComputation();
			/// Returns the list of files
			const QStringList& getFilenames();
			/// Called when the parent node has finished execution
			void finished();
			/// Checks if the given list of file names is valid
			bool fileNamesValid(const QStringList& files);
			/// Returns whether this vertex is ready
			bool isReady();
			/// Updates the status (if ready or not)
			void updateStatus();
			
		public slots:
		
			//documented in base class
			virtual void inEdgeHasChanged();
		
		signals:
			
			/// Emitted when an output file was written
			void outputFileWritten(const String& file);
		
		protected:
		
			/// The file names
			QStringList files_;
			/// Stores whether this vertex is ready (ready iff for every element of the list a name has been specified)
			bool ready_;
		
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
	};
}

#endif
