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

#ifndef OPENMS_VISUAL_TOPPASOUTPUTFILEVERTEX_H
#define OPENMS_VISUAL_TOPPASOUTPUTFILEVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
	/**
		@brief A vertex representing an output file
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASOutputFileVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			/// Default constructor
			TOPPASOutputFileVertex();
			/// Constructor
			TOPPASOutputFileVertex(const QString& file);
			/// Copy constructor
			TOPPASOutputFileVertex(const TOPPASOutputFileVertex& rhs);
			/// Destructor
			virtual ~TOPPASOutputFileVertex();
			/// Assignment operator
			TOPPASOutputFileVertex& operator= (const TOPPASOutputFileVertex& rhs);
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			/// Starts the workflow ending in this node
			void startComputation();
			/// Returns the file name
			const QString& getFilename();
			/// Called when the parent node has finished execution
			void finished();
			/// Checks if the given file name is valid
			bool fileNameValid(const QString& file);
			
		public slots:
		
			//documented in base class
			virtual void inEdgeHasChanged();
		
		signals:
			
			void outputFileWritten(const String& file);
		
		protected:
			
			/// The file name
			QString file_;
		
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
	};
}

#endif
