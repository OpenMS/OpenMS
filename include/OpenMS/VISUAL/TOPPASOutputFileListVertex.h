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
			/// Called when the parent node has finished execution
			void finished();
			/// Returns whether we are finished
			bool isFinished();
			/// Returns the directory where the output files are stored
			String getOutputDir();
			/// Creates all necessary directories (called by the scene before the pipeline is run)
			void createDirs(const QString& out_dir);
			/// Sets the topological sort number and removes invalidated tmp files
			virtual void setTopoNr(UInt nr);
			
		public slots:
		
			//documented in base class
			virtual void inEdgeHasChanged();
		
		signals:
			
			/// Emitted when an output file was written
			void outputFileWritten(const String& file);
			/// Emitted when the pipeline ending in this vertex is finished
			void iAmDone();
			
		protected:
		
			/// Stores whether the pipeline ending in this vertex has finished already
			bool finished_;
			/// The output file names
			QStringList files_;
		
			///@name reimplemented Qt events
      //@{
      void contextMenuEvent(QGraphicsSceneContextMenuEvent* event);
			//@}
	};
}

#endif
