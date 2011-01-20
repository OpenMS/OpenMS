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
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPASINPUTFILELISTVERTEX_H
#define OPENMS_VISUAL_TOPPASINPUTFILELISTVERTEX_H

#include <OpenMS/VISUAL/TOPPASVertex.h>

namespace OpenMS
{
	/**
		@brief A vertex representing an input file list
	
		@ingroup TOPPAS_elements
	*/
	class OPENMS_DLLAPI TOPPASInputFileListVertex
		: public TOPPASVertex
	{
		Q_OBJECT
		
		public:
			
			/// Default constructor
			TOPPASInputFileListVertex();
			/// Constructor
			TOPPASInputFileListVertex(const QStringList& files);
			/// Copy constructor
			TOPPASInputFileListVertex(const TOPPASInputFileListVertex& rhs);
			/// Destructor
			virtual ~TOPPASInputFileListVertex();
			/// Assignment operator
			TOPPASInputFileListVertex& operator= (const TOPPASInputFileListVertex& rhs);
			/// Returns the list of files
      const QStringList& getInputFilenames();
			/// Sets the list of files
			void setFilenames(const QStringList& files);
			/// Starts all tools below this node
			void startPipeline();
			// documented in base class
			virtual void paint(QPainter* painter, const QStyleOptionGraphicsItem* option, QWidget* widget);
			// documented in base class
			virtual QRectF boundingRect() const;
			// documented in base class
			virtual QPainterPath shape () const;
			// documented in base class
			virtual void checkListLengths(QStringList& unequal_per_round, QStringList& unequal_over_entire_run);
			/// Checks if the given list of file names is valid
			bool fileNamesValid(const QStringList& files);
			/// Shows the dialog for editing the files
			void showFilesDialog();
      /// Opens the folders of the input files
      void openContainingFolder();
      /// Opens the files in TOPPView
			void openInTOPPView();
			/// Returns the key (for applying resources from a resource file)
			const QString& getKey();
			/// Sets the key (for applying resources from a resource file)
			void setKey(const QString& key);
			
		protected:
		
			/// The file names
			QStringList files_;
			/// The key of this input node (for applying resources from a resource file)
			QString key_;
		
			///@name reimplemented Qt events
      //@{
      void mouseDoubleClickEvent(QGraphicsSceneMouseEvent* e);
			//@}
			
	};
}

#endif
