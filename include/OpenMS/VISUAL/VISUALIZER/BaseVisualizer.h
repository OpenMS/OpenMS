// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: stefan_heess $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_VISUALIZER_BASEVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_BASEVISUALIZER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/VISUAL/DataTable.h>

namespace OpenMS {
/**
@brief A base class for all visualizer classes

This class is a basic class for all visualizer classes. 
It provides some functions that are re-implemented in the subclasses. <br>
Increases ease of use to store data.
*/	
	class BaseVisualizer : public DataTable
	{
		Q_OBJECT

		public:
			/// Default constructor 
			BaseVisualizer(bool editable=FALSE, QWidget *parent =0);
			/// Returns the type of the visualizer class.
			String getType();
			///Defines a friend class that can use the functionality of the subclasses.
			friend class MSMetaDataExplorer;
			
		protected:
			/// Adds buttons common to all visualizers
			void finishAdding_();
			
			///The type of the object to be displayed.
			String type_;
			
	    ///Undo buttons.	
			QPushButton *undobutton_;		
			
			private slots:
			///Saves the changes made to the object.
			virtual void store()=0;
			
			///Undo the changes made to the object.
			virtual void reject()=0;
		
			
				
	};


}
#endif
