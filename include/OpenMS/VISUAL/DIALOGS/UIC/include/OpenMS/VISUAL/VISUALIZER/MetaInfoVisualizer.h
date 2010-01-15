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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

//STL
#include <vector>
#include <utility>

class QAbstractButton;
class QButtonGroup;

namespace OpenMS
{
	/**
		@brief MetaInfoVisualizer is a visualizer class for all classes that use one MetaInfo object as member.
		
		Meta information is an array of Type-Name-Value tupels. Classes that have a MetaInfo objects as a member can use this class to edit the MetaInfo object.
	*/
	class OPENMS_DLLAPI MetaInfoVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<MetaInfoInterface>
	{
		Q_OBJECT

		public:
		  ///Constructor
			MetaInfoVisualizer(bool editable = false, QWidget* parent = 0);
			
			//Docu in base class
			void load(MetaInfoInterface& m);

		public slots:
			
			//Docu in base class
			void store();

		protected slots:
			
		  /// Adds a new Type-Value pair to the MetaInfo Object.
			void add_();
			/// Clears out all fields.
			void clear_();
			/// Removes a selected Type-Value pair from the MetaInfo Object.
			void remove_(int);
			///Undo the changes made in the GUI.
			void undo_();
			
		protected: 
		  /// Loads all Type-Value pairs one after another. 
			void loadData_(UInt index);	
				
		  /** @name Edit fields for new Type-Value pair.
			*/
	    //@{
			QLineEdit* newkey_;
			QLineEdit* newvalue_;
			QLineEdit* newdescription_;
			//@}
			
			///@name Arrays of pointers to objects for temporary metaInfo data
			//@{
			std::vector< std::pair<UInt,QLineEdit*> > metainfoptr_;
			std::vector< std::pair<UInt,QLabel*> > metalabels_;
			std::vector< std::pair<UInt,QAbstractButton*> > metabuttons_;
			//@}		

			///@name Edit fields and buttons
	    //@{
			QPushButton* addbutton_;
			QPushButton* clearbutton_;
			QButtonGroup* buttongroup_;
			//@}
			
			/// Counter to keep track of the actual row in the layout.
			int nextrow_;
			
			/// The layout to display the Type-Value pairs.
			QGridLayout* viewlayout_;		
			
			/// Container for metainfo data.
			std::vector<UInt> keys_;
	};


}
#endif
