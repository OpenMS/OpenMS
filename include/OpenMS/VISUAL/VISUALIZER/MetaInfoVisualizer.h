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
 
#ifndef OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_METAINFOVISUALIZER_H

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/METADATA/MetaInfoInterface.h>

//STL
#include <vector>
#include <utility>

class QLabel;
class QLineEdit;
class QAbstractButton;
class QButtonGroup;

namespace OpenMS {
/**
@brief MetaInfoVisualizer is a visualizer class for all classes that use one MetaInfo object as member.

Meta information is an array of Type-Name-Value tupels. Classes that have a MetaInfo objects as a member can use this class to edit the MetaInfo object.
*/
	class MetaInfoVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		MetaInfoVisualizer(bool editable= FALSE, QWidget *parent =0);
		/// Loads the meta data from the object to the viewer.
		void load(MetaInfoInterface &m);
		
		
	private slots:
	  /// Adds a new Type-Value pair to the MetaInfo Object.
		void add();
		/// Clears out all fields.
		void clear();
		/// Removes a selected Type-Value pair from the MetaInfo Object.
		void remove(int);
		/// Saves the information to MetaInfo Object.
		void store();
		/// Deletes all changes made in the viewer and restores the original data.
		void reject();

	private: 
	  /// Loads all Type-Value pairs one after another. 
		void loadData_(UnsignedInt index);	
			
	  /** @name Edit fields for new Type-Value pair.
		*/
    //@{
		QLineEdit* newkey_;
		QLineEdit* newvalue_;
		QLineEdit* newdescription_;
		//@}
		
		/** @name Arrays of pointers to objects for temporary metaInfo data
		*/
		//@{
		std::vector< std::pair<UnsignedInt,QLineEdit*> > metainfoptr_;
		std::vector< std::pair<UnsignedInt,QLabel*> > metalabels_;
		std::vector< std::pair<UnsignedInt,QAbstractButton*> > metabuttons_;
		//@}		
		
		/** @name Some buttons.
		*/
    //@{
		QPushButton* savebutton_;
		QPushButton* cancelbutton_;
		QPushButton* addbutton_;
		QPushButton* clearbutton_;
		QButtonGroup* buttongroup_;
		//@}
		
		/// Counter to keep track of the actual row in the layout.
		int nextrow_;
		
		/// The layout to display the Type-Value pairs.
		QGridLayout* viewlayout_;		
		
		/// Pointer to current object.
		MetaInfoInterface* ptr_;
		/// Working-Copy of current object. 
		MetaInfoInterface tempmeta_;
		
		/// Container for metainfo data.
		std::vector<UnsignedInt> keys_;
		
	};


}
#endif
