// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free proteinIdentification; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free ProteinIdentification Foundation; either
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_PROTEINIDENTIFICATIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_PROTEINIDENTIFICATIONVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS
{
	class MSMetaDataExplorer;

	/**
	@brief Class that displays all meta information for ProteinIdentification objects

	This class provides all functionality to view the meta information of an object of type ProteinIdentification.
	*/
	
	class ProteinIdentificationVisualizer : public BaseVisualizer
	{
		Q_OBJECT

public: 
	  // Default constructor
		ProteinIdentificationVisualizer(bool editable= FALSE, QWidget *parent =0, MSMetaDataExplorer *caller=0);
		/// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
		void load(ProteinIdentification &s, int tree_item_id);
	  
	private slots:
		/// Save all changes
		void store();
		/// Restore all changes
		void reject();
		/** 
		@brief Updates the tree by calling MSMetaDataExplorer::updateProteinHits()
			
		Calls MSMetaDataExplorer::updateProteinHits().<br>
		Updates the tree depending of the protein significance threshold.<br>
		Only ProteinHits with a score superior or equal to the current threshold will be displayed.
		*/
		void updateTree();
		
	private:  
		/// Pointer to current object to keep track of the actual object
		ProteinIdentification *ptr_;
		/// Copy of current object for restoring the original values
		ProteinIdentification  tempproteinidentification_;
	  /// Pointer to MSMetaDataExplorer
		MSMetaDataExplorer *pidv_caller_;
		/// The id of the item in the tree
		int tree_id_;
		
		/** @name Edit fields
   */
    //@{
		QLineEdit *proteinidentification_date_;
		QLineEdit *proteinidentification_threshold_;
		//@}

     /** @name Some buttons.
		*/
    //@{
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		QPushButton *updatebutton_;
		//@}
		
		
					
	};
}
#endif
