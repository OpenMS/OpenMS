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
//  License as published by the Free Identification Foundation; either
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_IDENTIFICATIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_IDENTIFICATIONVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS
{
	class MSMetaDataExplorer;

	/**
		@brief Class that displays all meta information for Identification objects

		This class provides all functionality to view the meta information of an object of type Identification.
	*/
	
	class IdentificationVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		IdentificationVisualizer(bool editable= FALSE, QWidget *parent =0, MSMetaDataExplorer *caller=0);
		/// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
		void load(Identification &s, int tree_item_id);
	  
	private slots:
		/// Save all changes
		void store();
		/// Restore all changes
		void reject();
		/** 
		@brief Updates the tree by calling MSMetaDataExplorer::updatePeptideHits(Identification, int)
			
		Calls MSMetaDataExplorer::updatePeptideHits(Identification, int).<br>
		Updates the tree depending of the protein significance threshold.<br>
		Only ProteinHits with a score superior or equal to the current threshold will be displayed.
		*/
		void updateTree();
		
		/** 
		@brief Updates the tree by calling MSMetaDataExplorer::updatePeptideHits(Identification, int, String, String)
			
		Calls MSMetaDataExplorer::updatePeptideHits(Identification, int, String, String).<br>
		Updates the tree depending on the searched ProteinHit.<br>
		Only PeptideHits that reference the searched ProteinHit will be displayed.
		*/
		void searchRefPeptides();
		
		/** 
		@brief Updates the tree by calling MSMetaDataExplorer::updatePeptideHits(Identification, int, String, String)
			
		Calls MSMetaDataExplorer::updatePeptideHits(Identification, int, String, String).<br>
		Updates the tree depending on the existing ProteinHits.<br>
		Only PeptideHits that do not refererence any ProteinHit will be displayed.
		*/
		void searchNonRefPeptides();
		
	private:  
		/// Pointer to current object to keep track of the actual object
		Identification *ptr_;
		/// Copy of current object for restoring the original values
		Identification  tempidentification_;
	  /// Pointer to MSMetaDataExplorer
		MSMetaDataExplorer *pidv_caller_;
		/// The id of the item in the tree
		int tree_id_;
		
		/** @name Edit fields
   */
    //@{
		QLineEdit *identification_date_;
		QLineEdit *identification_ref_date_;
		QLineEdit *identification_threshold_;
		QLineEdit *identification_acc_;
		//@}

     /** @name Some buttons.
		*/
    //@{
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		QPushButton *updatebutton_;
		QPushButton *updatebutton2_;
		QPushButton *updatebutton3_;
		//@}
		
		
					
	};
}
#endif
