// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_PROTEINIDENTIFICATIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_PROTEINIDENTIFICATIONVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>


namespace OpenMS
{
	class MetaDataBrowser;

	/**
		@brief Class that displays all meta information for ProteinIdentification objects
	
		This class provides all functionality to view the meta information of an object of type Identification.	
	*/
	class OPENMS_DLLAPI ProteinIdentificationVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<ProteinIdentification>
	{
		Q_OBJECT
	
		public:
		  ///Constructor
			ProteinIdentificationVisualizer(bool editable= FALSE, QWidget* parent =0, MetaDataBrowser* caller=0);
		
			/// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
			void load(ProteinIdentification& s, int tree_item_id);
		  
		public slots:
			
			//Docu in base class
			void store();

		protected slots:
			
			/** 
				@brief Updates the tree by calling MetaDataBrowser::updateProteinHits()
					
				Calls MetaDataBrowser::updateProteinHits().<br>
				Updates the tree depending of the protein significance threshold.<br>
				Only ProteinHits with a score superior or equal to the current threshold will be displayed.
			*/
			void updateTree_();
			
			///Undo the changes made in the GUI.
			void undo_();
			
		protected:
			
		  /// Pointer to MetaDataBrowser
			MetaDataBrowser* pidv_caller_;
			/// The id of the item in the tree
			int tree_id_;
			
			///@name Edit fields and buttons
	    //@{
			QLineEdit* engine_;
			QLineEdit* engine_version_;
			QLineEdit* identification_date_;
			QLineEdit* identification_threshold_;
			QLineEdit* identifier_;
			QLineEdit* score_type_;
			QComboBox* higher_better_;
			
			QLineEdit* db_;
			QLineEdit* db_version_;
			QLineEdit* taxonomy_;
			QLineEdit* charges_;
			QLineEdit* missed_cleavages_;
			QLineEdit* peak_tolerance_;
			QLineEdit* precursor_tolerance_;
			QComboBox* mass_type_;
			QComboBox* enzyme_;
			//@}
	
	    /// Threshold for foltering by score
			QLineEdit* filter_threshold_;
	};
}
#endif //OPENMS_VISUAL_VISUALIZER_PROTEINIDENTIFICATIONVISUALIZER_H
