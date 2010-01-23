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
 
#ifndef OPENMS_VISUAL_VISUALIZER_PEPTIDEIDENTIFICATIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_PEPTIDEIDENTIFICATIONVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>


namespace OpenMS
{
	class MetaDataBrowser;

	/**
		@brief Class that displays all meta information for PeptideIdentification objects

		This class provides all functionality to view the meta information of an object of type PeptideIdentification.		
	*/
	class OPENMS_DLLAPI PeptideIdentificationVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<PeptideIdentification>
	{
		Q_OBJECT

		public:
		  ///Constructor
			PeptideIdentificationVisualizer(bool editable= FALSE, QWidget* parent =0, MetaDataBrowser* caller=0);
			
			/// Loads the meta data from the object to the viewer. Gets the id of the item in the tree as parameter.
			void load(PeptideIdentification& s, int tree_item_id);
			
		public slots:
			
			//Docu in base class
			void store();

		protected slots:
			
			///Undo the changes made in the GUI.
			void undo_();
			
			/** 
				@brief Updates the tree by calling MetaDataBrowser::updatePeptideHits(PeptideIdentification, int)
					
				Calls MetaDataBrowser::updatePeptideHits(PeptideIdentification, int).<br>
				Updates the tree depending of the protein significance threshold.<br>
				Only ProteinHits with a score superior or equal to the current threshold will be displayed.
			*/
			void updateTree_();
			
		protected:  	
		  /// Pointer to MetaDataBrowser
			MetaDataBrowser* pidv_caller_;
			/// The id of the item in the tree
			int tree_id_;
			
			///@name Edit fields and buttons
	    //@{
			QLineEdit* identifier_;
			QLineEdit* score_type_;
			QComboBox* higher_better_;
			QLineEdit* identification_threshold_;
			//@}
	
	    /// Threshold for foltering by score
			QLineEdit* filter_threshold_;
	};
}
#endif
