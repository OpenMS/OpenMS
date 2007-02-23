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

 
#ifndef OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/METADATA/Gradient.h>

//STL
#include <vector>

class QLabel;
class QLineEdit;
class QIntValidator;

namespace OpenMS {
/**
@brief GradientVisualizer is a visualizer class for objects of type gradient.

Each HPLC objects contains a gradient object. A gradient objects contains a list of eluents, timepoints and percentage values. Values can be added to the list, or the whole list can be deleted.
*/
	class GradientVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		GradientVisualizer(bool editable= FALSE, QWidget *parent =0);
		/// Loads the meta data from the object to the viewer.
		void load(Gradient &g);
		
		
	private slots:
	 
		/// Saves the information to Gradient Object.
		void store();
		/// Deletes all changes made in the viewer and restores the original data.
		void reject();
		/// Add new timepoint to the list
		void addTimepoint();
		/// Add new eluent to the list
		void addEluent();
		///Delete all data from gradient
		void deleteData();
		
	private: 
	  /// Loads a list of eluent, timepoint and percentage triplets.
		void loadData_();	
		/// Remove all data from layout
		void removeData_();
		/// Updates GUI with new data
		void update_();
		
		/// A validator to check the input for the new timepoint.
		QIntValidator *timepoint_vali_;	
		
					
	  /** @name Edit fields for new eluent-timepoint-percentage-triplets.
		*/
    //@{
		QLineEdit* new_eluent_;
		QLineEdit* new_timepoint_;
		//@}
		
		/** @name Arrays of string values containing eluent, timepoint and percentage values.
		*/
		//@{
		std::vector< String > eluents_;
		std::vector< SignedInt > timepoints_;
		//@}
		
	
		/// Array of temporary pointers to gradient edit fields
		std::vector< QLineEdit* > gradientdata_;
			
		/// Array of temporary pointers to gradient labels
		std::vector< QLabel* > gradientlabel_;
		
		/// Pointer to fields with actual data
		QLineEdit* percentage_;
		
		
		/** @name Some buttons.
		*/
    //@{
		QPushButton* savebutton_;
		QPushButton* cancelbutton_;
		QPushButton* add_eluent_button_;
		QPushButton* add_timepoint_button_;
		QPushButton* removebutton_;
		//@}
		
		/// Counter to keep track of the actual row in the layout.
		int nextrow_;
		
			
		/// The layout to display the eluents, timepoints and percentages.
		QGridLayout* viewlayout_;		
				
		/// Pointer to current object.
		Gradient* ptr_;
		/// Working-Copy of current object. 
		Gradient tempgradient_;
		
	};


}
#endif
