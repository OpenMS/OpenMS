// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

 
#ifndef OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_GRADIENTVISUALIZER_H

//OpenMS
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>
#include <OpenMS/METADATA/Gradient.h>

//STL
#include <vector>

class QIntValidator;

namespace OpenMS
{
	/**
		@brief GradientVisualizer is a visualizer class for objects of type gradient.
		
		Each HPLC objects contains a gradient object. A gradient objects contains a list of eluents, timepoints and percentage values. Values can be added to the list, or the whole list can be deleted.
	*/
	class OPENMS_GUI_DLLAPI GradientVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<Gradient>
	{
		Q_OBJECT

		public:
			
		  ///Constructor
			GradientVisualizer(bool editable = false, QWidget* parent = 0);
			
			//Docu in base class
			void load(Gradient& g);
		
		public slots:
			
		  //Docu in base class
			void store();
		
		protected slots:
			
			/// Add new timepoint to the list
			void addTimepoint_();
			/// Add new eluent to the list
			void addEluent_();
			///Delete all data from gradient
			void deleteData_();
			///Undo the changes made in the GUI.
			void undo_();
		
		protected: 
		  /// Loads a list of eluent, timepoint and percentage triplets.
			void loadData_();	
			/// Remove all data from layout
			void removeData_();
			
						
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
			std::vector< Int > timepoints_;
			//@}
			
			/** @name Some buttons.
			*/
	    //@{
			QPushButton* add_eluent_button_;
			QPushButton* add_timepoint_button_;
			QPushButton* removebutton_;
			//@}
		
			/// Array of temporary pointers to gradient edit fields
			std::vector< QLineEdit* > gradientdata_;
				
			/// Array of temporary pointers to gradient labels
			std::vector< QLabel* > gradientlabel_;
			
			/// Pointer to fields with actual data
			QLineEdit* percentage_;
			
			/// A validator to check the input for the new timepoint.
			QIntValidator* timepoint_vali_;		
			
			/// Counter to keep track of the actual row in the layout.
			int nextrow_;
				
			/// The layout to display the eluents, timepoints and percentages.
			QGridLayout* viewlayout_;		
			
			//Docu in base class
			void update_();
	};


}
#endif
