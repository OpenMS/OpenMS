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
//  License as published by the Free MassAnalyzer Foundation; either
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_MASSANALYZERVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_MASSANALYZERVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>


namespace OpenMS
{
	/**
		@brief Class that displays all meta information for MassAnalyzer objects
		
		This class provides all functionality to view the meta information of an object of type MassAnalyzer.
	*/
	class OPENMS_DLLAPI MassAnalyzerVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<MassAnalyzer>
	{
		Q_OBJECT
	
		public:
			
		  ///Constructor
			MassAnalyzerVisualizer(bool editable = false, QWidget* parent = 0);
			
		public slots:
			
		  //Docu in base class
			void store();
		
		protected slots:
			
			///Undo the changes made in the GUI.
			void undo_();
			
		protected:  
				
			///@name edit fields to modify properties
	    //@{
	    QLineEdit* order_;
			QLineEdit* res_;
			QLineEdit* acc_;
			QLineEdit* scan_rate_;
			QLineEdit* scan_time_;
			QLineEdit* TOF_;
			QLineEdit* iso_;
			QLineEdit* final_MS_;
			QLineEdit* magnetic_fs_;
			QComboBox* type_;
			QComboBox* res_method_;
			QComboBox* res_type_;
			QComboBox* scan_dir_;
			QComboBox* scan_law_;
			QComboBox* reflectron_state_;
			//@}
			
			//Docu in base class
			void update_();
	};
}
#endif
