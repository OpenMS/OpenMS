// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_VISUAL_VISUALIZER_MASSANALYZERVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_MASSANALYZERVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/MassAnalyzer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS {
/**
@brief Class that displays all meta information for MassAnalyzer objects

This class provides all functionality to view the meta information of an object of type MassAnalyzer.
*/
	
	class MassAnalyzerVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		MassAnalyzerVisualizer(bool editable= FALSE, QWidget *parent =0);
		
		/// Loads the meta data from the object to the viewer.
		void load(MassAnalyzer &s);
	  
	private slots:
		/// Saves the changes made to the meta data into the object.
		void store_();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject_();
		
	private:  
		/// Pointer to current object to keep track of the actual object
		MassAnalyzer *ptr_;
		/// Copy of current object for restoring the original values
		MassAnalyzer  tempmassanalyzer_;
		/// Fills the comboboxes with current values
		void update_();
	  
		
			
		/** @name edit fields to modify properties
   */
    //@{
		QLineEdit *massanalyzer_res_;
		QLineEdit *massanalyzer_acc_;
		QLineEdit *massanalyzer_scan_rate_;
		QLineEdit *massanalyzer_scan_time_;
		QLineEdit *massanalyzer_TOF_;
		QLineEdit *massanalyzer_iso_;
		QLineEdit *massanalyzer_final_MS_;
		QLineEdit *massanalyzer_magnetic_fs_;
		//@}

		
		/** @name Comboboxes to choose properties
   */
    //@{
		QComboBox *massanalyzer_type_;
		QComboBox *massanalyzer_res_method_;
		QComboBox *massanalyzer_res_type_;
		QComboBox *massanalyzer_scan_func_;
		QComboBox *massanalyzer_scan_dir_;
		QComboBox *massanalyzer_scan_law_;
		QComboBox *massanalyzer_tandem_scan_method_;
		QComboBox *massanalyzer_reflectron_state_;
		//@}
		
		
					
	};
}
#endif
