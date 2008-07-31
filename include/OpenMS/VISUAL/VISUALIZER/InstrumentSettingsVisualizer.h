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
//  License as published by the Free InstrumentSettings Foundation; either
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_INSTRUMENTSETTINGSVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_INSTRUMENTSETTINGSVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/InstrumentSettings.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS {
/**
@brief Class that displays all meta information for InstrumentSettings objects

This class provides all functionality to view the meta information of an object of type InstrumentSettings.
*/
	
	class InstrumentSettingsVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		InstrumentSettingsVisualizer(bool editable= FALSE, QWidget *parent =0);
		
		/// Loads the meta data from the object to the viewer.
		void load(InstrumentSettings &s);
	  
	private slots:
		 /// Saves the changes made to the meta data into the object.
		void store_();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject_();

	private:  
		/// Pointer to current object to keep track of the actual object
		InstrumentSettings *ptr_;
		/// Copy of current object for restoring the original values
		InstrumentSettings  tempinstrumentsettings_;
		/// Fills the comboboxes with current values
		void update_();
	  
		/** @name Edit fields
   */
    //@{
   	QComboBox *instrumentsettings_scan_mode_;
		QComboBox *instrumentsettings_polarity_;
		QLineEdit *instrumentsettings_mz_range_start_;
		QLineEdit *instrumentsettings_mz_range_stop_;
		//@}
		
					
	};
}
#endif
