// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2005 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free instrument; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Instrument Foundation; either
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
 
#ifndef OPENMS_VISUAL_VISUALIZER_INSTRUMENTVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_INSTRUMENTVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS {
/**
@brief Class that displays all meta information for HPLC objects

This class provides all functionality to view the meta information of an object of type HPLC.
*/
	
	class InstrumentVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		InstrumentVisualizer(bool editable= FALSE, QWidget *parent =0);
		/// Loads the meta data from the object to the viewer.
		void load(Instrument &s);
	  
	private slots:
		/// Saves the changes made to the meta data into the object.
		void store();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject();

	private:  
		/// Pointer to current object to keep track of the actual object
		Instrument *ptr_;
		/// Copy of current object for restoring the original values
		Instrument  tempinstrument_;
	  
		/** @name Edit fields and buttons
    */
    //@{
		QLineEdit *instrument_name_;
		QLineEdit *instrument_vendor_;
		QLineEdit *instrument_model_;
		QTextEdit *instrument_customizations_;
    //@}

     /** @name Some buttons.
		*/
    //@{
		QPushButton *savebutton_;
		QPushButton *cancelbutton_;
		//@}
		
		
					
	};
}
#endif
