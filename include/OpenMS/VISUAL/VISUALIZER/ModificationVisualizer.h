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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_VISUAL_VISUALIZER_MODIFICATIONVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_MODIFICATIONVISUALIZER_H

#include <OpenMS/METADATA/Modification.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLabel;
class QLineEdit;
class QComboBox;

namespace OpenMS {
/**
@brief Class that displays all meta information of modification objects.

This class provides all functionality to view the meta information of an object of type Modification.
*/
	class ModificationVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		ModificationVisualizer(bool editable= FALSE, QWidget *parent =0);
		
		/// Loads the meta data from the object to the viewer.
		void load(Modification &m);

	private slots:
	  /// Saves the changes made to the meta data into the object.
		void store_();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject_();

	private:  
		
				
		/// Sets the comboboxes with current values
		void updateMod_();
		
		/** @name Edit fields 
   */
    //@{
		QLineEdit *treatmenttype_;
		QTextEdit *treatmentcomment_;
		QLineEdit *modificationname_;
		QLineEdit *modificationmass_;
		QComboBox *modificationspecificity_;
		QLineEdit *modificationAA_;
		//@}
		
		/// Pointer to current object	 to keep track of the actual object.
		Modification *ptr_;
		/// Copy of current object for restoring the original values.
		Modification tempmod_;
		
		
	};


}
#endif
