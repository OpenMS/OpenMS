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
 
#ifndef OPENMS_VISUAL_VISUALIZER_SAMPLEVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_SAMPLEVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/Sample.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>

class QLineEdit;
class QComboBox;

namespace OpenMS {
/**
@brief Class that displays all meta information of sample objects.

This class provides all functionality to view the meta information of an object of type Sample.
*/

	
	class SampleVisualizer : public BaseVisualizer
	{
		Q_OBJECT

	public: 
	  /// Default constructor
		SampleVisualizer(bool editable= FALSE, QWidget *parent =0);
		
			/// Loads the meta data from the object to the viewer.
		void load(Sample &s);
	  
	private slots:
	  /// Saves the changes made to the meta data into the object.
		void store_();
		/// Deletes all changes made in the viewer and restores the original meta data.
		void reject_();

	private:  
		/// Pointer to current object	 to keep track of the actual object.
		Sample *ptr_;
		/// Copy of current object for restoring the original values.
		Sample  tempsample_;
	  	
		/// Sets the comboboxes with current values
		void update_();
		
		/** @name Edit fields
   */
    //@{
		QLineEdit *samplename_;
		QLineEdit *samplenumber_;
		QLineEdit *sampleorganism_;
		QTextEdit *samplecomment_;
		QComboBox *samplestate_;
		
		QLineEdit  *samplemass_;
		QLineEdit  *samplevolume_;
		QLineEdit  *sampleconcentration_;
		//@}
		
		
	};


}
#endif
