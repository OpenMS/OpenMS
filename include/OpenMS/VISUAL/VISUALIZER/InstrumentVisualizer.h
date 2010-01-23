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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_VISUAL_VISUALIZER_INSTRUMENTVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_INSTRUMENTVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/Instrument.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>


namespace OpenMS
{
	/**
		@brief Class that displays all meta information for an MS instrument
		
		This class provides all functionality to view the meta information of an object of type Instrument.
	*/
	class OPENMS_DLLAPI InstrumentVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<Instrument>
	{
		Q_OBJECT

		public:
			
		  ///Constructor
			InstrumentVisualizer(bool editable = false, QWidget* parent = 0);
			
		public slots:
			
			//Docu in base class
			void store();
		
		protected slots:
			
			///Undo the changes made in the GUI.
			void undo_();
	
		protected:  

			///@name Edit fields and buttons
	    //@{
			QLineEdit* name_;
			QLineEdit* vendor_;
			QLineEdit* model_;
			QTextEdit* customizations_;
			QComboBox* ion_optics_;
	    //@}
	    
	    //Docu in base class
			void update_();
	};
}
#endif
