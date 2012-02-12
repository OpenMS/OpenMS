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
//  License as published by the Free SpectrumSettings Foundation; either
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------
 
#ifndef OPENMS_VISUAL_VISUALIZER_SPECTRUMSETTINGSVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_SPECTRUMSETTINGSVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

class QTextEdit;

namespace OpenMS 
{
	/**
		@brief Class that displays all meta information for SpectrumSettings objects
		
		This class provides all functionality to view the meta information of an object of type SpectrumSettings.
	*/
	class OPENMS_GUI_DLLAPI SpectrumSettingsVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<SpectrumSettings>
	{
		Q_OBJECT

		public:
			
		  ///Constructor
			SpectrumSettingsVisualizer(bool editable = false, QWidget* parent = 0);
		  
		public slots:
			
		  //Docu in base class
			void store();
		
		protected slots:
			
			///Undo the changes made in the GUI.
			void undo_();
	
		protected:  
			/// The date of this experiment
			QLineEdit* native_id_;
			/// The type of this experiment
	   	QComboBox* type_;
			/// The date of this experiment
			QTextEdit* comment_;
			
			//Docu in base class
			void update_();
	};
}
#endif
