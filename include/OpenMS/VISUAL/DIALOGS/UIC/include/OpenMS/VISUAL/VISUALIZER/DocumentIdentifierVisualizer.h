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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_VISUALIZER_DOCUMENTIDENTIFIERVISUALIZER_H
#define OPENMS_VISUAL_VISUALIZER_DOCUMENTIDENTIFIERVISUALIZER_H

//OpenMS
#include <OpenMS/METADATA/DocumentIdentifier.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizer.h>
#include <OpenMS/VISUAL/VISUALIZER/BaseVisualizerGUI.h>

//QT

namespace OpenMS
{
	/**
		@brief Class that displays all meta information for DocumentIdentifier objects

		This class provides all functionality to view the meta information of an object of type DocumentIdentifier.
	*/
	class OPENMS_DLLAPI DocumentIdentifierVisualizer
		: public BaseVisualizerGUI,
			public BaseVisualizer<DocumentIdentifier>
	{
			Q_OBJECT

		public:

		   ///Constructor
			DocumentIdentifierVisualizer(bool editable = false, QWidget* parent = 0);

		public slots:

		  //Docu in base class
			void store();

		protected slots:

			///Undo the changes made in the GUI.
			void undo_();

		protected:
			///@name Edit fields and buttons
	    //@{
			QLineEdit* identifier_;
			QLineEdit* file_path_;
			QLineEdit* file_type_;
			//@}

			//Docu in base class
			void update_();
	};
}
#endif
