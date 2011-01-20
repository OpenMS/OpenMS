// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Johannes Junker $
// $Authors: Johannes Junker $
// --------------------------------------------------------------------------


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASIOMAPPINGDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASIOMAPPINGDIALOG_H

#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASIOMappingDialog.h>
#include <OpenMS/VISUAL/TOPPASToolVertex.h>

#include <QtCore/QVector>


namespace OpenMS 
{
	class TOPPASEdge;
	
	/**
		@brief Dialog which allows to configure the input/output parameter mapping of an edge.
		
		This dialog allows to select an output parameter of the source vertex and an input
		parameter of the target vertex. Only valid selections are allowed, i.e. the type
		(file or list of files) and at least one valid file type of either vertex must match. 
		
		@ingroup TOPPAS_elements
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI TOPPASIOMappingDialog
		: public QDialog,
			public Ui::TOPPASIOMappingDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASIOMappingDialog(TOPPASEdge* parent);
			
		public slots:
	
			/// Called instead of exec() after edge is constructed (in order to avoid showing the dialog if not necessary)
			int firstExec();

		protected:
		
			/// Fills the table
			void fillComboBoxes_();
			
			/// The edge we are configuring
			TOPPASEdge* edge_;
			
			/// Vector storing the mapping of the target input combobox indices to param indices of edges
			QVector<int> target_input_param_indices_;
		
		protected slots:
		
			/// Called when OK is pressed; checks if the selected parameters are valid
			void checkValidity_();
		
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
