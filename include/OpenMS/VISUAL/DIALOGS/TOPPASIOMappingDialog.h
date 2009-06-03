// -*- Mode: C++; tab-width: 2; -*-
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
		@brief Dialog which allows to specify a list of input files
		
		@ingroup Dialogs
	*/
	class TOPPASIOMappingDialog
		: public QDialog,
			public Ui::TOPPASIOMappingDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASIOMappingDialog(TOPPASEdge* parent);
			
		public slots:
		
		protected:
		
		enum EdgeType_
		{
			ET_FILE_TO_TOOL,
			ET_LIST_TO_TOOL,
			ET_TOOL_TO_FILE,
			ET_TOOL_TO_LIST,
			ET_TOOL_TO_TOOL,
			ET_INVALID
		};
		
		///@name reimplemented Qt events
		//@{
		void resizeEvent(QResizeEvent* /*event*/);
		
		/// Determines the type of this edge
		void determineEdgeType_();
		/// Fills the table
		void fillTable_();
		
		/// The edge we are configuring
		TOPPASEdge* edge_;
		/// The type of this edge
		EdgeType_ edge_type_;
		/// The output parameters of the source tool
		QVector<TOPPASToolVertex::IOInfo> source_output_files_;
		/// The input parameters of the target tool
		QVector<TOPPASToolVertex::IOInfo> target_input_files_;
		
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
