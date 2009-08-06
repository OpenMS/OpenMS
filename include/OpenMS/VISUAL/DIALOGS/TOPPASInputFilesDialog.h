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


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASInputFilesDialog.h>

namespace OpenMS 
{
	/**
		@brief Dialog which allows to specify a list of input files
		
		@ingroup TOPPAS_elements
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI TOPPASInputFilesDialog
		: public QDialog,
			public Ui::TOPPASInputFilesDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASInputFilesDialog(const QStringList& list);
			
			/// Stores the list of all filenames in the list widget in @p files
			void getFilenames(QStringList& files);
			
		public slots:
		
			/// Lets the user select files via a file dialog
			void showFileDialog();
			/// Removes all currently selected files from the list
			void removeSelected();
			/// Shows a TOPPASInputFileDialog which edits the current item
			void editCurrentItem();
		
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
