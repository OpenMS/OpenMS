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


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASOutputFilesDialog.h>

namespace OpenMS 
{
	/**
		@brief Dialog which allows to specify the directory for the output files
		
		@ingroup TOPPAS_elements
		@ingroup Dialogs
	*/
	class OPENMS_GUI_DLLAPI TOPPASOutputFilesDialog
		: public QDialog,
			public Ui::TOPPASOutputFilesDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASOutputFilesDialog(const QString& dir_name);
			
			/// Returns the name of the directory
			QString getDirectory();
			
			/// Returns if the directory is valid (is a directory and writable)
			static bool dirNameValid(const QString& dir_name);
			
		public slots:
		
			/// Lets the user select the directory via a file dialog
			void showFileDialog();
		
		protected slots:
		
			/// Called when OK is pressed; checks if the selected file is valid
			void checkValidity_();
			
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASOUTPUTFILESDIALOG_H
