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


#ifndef OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILEDIALOG_H
#define OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILEDIALOG_H

#include <OpenMS/config.h>
#include <OpenMS/VISUAL/DIALOGS/UIC/ui_TOPPASInputFileDialog.h>

namespace OpenMS 
{
	/**
		@brief Dialog which allows to specify an input file
		
		@ingroup TOPPAS_elements
		@ingroup Dialogs
	*/
	class OPENMS_DLLAPI TOPPASInputFileDialog
		: public QDialog,
			public Ui::TOPPASInputFileDialogTemplate
	{
		Q_OBJECT
				
		public:
			
			/// Constructor
			TOPPASInputFileDialog(const QString& file_name);
			
			/// Returns the filename
			QString getFilename();
			
			/// Returns if the file name is valid (exists, is readable, and is not a directory)
			static bool fileNameValid(const QString& file_name);
			
		public slots:
		
			/// Lets the user select the file via a file dialog
			void showFileDialog();
		
		protected slots:
		
			/// Called when OK is pressed; checks if the selected file is valid
			void checkValidity_();
			
		protected:
		
			/// The parent
			QObject* parent_;
	};
	
}
#endif // OPENMS_VISUAL_DIALOGS_TOPPASINPUTFILESDIALOG_H
