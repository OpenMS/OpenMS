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
// $Maintainer: Marc Sturm $
// $Authors: $
// --------------------------------------------------------------------------

#ifndef OPENMS_APPLICATIONS_INIFILEEDITORWINDOW_H
#define OPENMS_APPLICATIONS_INIFILEEDITORWINDOW_H

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/DATASTRUCTURES/Param.h>

#include <QtGui/QMainWindow>

class QToolBar;
class QAction;
class QString;
class QFileDialog;

namespace OpenMS
{
	/**
		@brief shows the ParamEditor widget in a QMainWindow with a toolbar
	*/
	class OPENMS_DLLAPI INIFileEditorWindow
		: public QMainWindow
	{
		Q_OBJECT
		
		public:
			/// menu is created here
			INIFileEditorWindow(QWidget *parent = 0);
			/// when user closes window a message box asks the user if he wants to save
			void closeEvent(QCloseEvent *event);
			
		public slots:
			///loads the xml-file into a Param object and loads Param into ParamEditor
			bool openFile(const String& filename="");
			/// saves the users changes in a xml-file if the Param object is valid
			bool saveFile();
			/// like saveFile but with a file dialog to choose a filename
			bool saveFileAs();
			/// if the user changes data in ParamEditor the title shows a '*'
			void updateWindowTitle(bool);
		
		private:
			/// ParamEditor object for visualization
			ParamEditor* editor_;
			/// Param object for storing data
			Param param_;
			/// filename of xml-file to store the Param object
			QString filename_;
			/// path used as next default location of the load/store dialogs
			String current_path_;
	};
}

#endif //OPENMS_APPLICATIONS_INIFILEEDITORWINDOW_H
