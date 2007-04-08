// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stefan Rink $
// --------------------------------------------------------------------------

#ifndef INIFILEEDITORWINDOW_H
#define INIFILEEDITORWINDOW_H

#include <OpenMS/VISUAL/ParamEditor.h>
#include <OpenMS/FORMAT/Param.h>

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

	class INIFileEditorWindow
		: public QMainWindow
	{
		Q_OBJECT
		
		public:
			 INIFileEditorWindow(QWidget *parent = 0);
		
		public slots:
			bool openFile(const String& filename="");
			bool saveFile();
			bool saveFileAs();
		
		private:
			QToolBar* toolbar_;
			ParamEditor* editor_;
			Param param_;
			QString filename_;
	};
}

#endif
