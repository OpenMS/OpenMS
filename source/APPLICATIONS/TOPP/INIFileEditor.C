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
//  License as editor_windowpublished by the Free Software Foundation; either
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


#include <OpenMS/APPLICATIONS/INIFileEditorWindow.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>

using namespace OpenMS;
using namespace std;

/**
	@page INIFileEditor INIFileEditor
	
	@brief Can be used to visually edit INI files of TOPP tools.
	
	@image html INIFileEditor.png
	
*/

int main(int argc, char** argv)
{
	if (argc>2)
	{
		cout << "Usage: " << argv[0] << " [file ]" << endl;
		return 0;
	}
	
	//Create window
	QApplication app(argc,argv);

	  //set plastique style unless windows / mac style is available
	  QStringList styles = QStyleFactory::keys();
	  
	  if (styles.contains("windowsxp",Qt::CaseInsensitive))
	  {
			app.setStyle("windowsxp");
	  }
	  else if (styles.contains("macintosh",Qt::CaseInsensitive))
	  {
			app.setStyle("macintosh");
	  }
	  else if (styles.contains("plastique",Qt::CaseInsensitive))
	  {
			app.setStyle("plastique");
	  }


	INIFileEditorWindow editor_window;
	
	//Open passed file
	if (argc==2)
	{
		//cout << "OPEN: "  << argv[1] << endl;
		editor_window.openFile(argv[1]);
	}
	
	editor_window.showMaximized();
	return app.exec();
}
