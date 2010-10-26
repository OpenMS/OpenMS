// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------


#include <OpenMS/APPLICATIONS/INIFileEditorWindow.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/SYSTEM/File.h>

#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>

#ifdef OPENMS_WINDOWSPLATFORM
#	ifndef _WIN32_WINNT
#		define _WIN32_WINNT 0x0501 // Win XP (and above)
#	endif
#	include <Windows.h>
#endif


using namespace OpenMS;
using namespace std;

/**
	@page TOPP_INIFileEditor INIFileEditor

	@brief Can be used to visually edit INI files of TOPP tools.

	The values can be edited by double-clicking or pressing F2.

	The documentation of each value is shown in the text area on the bottom of the widget.

	@image html INIFileEditor.png

	More information about TOPPAS can be found in the @ref TOPP_tutorial.
*/

int main(int argc, const char** argv)
{

	Map<String,String> option_lists;
	Map<String,String> options;
	options["-print"] = "print";
	Map<String,String> flags;
	flags["--help"] = "help";
	Param param;
	param.parseCommandLine(argc, argv, options, flags, option_lists);

	//catch command line errors
	if (param.exists("help") //help requested
		  || argc>3 //too many arguments
		  || ( argc==3 && !param.exists("print")) //three argument but no -print
		  || (param.exists("print") && param.getValue("print")=="") //-print but no file given
		 )
	{
		cerr << endl
	       << "INIFileEditor -- An editor for OpenMS configuration files." << endl
	       << endl
	       << "Usage:" << endl
				 << " INIFileEditor [options] [file]" << endl
				 << endl
				 << "Options are:" << endl
		     << " --help         Shows this help and exits" << endl
				 << " -print <file>  Prints the content of the file to the command line and exits" << endl
				 << endl;
		return 0;
	}

	//print a ini file as text
	if (param.exists("print"))
	{
		Param data;
		try
		{
			data.load(param.getValue("print"));
			for (Param::ParamIterator it=data.begin(); it!=data.end(); ++it)
			{
				cout << it.getName() << " = " << it->value << endl;
			}
		}
		catch (Exception::BaseException &e)
		{
			LOG_ERROR << "Error while parsing file '" << param.getValue("print") << "'\n";
			LOG_ERROR << e << "\n";
		}

		return 0;
	}

	//Create window
	QApplication app(argc,const_cast<char**>(argv));

  //set plastique style unless windows / mac style is available
  if (QStyleFactory::keys().contains("windowsxp",Qt::CaseInsensitive))
  {
		app.setStyle("windowsxp");
  }
  else if (QStyleFactory::keys().contains("macintosh",Qt::CaseInsensitive))
  {
		app.setStyle("macintosh");
  }
  else if (QStyleFactory::keys().contains("plastique",Qt::CaseInsensitive))
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

#ifdef OPENMS_WINDOWSPLATFORM
  FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
  AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif

	editor_window.show();
	return app.exec();
}
