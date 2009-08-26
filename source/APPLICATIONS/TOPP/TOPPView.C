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

/**
  @page TOPP_TOPPView TOPPView
  
  TOPPView is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzML, mzData, mzXML, ANDI/MS
  and several other file formats. It also supports viewing data from an OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC/MS map and a MS raw spectrum:
 	
 	@image html TOPPView.png
	
	More information about TOPPView can be found in the @ref TOPP_tutorial.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_TOPPView.cli
*/

//QT
#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>
#include <QtGui/QSplashScreen>

//OpenMS
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>
#include <OpenMS/SYSTEM/StopWatch.h>


using namespace OpenMS;
using namespace std;

//STL
#include <iostream>
#include <map>
#include <vector>

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "TOPPView";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- A viewer for mass spectrometry data." << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options] [files]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  --help           Shows this help" << endl
			 << "  -ini <File>      Sets the INI file (default: ~/.TOPPView.ini)" << endl
			 << endl
			 << "Hints:" << endl
			 << " - To open several files in one window put a '+' in between the files." << endl
			 << " - '@bw' after a map file displays the dots in a white to black gradient." << endl
			 << " - '@bg' after a map file displays the dots in a grey to black gradient." << endl
			 << " - '@b'  after a map file displays the dots in black." << endl
			 << " - '@r'  after a map file displays the dots in red." << endl
			 << " - '@g'  after a map file displays the dots in green." << endl
			 << " - '@m'  after a map file displays the dots in magenta." << endl
			 << " - Example: 'TOPPView 1.mzML + 2.mzML @bw + 3.mzML @bg'" << endl
			 << endl ;
}

int main( int argc, const char** argv )
{
	//list of all the valid options
	Map<String,String> valid_options, valid_flags, option_lists;
	valid_flags["--help"] = "help";
	valid_options["-ini"] = "ini";
	
	Param param;
	param.parseCommandLine(argc, argv, valid_options, valid_flags, option_lists);

	// '--help' given
	if (param.exists("help"))
	{
		print_usage();
		return 0;
	}	

	// test if unknown options were given
	if (param.exists("unknown"))
	{    
    // if TOPPView is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if(!(param.getValue("unknown").toString().hasSubstring("-psn") && !param.getValue("unknown").toString().hasSubstring(", ")))
    {
      cout << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage();
      return 1;
    }
	}

#ifndef DEBUG_TOPP	
	try
	{
#endif
	  QApplication a( argc, const_cast<char**>(argv));
	  a.connect( &a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()) );
	  		
	  //set plastique style unless windows / mac style is available
	  if (QStyleFactory::keys().contains("windowsxp",Qt::CaseInsensitive))
	  {
			a.setStyle("windowsxp");
	  }
	  else if (QStyleFactory::keys().contains("macintosh",Qt::CaseInsensitive))
	  {
			a.setStyle("macintosh");
	  }
	  else if (QStyleFactory::keys().contains("plastique",Qt::CaseInsensitive))
	  {
			a.setStyle("plastique");
	  }

	  TOPPViewBase* mw = new TOPPViewBase();
	  mw->show();

		// Create the splashscreen that is displayed while the application loads
		QSplashScreen* splash_screen = new QSplashScreen(QPixmap(":/TOPPView_Splashscreen.png"));
		splash_screen->show();
		splash_screen->showMessage("Loading parameters");
		QApplication::processEvents();
		StopWatch stop_watch;
		stop_watch.start();

	  if (param.exists("ini"))
	  {
	  	mw->loadPreferences((String)param.getValue("ini"));
	  }

	  //load command line files
	  if (param.exists("misc"))
	  {
	  	mw->loadFiles((StringList)(param.getValue("misc")), splash_screen);
	  }

		// We are about to show the application. 
		// Proper time to  remove the splashscreen, if at least 1.5 seconds have passed...
		while(stop_watch.getClockTime()<1.5) {/*wait*/};
		stop_watch.stop();
		splash_screen->close();
		delete splash_screen;
		
	  int result = a.exec();
	  delete(mw);
	  return result;
#ifndef DEBUG_TOPP
	}
	//######################## ERROR HANDLING #################################

	catch(Exception::UnableToCreateFile& e)
	{
		cout << String("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}	
	catch(Exception::FileNotFound& e)
	{
		cout << String("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
	catch(Exception::FileNotReadable& e)
	{
		cout << String("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
	catch(Exception::FileEmpty& e)
	{
		cout << String("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
	catch(Exception::ParseError& e)
	{
		cout << String("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
	catch(Exception::InvalidValue& e)
	{
		cout << String("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
	catch(Exception::BaseException& e)
	{
		cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
#endif
	
	return 1;
}

