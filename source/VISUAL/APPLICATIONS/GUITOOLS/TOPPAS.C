// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// $Authors: Johannes Junker, Chris Bielow $
// --------------------------------------------------------------------------

/**
	@page TOPP_TOPPAS TOPPAS
	
	@brief An assistant for GUI-driven TOPP workflow design.
	
  TOPPAS allows to create, edit, open, save, and run TOPP workflows. Pipelines
  can be created conveniently in a GUI by means of mouse interactions. The
  parameters of all involved tools can be edited within the application
  and are also saved as part of the pipeline definition in the @em .toppas file.
  Furthermore, TOPPAS interactively performs validity checks during the pipeline
  editing process, in order to make it more difficult to create an invalid workflow.
  Once set up and saved, a workflow can also be run without the GUI using
  the @em ExecutePipeline TOPP tool.
  
  The following figure shows a simple example pipeline that has just been created
  and executed successfully:
  
  @image html TOPPAS_simple_example.png
  
  More information about TOPPAS can be found in the @ref TOPPAS_tutorial.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_TOPPAS.cli
*/

//QT
#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>
#include <QtGui/QSplashScreen>
#include <QtCore/QDir>

//OpenMS
#include <OpenMS/VISUAL/APPLICATIONS/TOPPASBase.h>
#include <OpenMS/SYSTEM/StopWatch.h> 
#include <OpenMS/CONCEPT/LogStream.h>

using namespace OpenMS;
using namespace std;

//STL
#include <iostream>
#include <map>
#include <vector>

#ifdef OPENMS_WINDOWSPLATFORM
#	ifndef _WIN32_WINNT
#		define _WIN32_WINNT 0x0501 // Win XP (and above)
#	endif
#	include <Windows.h>
#endif


//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "TOPPAS";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage(Logger::LogStream& stream = Log_info)
{
	stream << "\n"
         << tool_name << " -- An assistant for GUI-driven TOPP workflow design." << "\n"
         << "\n"
         << "Usage:" << "\n"
			   << " " << tool_name << " [options] [files]" << "\n"
			   << "\n"
			   << "Options are:" << "\n"
			   << "  --help           Shows this help" << "\n"
			   << "  -ini <File>      Sets the INI file (default: ~/.TOPPAS.ini)" << "\n"
			   << endl ;
}

int main( int argc, const char** argv )
{
	//list of all the valid options
	Map<String,String> valid_options, valid_flags, option_lists;
	valid_flags["--help"] = "help";
	valid_options["-ini"] = "ini";
	//invalid, but keep for now in order to inform users where to find this functionality now
	valid_options["-execute"] = "execute";
	valid_options["-out_dir"] = "out_dir";
	
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
    // if TOPPAS is packed as Mac OS X bundle it will get a -psn_.. parameter by default from the OS
    // if this is the only unknown option it will be ignored .. maybe this should be solved directly
    // in Param.h
    if(!(param.getValue("unknown").toString().hasSubstring("-psn") && !param.getValue("unknown").toString().hasSubstring(", ")))
    {
      LOG_ERROR << "Unknown option(s) '" << param.getValue("unknown").toString() << "' given. Aborting!" << endl;
      print_usage(Log_error);
      return 1;
    }
	}

	try
	{

	 	if (param.exists("execute") || param.exists("out_dir"))
		{
			LOG_ERROR << "The parameters '-execute' and '-out_dir' are not valid anymore. This functionality has been moved to the ExecutePipeline tool." << endl; 
			return 1;
		}
		
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

	  TOPPASBase* mw = new TOPPASBase();
	  mw->show();

		// Create the splashscreen that is displayed while the application loads
		QSplashScreen* splash_screen = new QSplashScreen(QPixmap(":/TOPPAS_Splashscreen.png"));
		splash_screen->show();
		splash_screen->showMessage("Loading parameters");
		QApplication::processEvents();
		StopWatch stop_watch;
		stop_watch.start();

	  if (param.exists("ini"))
	  {
	  	mw->loadPreferences((String)param.getValue("ini"));
	  }

 	  if (param.exists("misc"))
 	  {
 	  	mw->loadFiles((StringList)(param.getValue("misc")), splash_screen);
 	  }
		else
		{ // remember this new window as obsolete once a real workflow is loaded without this window being touched
      // if this is not desired, simply call newPipeline() without arguments
			mw->newPipeline(mw->IDINITIALUNTITLED);
		}

		// We are about to show the application. 
		// Proper time to  remove the splash screen, if at least 1.5 seconds have passed...
		while(stop_watch.getClockTime()<1.5) {/*wait*/};
		stop_watch.stop();
		splash_screen->close();
		delete splash_screen;

#ifdef OPENMS_WINDOWSPLATFORM
    FreeConsole(); // get rid of console window at this point (we will not see any console output from this point on)
    AttachConsole(-1); // if the parent is a console, reattach to it - so we can see debug output - a normal user will usually not use cmd.exe to start a GUI)
#endif
		
	  int result = a.exec();
	  delete(mw);
	  return result;
	}
	//######################## ERROR HANDLING #################################

	catch(Exception::UnableToCreateFile& e)
	{
		cout << String("Error: Unable to write file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}	
	catch(Exception::FileNotFound& e)
	{
		cout << String("Error: File not found (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	catch(Exception::FileNotReadable& e)
	{
		cout << String("Error: File not readable (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	catch(Exception::FileEmpty& e)
	{
		cout << String("Error: File empty (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	catch(Exception::ParseError& e)
	{
		cout << String("Error: Unable to read file (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	catch(Exception::InvalidValue& e)
	{
		cout << String("Error: Invalid value (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	catch(Exception::BaseException& e)
	{
		cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
	}
	
	return 1;
}

