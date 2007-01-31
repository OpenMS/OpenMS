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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

/**
  @page TOPPView TOPPView
  
  TOPPView is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzData, mzXML, ANDI/MS
  and several other text-based file formats. It also supports viewing data from an OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC/MS map and a MS raw spectrum:
 	
 	\image html TOPPView.png

  <HR>
  
  Short description of the main features and options of TOPPView: 
  
  <b>Intensity display modes:</b>
  <BR>
  Intensity display modes determine the way peak intensities are displayed.
  <UL>
  <LI><b>Linear:</b> <BR> Normal display mode.
  <LI><b>Logarithmic:</b> <BR> Log10 of the intensity is displayed.
  <LI><b>Percentage:</b> <BR> In this display mode the intensities of each dataset are normalized with the maximum 
                         intensity of the dataset. This is espoecially useful to visualize several datasets that have
                         large intensity differences. When only one dataset is is opened it corresponds to the normal mode.
  <LI><b>Snap to maximum intensity:</b> <BR> In this mode the maxiumum currently displayed intensity is treated as if it was
                                        maxium overall intensity.
  </UL>

  <B>Action modes:</b>
  <BR>
  Action modes determine the mouse actions. Action modes not supported in the chosen spectrum display mode are displayed in gray.
  <UL>
  <LI><b>Zoom:</b> <BR> Allows displaying a specific data area.
  <LI><b>Translate:</b> <BR> In this mode it is possible to translate the visible area.
  <LI><b>Select:</b> <BR> The m/z, RT and intensity of a selected peak are displayed in the sttatus bar in this mode.
  <LI><b>Measure</b> <BR> This mode is used to determine the difference in m/z and RT, and intensity ratio of the seleced
                          peaks.
  </UL>

  <B>Context menu:</b>
  <BR>
  Some important features of TOPPView are only accessable through the context menu of the spectrum.
  <UL>
  <LI><b>Intensity distribution:</b> <BR> The displayed intensity range can be adapted by this dialog. Use the 
                                          left and right slider to set the displayed intensity range.
  <LI><b>Preferences:</b> <BR> The properties of each view (colors, ...) can be set through this menu. For each open window 
                               a properties page is displayed. Additionally the default settings of TOPPView can be set.                               
  </UL>

  <B>Open menu options:</b>
  <BR>
  The open menu of TOPPView offers several options which are explained here.
  <UL>
  <LI><b>Source:</b> <BR> Determines if the file system or database is browsed.
  <LI><b>Force file type:</b> <BR> Forces the file type of the chosen files. 
                                   Use this option only when the file extension is ambigous e.g. '.xml'.
  <LI><b>Open in:</b> <BR> Determines if the data of the chosen file(s) is displayed in the currently active spectrum
                           window or in a new window.
  <LI><b>Open map as:</b> <BR> Determines the display mode for HPLC-MS maps: '1D' or '2D'.
  <LI><b>Low intensity cutoff:</b> <BR> Suppresses displaying low intensity peaks by estimating the global noise level
                                        of the data. This is especially usefull for large datasets. Low intensity peaks
                                        are not removed from the data however. They can be displayed through the 
                                        'intensity distribution' menu of the context menu.
  </UL>
*/

//QT
#include <qapplication.h>
#include <qwindowsstyle.h>

//OpenMS
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>

using namespace OpenMS;
using namespace std;

//STL
#include <iostream.h>
#include <map.h>
#include <vector.h>

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
       << tool_name << " -- A viewer for MS data." << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options] [files]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  --help           Shows this help" << endl
			 << "  -ini <File>      Sets the INI file (default: ~/.TOPPView.ini)" << endl
			 << endl
			 << "To open several files in one window put a '+' in between the files." << endl
			 << "Example: 'TOPPView 1.dta + 2.dta + 3.dta'" << endl
			 << endl ;
}

int main( int argc, char ** argv )
{
	//list of all the valid options
	map<string,string> valid_options, valid_flags;
	valid_flags["--help"] = "help";
	valid_options["-ini"] = "ini";
	
	Param param;
	param.parseCommandLine(argc, argv, valid_options, valid_flags, "misc", "unkonwn");

	// '--help' given
	if (!(param.getValue("help").isEmpty()))
	{
		print_usage();
		return 0;
	}	

	// test if unknown options were given
	if (!param.getValue("unknown").isEmpty())
	{
		cout << "Unknown option '" << (string)(param.getValue("unknown")) << "' given. Aborting!" << endl;
		print_usage();
		return 1;
	}

#ifndef DEBUG_TOPP	
	try
	{
#endif
	  QApplication a( argc, argv );
	  TOPPViewBase* mw = TOPPViewBase::instance();
	  a.setMainWidget(mw);
	  if (!param.getValue("ini").isEmpty())
	  {
	  	mw->loadPreferences((String)param.getValue("ini"));
	  }
	  mw->setCaption( "TOPPView" );
	  mw->show();
	  
	  //load command line files
	  if (!param.getValue("misc").isEmpty())
	  {
	  	vector<String> filelist;
	  	bool several=((String)(param.getValue("misc"))).split(' ',filelist);
	  	if(!several)
	  	{
	  		filelist.push_back((String)(param.getValue("misc")));
	  	}
	  	mw->loadFiles(filelist.begin(),filelist.end());
	  }
	  
	  a.connect( &a, SIGNAL(lastWindowClosed()), &a, SLOT(quit()) );
	
	  int res = a.exec();
		mw->savePreferences();
	  return res;
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
	catch(Exception::Base& e)
	{
		cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
#endif
	
	return 0;
}

