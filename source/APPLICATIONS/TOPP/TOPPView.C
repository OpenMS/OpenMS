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
 	
 	@image html TOPPView.png

  <HR>
  
  Short description of the main features and options of TOPPView: 
  
  <b>Intensity display modes:</b>
  <BR>
  Intensity display modes determine the way peak intensities are displayed.
  <UL>
  <LI><b>Linear:</b> <BR> Normal display mode.
  <LI><b>Logarithmic:</b> <BR> Log10 of the intensity is displayed.
  <LI><b>Percentage:</b> <BR> In this display mode the intensities of each dataset are normalized with the maximum 
                         intensity of the dataset. This is especially useful to visualize several datasets that have
                         large intensity differences. When only one dataset is is opened it corresponds to the normal mode.
  <LI><b>Snap to maximum intensity:</b> <BR> In this mode the maxiumum currently displayed intensity is treated as if it was
                                        maxium overall intensity.
  </UL>

  <B>Action modes:</b>
  <BR>
  Action modes determine the mouse actions. Action modes not supported in the chosen spectrum display mode are displayed in gray.
  <UL>
  <LI><b>Zoom + Translate:</b> <BR> Allows zooming to a specific data area.
  																	When pressing the CTRL key you can translate the displayed area.
  <LI><b>Select + Measure:</b> <BR> The m/z, RT and intensity of a selected peak are displayed in the sttatus bar in this mode.
  														 <BR> When pressing the CTRL key you can determine the difference in m/z and RT, and intensity ratio of 
  														 			the seleced peaks.
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
  <LI><b>Open map as:</b> <BR> Determines the display mode for HPLC-MS maps: '2D' or '3D (experimental)' .
  <LI><b>Low intensity cutoff:</b> <BR> Suppresses displaying low intensity peaks by estimating the global noise level
                                        of the data. This is especially usefull for large datasets. Low intensity peaks
                                        are not removed from the data however. They can be displayed through the 
                                        'Intensity distribution' tool in the 'Layer' menu.
  </UL>

  <B>Context menu options:</b>
  <UL>
  	<LI> Show the currently displayed data in 3D
  	<LI> Extract a scan in 1D
  	<LI> Edit meta data of a scan/feature
  </UL>
*/

//QT
#include <QtGui/QApplication>
#include <QtGui/QStyleFactory>

//OpenMS
#include <OpenMS/APPLICATIONS/TOPPViewBase.h>

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
			 << "To open several files in one window put a '+' in between the files." << endl
			 << "Example: 'TOPPView 1.dta + 2.dta + 3.dta'" << endl
			 << endl ;
}

int main( int argc, char ** argv )
{
	//list of all the valid options
	map<String,String> valid_options, valid_flags;
	valid_flags["--help"] = "help";
	valid_options["-ini"] = "ini";
	
	Param param;
	param.parseCommandLine(argc, argv, valid_options, valid_flags, "misc", "unkonwn");

	// '--help' given
	if (param.exists("help"))
	{
		print_usage();
		return 0;
	}	

	// test if unknown options were given
	if (param.exists("unknown"))
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
	  
	  //set plastique style unless windows / mac style is available
	  QStringList styles = QStyleFactory::keys();
	  
	  if (styles.contains("windowsxp",Qt::CaseInsensitive))
	  {
			a.setStyle("windowsxp");
	  }
	  else if (styles.contains("macintosh",Qt::CaseInsensitive))
	  {
			a.setStyle("macintosh");
	  }
	  else if (styles.contains("plastique",Qt::CaseInsensitive))
	  {
			a.setStyle("plastique");
	  }
	  
	  TOPPViewBase* mw = new TOPPViewBase();
	  if (param.exists("ini"))
	  {
	  	mw->loadPreferences((String)param.getValue("ini"));
	  }
	  mw->show();
	  
	  //load command line files
	  if (param.exists("misc"))
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
	catch(Exception::Base& e)
	{
		cout << String("Error: Unexpected error (") << e.what() << ")" << endl << "Code location: " << e.getFile() << ":" << e.getLine() << endl;
		return 1;
	}
#endif
	
	return 1;
}

