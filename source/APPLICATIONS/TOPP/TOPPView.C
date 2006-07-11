// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
  @defgroup TOPPView TOPPView - A MS data viewer
  
  @brief TOPPView documentation
  
  Specview is a viewer for MS and HPLC-MS data. It can be used to inspect files in mzData, mzXML, ANDI/MS
  and several other text-based file formats. It also supports viewing data from an OpenMS database.
  The following figure shows two instances of TOPPView displaying a HPLC/MS map and a MS raw spectrum:
 	
 	\image html TOPPView.png
 	
  Specview can be found in the <tt>OpenMS/source/APPLICATIONS/TOPP/</tt> directory.
  
  Use the command <tt>'make TOPPView'</tt> to build it after you built OpenMS.
  
  @ingroup Applications
*/

//QT
#include <qapplication.h>
#include <qwindowsstyle.h>

//OpenMS
#include <OpenMS/VISUAL/SpectrumMDIWindow.h>

using namespace OpenMS;

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
			 << "  --help            shows this help" << endl
			 << "  --ini <File>      Sets the INI file (default: ~/.TOPPView.ini)" << endl
			 << endl ;
}

int main( int argc, char ** argv )
{
	//list of all the valid options
	map<string,string> valid_options, valid_flags;
	valid_flags["--help"] = "help";
	valid_options["--ini"] = "ini";
	
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
	
	try
	{
	  QApplication a( argc, argv );
	  SpectrumMDIWindow* mw = SpectrumMDIWindow::instance();
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
	}
	catch(Exception::Base& e)
	{
		cout << "Error: Unexpected error (" << e.what() <<")"<< endl;
		return 1;
	}
	
	return 0;
}

