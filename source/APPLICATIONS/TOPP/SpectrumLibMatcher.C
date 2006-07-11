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

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <TOPPCommon.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "TOPP_Skeleton";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
			 << tool_name << " usage: [-in <file>] [-out <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]" << endl
			 << "  -in <file>   input file in MzData format (default read from INI file)" << endl
			 << "  -out <file>  output file in analysisXML format (default read from INI file)" << endl
			 << "  -ini <file>  TOPP INI file (default: TOPP.ini)" << endl
			 << "  -log <file>  log file (default: TOPP.log)" << endl
			 << "  -n <int>     instance number (default: 1)" << endl
			 << "  -d <level>   sets debug level (default: 0)" << endl
			 << "  --help       shows this help" << endl
			 << endl
		;
}

//-------------------------------------------------------------
// main program
//-------------------------------------------------------------

int main( int argc, char ** argv )
{
	// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
	String ini_location;
	// path to the log file
	String logfile = "";
	// debug level
	int debug_level = 0;
	// log filestream
	ofstream log;
	log.open("TOPP.log", ofstream::out | ofstream::app);
	
	//-------------------------------------------------------------
	// command line parsing
	//-------------------------------------------------------------
	
	//list of all the valid options
	map<string,string> valid_options;
	valid_options["-out"] = "out";
	valid_options["-in"] = "in";
	valid_options["-ini"] = "ini";
	valid_options["-log"] = "log";
	valid_options["-n"] = "instance";
	valid_options["-d"] = "debug";
	valid_options["--help"] = "help";
	//for debugging
	valid_options["unknown"] = "unknown";
	valid_options["misc"] = "misc";
	
	Param param;
	param.parseCommandLine(argc, argv, valid_options);

	//-------------------------------------------------------------
	// read debug level from command line if set 
	//-------------------------------------------------------------
	if (!param.getValue("debug").isEmpty())
	{
		debug_level = (int)(param.getValue("debug"));
	}
	
	//-------------------------------------------------------------
	// determine instance number
	//-------------------------------------------------------------

	if (param.getValue("instance").isEmpty())
	{
		param.setValue("instance",1);
	}
	ini_location = String(tool_name) + ":" + param.getValue("instance").toString() + ":";
	if (debug_level>0) cout << "Instance number: " << param.getValue("instance") << endl;

	//-------------------------------------------------------------
	// check command line options
	//-------------------------------------------------------------
	
	// '--help' given
	if (!(param.getValue("help").isEmpty()))
	{
		print_usage();
		return OK;
	}
	
	// test if unknown options were given
	if (!param.getValue("unknown").isEmpty())
	{
		log << ini_location << " Unknown option '" << (string)(param.getValue("unknown")) << "' given. Aborting!" << endl;
		print_usage();
		return ILLEGAL_PARAMETERS;
	}
	
	// test if unknown text argument were given (we do not use them)
	if (!param.getValue("misc").isEmpty())
	{
		log << ini_location << " Trailing text argument '" << (string)(param.getValue("misc")) << "' given. Aborting!" << endl;
		print_usage();
		return ILLEGAL_PARAMETERS;
	}

	try
	{
	
	//-------------------------------------------------------------
	// loading INI file
	//-------------------------------------------------------------
	if (param.getValue("ini").isEmpty())
	{
		param.setValue("ini",string("TOPP.ini"));
	}
	param.load((string)(param.getValue("ini")));
	if (debug_level>0) cout << "INI file: " << param.getValue("ini") << endl;

	
	//-------------------------------------------------------------
	// determine and open log file
	//-------------------------------------------------------------
	if (!(param.getValue("log").isEmpty()))
	{
		logfile = (string)(param.getValue("log"));
	}
	if (param.getValue("log").isEmpty() && !(param.getValue(ini_location+"log").isEmpty()))
	{
		logfile = (string)(param.getValue(ini_location+"log"));
	}
	if (param.getValue("log").isEmpty() && !(param.getValue("common:log").isEmpty()))
	{
		logfile = (string)(param.getValue("common:log"));
	}
	if (debug_level>0) cout << "log file: " << logfile << endl;
	log.close();
	log.open(logfile.c_str(), ofstream::out | ofstream::app);
	
	//-------------------------------------------------------------
	// calculations
	//-------------------------------------------------------------
	
	//...
	
	//-------------------------------------------------------------
	// writing files
	//-------------------------------------------------------------
	
	//...
	
	}
	catch(Exception::UnableToCreateFile& e)
	{
		log << ini_location << " Error: Unable to write file (" << e.what() <<")"<< endl;
		return CANNOT_WRITE_OUTPUT_FILE;
	}	
	catch(Exception::FileNotFound& e)
	{
		log << ini_location << " Error: File not found (" << e.what() <<")"<< endl;
		return INPUT_FILE_NOT_FOUND;
	}
	catch(Exception::ParseError& e)
	{
		log << ini_location << " Error: Unable to read file (" << e.what() <<")"<< endl;
		return INPUT_FILE_CORRUPT;
	}
	catch(Exception::Base& e)
	{
		log << ini_location << " Error: Unknown error (" << e.what() <<")"<< endl;
		return UNKNOWN_ERROR;
	}
  
	log.close();
	
	return OK;
}




