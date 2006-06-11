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
// $Id: FeaturePairSplitter.C,v 1.1 2006/03/16 12:47:10 groepl Exp $
// $Author: groepl $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <TOPPCommon.h>

#include <map>
#include <iostream>
#include <fstream>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------


/**
	@page FeaturePairSplitter FeaturePairSplitter
	
	@brief Splits a feature pair file into two feature files.

	This is just a small utility.  The features are copied from the pairs.  The
	relative order of features is preserved.  For example the first two features
	of each output file belong to each other, then the second two, and so on.
	However, the <i>quality</i> information of the feature pairs is lost.  (Thus
	you cannot restore the pairs from the two files!)
	
  The input file is parsed by DFeaturePairsFile; a typical file name extension
  would be ".fp".

  The two output files are written by DFeatureMapFile; a typical file name
  extension would be '.feat'.
	
	@note All options needed to operate this tool can be given in the command
	line, so no INI file is needed.
	
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "FeaturePairSplitter";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- split a feature pairs file into two feature files." << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>        input file" << endl
			 << "  -out1 <file>      first output file" << endl
			 << "  -out2 <file>      second output file" << endl
			 << "Common TOPP options are:" << endl
			 << "  -ini <file>       TOPP INI file (default: TOPP.ini)" << endl
			 << "  -log <file>       log file (default: TOPP.log)" << endl
			 << "  -n <int>          instance number (default: 1)" << endl
			 << "  -d <level>        sets debug level (default: 0)" << endl
			 << "  --help            shows this help" << endl
			 << endl ;
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
	valid_options["-in"] = "in";
	valid_options["-out1"] = "out1";
	valid_options["-out2"] = "out2";
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
	if (debug_level>0) log << Date::now() << " " << ini_location << " Instance number: " << param.getValue("instance") << endl;

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
		log << Date::now() << " " << ini_location << " Unknown option '" << (string)(param.getValue("unknown")) << "' given. Aborting!" << endl;
		cout << "Unknown option '" << (string)(param.getValue("unknown")) << "' given. Aborting!" << endl;
		print_usage();
		return ILLEGAL_PARAMETERS;
	}
	
	// test if unknown text argument were given (we do not use them)
	if (!param.getValue("misc").isEmpty())
	{
		log << Date::now() << " " << ini_location << " Trailing text argument '" << (string)(param.getValue("misc")) << "' given. Aborting!" << endl;
		cout << "Trailing text argument '" << (string)(param.getValue("misc")) << "' given. Aborting!" << endl;
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
		if (debug_level>0) log << Date::now() << " " << ini_location << " INI file: " << param.getValue("ini") << endl;
		//Try to load ini file
		try
		{
			param.load((string)(param.getValue("ini")));
		}
		catch(Exception::FileNotFound& e)
		{
			if (debug_level>0) log << Date::now() << " " << ini_location << " INI file not found!" << endl;
		}
		
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
		if (param.getValue("log").isEmpty() && param.getValue("common:log").isEmpty() && param.getValue(ini_location+"log").isEmpty())
		{
			logfile = "TOPP.log";
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " log file: " << logfile << endl;
		log.close();
		log.open(logfile.c_str(), ofstream::out | ofstream::app);
		
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
		
		//File names
		String in, out1, out2;
		
		
		// determine input file name and type
		if (!(param.getValue("in").isEmpty())) //from command line
		{
			in = (string)(param.getValue("in"));
		}
		else if (!(param.getValue(ini_location+"in").isEmpty())) //from INI file
		{
			in = (string)(param.getValue(ini_location+"in"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " input file: `" << in << '\'' << endl;	
		
		// determine first output file name and type
		if (!(param.getValue("out1").isEmpty())) //from command line
		{
			out1 = (string)(param.getValue("out1"));
		}
		else if (!(param.getValue(ini_location+"out1").isEmpty())) //from INI file
		{
			out1 = (string)(param.getValue(ini_location+"out1"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " first output file: `" << out1 << '\'' << endl;	

		// determine second output file name and type
		if (!(param.getValue("out2").isEmpty())) //from command line
		{
			out2 = (string)(param.getValue("out2"));
		}
		else if (!(param.getValue(ini_location+"out2").isEmpty())) //from INI file
		{
			out2 = (string)(param.getValue(ini_location+"out2"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " second output file: `" << out2 << '\'' << endl;	

		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
	
		//load input file data.
		DFeaturePairVector<2> feature_pairs;
		DFeaturePairsFile feature_pairs_file;
		feature_pairs_file.load(in,feature_pairs);

		//-------------------------------------------------------------
		// Do the transformation, create the feature maps.
		//-------------------------------------------------------------
		DFeatureMap<2> first_feature_map, second_feature_map;
		for ( DFeaturePairVector<2>::ConstIterator iter = feature_pairs.begin();
					iter != feature_pairs.end();
					++iter
				)
		{
			first_feature_map.push_back(iter->getFirst());
			second_feature_map.push_back(iter->getSecond());
		}

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
	
		DFeatureMapFile f;
		f.store(out1,first_feature_map);
		f.store(out2,second_feature_map);

	}
	catch(Exception::UnableToCreateFile& e)
	{
		cout << "Error: Unable to write file (" << e.what() <<")"<< endl;
		log << Date::now() << " " << ini_location << " Error: Unable to write file (" << e.what() <<")"<< endl;
		return CANNOT_WRITE_OUTPUT_FILE;
	}	
	catch(Exception::FileNotFound& e)
	{
		cout << "Error: File not found (" << e.what() <<")"<< endl;
		log << Date::now() << " " << ini_location << " Error: File not found (" << e.what() <<")"<< endl;
		return INPUT_FILE_NOT_FOUND;
	}
	catch(Exception::ParseError& e)
	{
		cout << "Error: Unable to read file (" << e.what() <<")"<< endl;
		log << Date::now() << " " << ini_location << " Error: Unable to read file (" << e.what() <<")"<< endl;
		return INPUT_FILE_CORRUPT;
	}
	catch(Exception::Base& e)
	{
		cout << "Error: Unexpected error (" << e.what() <<")"<< endl;
		log << Date::now() << " " << ini_location << " Error: Unexpected error (" << e.what() <<")"<< endl;
		return UNKNOWN_ERROR;
	}
  
	log.close();
	
	return OK;
}




