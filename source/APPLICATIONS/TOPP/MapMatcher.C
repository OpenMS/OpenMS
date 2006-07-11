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
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DMapMatcherRegression.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePair.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DFeaturePairVector.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>

#include <OpenMS/FORMAT/DGridFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>

#include <TOPPCommon.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapMatcher MapMatcher
	
	@brief Computes a transformation for a list of feature pairs.
	
	This is the second step in the map matching workflow. This application
	takes a list of feature pairs as computed by the FeatureMatcher and
	a grid (partially) covering the LC/MS map. For each grid cell, a 
	transformation is computed that maps the feature partners on each
	other. Currently, this transformation is linear.
	
	The output of this application is the list of grid cells with the
	estimated transformation.
	
		
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "MapMatcher";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- estimates a transformation for pairs of features in different LC/MS maps" << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << " -grid <file>   grid covering the map to be transformed (default read from INI file)" << endl
			 << " -pairs <file>  feature pairs (default read from INI file)" << endl
			 << " -q <float>  	 minimum quality of pairs considered (default read from INI file)" << endl
			 << " -out <file>  	 output file (default read from INI file)" << endl
			 << endl
			 << "Common TOPP options are:" << endl
			 << "  -ini <file>       TOPP INI file (default: TOPP.ini)" << endl
			 << "  -log <file>       log file (default: TOPP.log)" << endl
			 << "  -n <int>          instance number (default: 1)" << endl
			 << "  -d <level>        sets debug level (default: 0)" << endl
			 << "  --help            shows this help" << endl
       << "  --help-opt        shows help on the INI options accepted" << endl
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
	// log filestream (as long as the real logfile is not setermined yet)
	ofstream log;
	log.open("TOPP.log", ofstream::out | ofstream::app);
	
	//-------------------------------------------------------------
	// command line parsing
	//-------------------------------------------------------------
	
	//list of all the valid options
	map<string,string> valid_options;
	valid_options["-grid"] = "grid";
	valid_options["-pairs"] = "pairs";
	valid_options["-q"] = "min_quality";
	valid_options["-out"] = "out";
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
		param.load((string)(param.getValue("ini")));
		if (debug_level>0) log << Date::now() << " " << ini_location << " INI file: " << param.getValue("ini") << endl;
	
		
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
		if (debug_level>0) log << Date::now() << " " << ini_location << " log file: " << logfile << endl;
		log.close();
		log.open(logfile.c_str(), ofstream::out | ofstream::app);

		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
		//File names
		String gridfile, pairsfile, outfile;
		
		if (!(param.getValue("grid").isEmpty())) //from command line
		{
			gridfile = (string)(param.getValue("grid"));
		}
		else if (!(param.getValue(ini_location+"grid").isEmpty())) //from INI file
		{
			gridfile = (string)(param.getValue(ini_location+"grid"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " grid file: " << gridfile << endl;	
		
		if (!(param.getValue("pairs").isEmpty())) //from command line
		{
			pairsfile = (string)(param.getValue("pairs"));
		}
		else if (!(param.getValue(ini_location+"pairs").isEmpty())) //from INI file
		{
			pairsfile = (string)(param.getValue(ini_location+"pairs"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " pairs file: " << pairsfile << endl;	
	
		//determine output file name
		if (!(param.getValue("out").isEmpty())) //from command line
		{
			outfile = (string)(param.getValue("out"));
		}
		else if (!(param.getValue(ini_location+"out").isEmpty())) //from INI file
		{
			outfile = (string)(param.getValue(ini_location+"out"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " output file: " << outfile << endl;	
		
		//parameters
		double min_quality = 0;
			
		if (!(param.getValue("min_quality").isEmpty())) //from command line
		{
			min_quality = (double)(param.getValue("min_quality"));
		}
		else if (!(param.getValue(ini_location+"min_quality").isEmpty())) //from INI file
		{
			min_quality = (double)(param.getValue(ini_location+"min_quality"));
		}
		if (debug_level>1) log << Date::now() << " " << ini_location << " min_quality: " << min_quality << endl;	

		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		
		DGridFile grid_file;
		DGrid<2> the_grid;
		grid_file.load(gridfile,the_grid);
				
		DFeaturePairsFile pairs_file;
		DFeaturePairVector<2> pairs_vector;
		pairs_file.load(pairsfile,pairs_vector);
		
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		DMapMatcherRegression<2> map_matcher;
		map_matcher.setFeaturePairs(pairs_vector);
		map_matcher.setGrid(the_grid);
		map_matcher.setMinQuality(min_quality);
		
		map_matcher.estimateTransform();
		
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		
		DGrid<2> grid_with_transform;
		grid_with_transform = map_matcher.getGrid();
		
		grid_file.store(outfile,grid_with_transform);
	
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




