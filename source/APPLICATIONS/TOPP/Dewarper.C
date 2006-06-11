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
// $Id: Dewarper.C,v 1.7 2006/03/28 08:03:30 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/DMapDewarper.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGrid.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/DGridCell.h>

#include <OpenMS/FORMAT/DGridFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>

#include <OpenMS/KERNEL/DFeatureMap.h>

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
	@page Dewarper Dewarper
	
	@brief Dewarps a feature map by applying a transform to
	the coordinates of each feature.
	
	The dewarping is the last and optional step in a map
	matching workflow. The transform was computed in the map matching
	step of the workflow. Currently, we use a piecewise
	linear transform, but others can be implemented easily. This
	module simply applies this transform to the coordinates
	of each feature contained in the corresponding grid cells.
		
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "MapDewarper";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- dewarps a feature map" << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -grid <file>   grid covering the map to be transformed (default read from INI file)" << endl
			 << "  -feat <file>   feature pairs (default read from INI file)" << endl
			 << "  -out <file>  	dewarped feature map (default read from INI file)" << endl
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
	valid_options["-out"] = "out";
	valid_options["-grid"] = "grid";
	valid_options["-feat"] = "feat";
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
	std::cout << "Debug level: "  << debug_level << std::endl;
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
		String gridfile, features_file, outfile;
		
		if (!(param.getValue("grid").isEmpty())) //from command line
		{
			gridfile = (string)(param.getValue("grid"));
		}
		else if (!(param.getValue(ini_location+"grid").isEmpty())) //from INI file
		{
			gridfile = (string)(param.getValue(ini_location+"grid"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " grid file: " << gridfile << endl;	
		
		if (!(param.getValue("feat").isEmpty())) //from command line
		{
			features_file = (string)(param.getValue("feat"));
		}
		else if (!(param.getValue(ini_location+"feat").isEmpty())) //from INI file
		{
			features_file = (string)(param.getValue(ini_location+"feat"));
		}
		if (debug_level>0) log << Date::now() << " " << ini_location << " feature file: " << features_file << endl;	
	
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
				
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
		
		DGridFile grid_file;
		DGrid<2> the_grid;
		grid_file.load(gridfile,the_grid);
		
		DFeatureMapFile fmap_file;
		DFeatureMap<2> feature_map;
		fmap_file.load(features_file,feature_map);
				
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		DMapDewarper<2> map_dewarper;
		map_dewarper.setFeatures(feature_map);
		map_dewarper.setGrid(the_grid);
				
		map_dewarper.dewarp();
		
		DFeatureMap<2> dewarped_features;
		dewarped_features = map_dewarper.getFeatures();
		
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
		fmap_file.store(outfile,dewarped_features);
		
		// For debugging only: write feature in dta2d format
		/*DFeatureMap<2>::const_iterator cit = dewarped_features.begin();
		ofstream featurefile("features1_dewarped.dta2d");
		while (cit != dewarped_features.end() )
		{
			featurefile << cit->getPosition()[0] << " " << cit->getPosition()[1] << " " << cit->getIntensity() << std::endl;
			cit++;
		}
		featurefile.close();*/
		
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




