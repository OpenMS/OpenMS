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
// $Id: UnlabeledMatcher.C,v 1.5 2006/05/02 09:39:07 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DSimpleFeatureMatcher.h>

#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/DGridFile.h>

#include <fstream>
#include <iostream>
#include <map>
#include <string>

#include <TOPPCommon.h>

using namespace OpenMS;
using namespace std;

typedef DFeature<2,KernelTraits> Feature;
typedef DFeatureMap<2,KernelTraits,Feature> FeatureMap;
typedef DFeatureMapFile FeatureMapFile;
typedef DFeaturePair<2,Feature> FeaturePair;
typedef DFeaturePairVector<2,Feature> FeaturePairVector;
typedef DFeaturePairsFile FeaturePairVectorFile;
typedef DSimpleFeatureMatcher<2,KernelTraits,Feature> FeatureMatcher;
typedef DGrid<2> GridType;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UnlabeledMatcher UnlabeledMatcher
	
	@brief For each feature in a given map, this
	module tries to find its partner in the second map.
	
	This module is the first step in the map matching
	workflow. It identifies pairs of features in two
	feature map. Currently, two different approaches
	can be used: if there is only a slight shift
	between the feature positions in the two maps,
	a simple pairwise matching procedure can be used.
	For more complex situations, an algorithm based
	on geometric hashing can be used to estimate
	a transform and to compute feature pairs based
	on this transform.	
		
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "UnlabeledMatcher";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- match common two-dimensional features of two LC/MS data sets\n" 
		"\n"
		"Usage:\n"
		"  " << tool_name << 
		" [-in1 <file>] [-in2 <file>] [-grid <file>] [-pairs <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]\n" 
		"  -in1 <file>   input file 1 in xml format (default read from INI file)\n"
		"  -in2 <file>   input file 2 in xml format (default read from INI file)\n"
		"  -pairs <file> XML formatted list of feature pairs (default read from INI file)\n"
		"  -grid <file>  grid covering the feature map (default read from INI file)\n"
		" Common TOPP options are:\n" 
		"  -ini <file>  TOPP INI file (default: TOPP.ini)\n"
		"  -log <file>  log file (default: TOPP.log)\n"
		"  -n <int>     instance number (default: 1)\n"
		"  -d <level>   sets debug level (default: 0)\n"
		"  --help       shows this help\n"
		;
}

//-------------------------------------------------------------
// main program
//-------------------------------------------------------------

int main( int argc, char ** argv )
{
	// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
	String ini_location = tool_name;
	// path to the log file
	String logfile = "";
	// debug level
	int debug_level = 0;
	// log filestream
	ofstream log;
	log.open("TOPP.log", ofstream::out | ofstream::app);
	
	// input files to be read
	String inputfile[2];
	
	// output files to be written
	String gridfilename  = "";
	String pairsfile = "";
	
	//-------------------------------------------------------------
	// command line parsing
	//-------------------------------------------------------------
	
	//list of all the valid options
	map<string,string> valid_options;
	valid_options["--help"] = "help";
	valid_options["-d"] = "debug";
	valid_options["-in1"] = "in1";
	valid_options["-in2"] = "in2";
	valid_options["-ini"] = "ini";
	valid_options["-log"] = "log";
	valid_options["-n"] = "instance";
	valid_options["-grid"] = "grid";
	valid_options["-pairs"] = "pairs";
	//for debugging the parameters
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
		
		/// determine names of input files
		for ( Size index = 0; index < 2; ++index )
		{
			string const inputfile_key = string("in") + ('1'+index);
			if (!(param.getValue(inputfile_key).isEmpty()))
			{
				inputfile[index] = (string)(param.getValue(inputfile_key));
			}
			else if (!(param.getValue(ini_location+inputfile_key).isEmpty()))
			{
				inputfile[index] = (string)(param.getValue(ini_location+inputfile_key));
			}
			else 
			{
				cout <<  " Could not find input file " << index+1 << ". Aborting!" << endl;
				return INPUT_FILE_NOT_FOUND;
			}
		}

		//-------------------------------------------------------------
		
		// determine name ouf grid file
		if (!(param.getValue("grid").isEmpty()))
		{
			gridfilename = (string)(param.getValue("grid"));
		}
		else if (!(param.getValue(ini_location+"grid").isEmpty()))
		{
			gridfilename = (string)(param.getValue(ini_location+"grid"));
		}
		else 
		{
			log << ini_location << " No file name for grid file given. Aborting!" << endl;
			return CANNOT_WRITE_OUTPUT_FILE;
		}
		
		// determine name ouf feature pairs file
		if (!(param.getValue("pairs").isEmpty()))
		{
			pairsfile = (string)(param.getValue("pairs"));
		}
		else if (!(param.getValue(ini_location+"pairs").isEmpty()))
		{
			pairsfile = (string)(param.getValue(ini_location+"pairs"));
		}
		else 
		{
			log << ini_location << " No file name for pairs file given. Aborting!" << endl;
			return CANNOT_WRITE_OUTPUT_FILE;
		}

		//-------------------------------------------------------------

		// read input files
		FeatureMapFile feature_file[2];
		FeatureMap feature_map[2];
		for ( Size index = 0; index < 2; ++index )
		{
			log << ini_location << " Reading input file " << index+1 << ", `" << inputfile[index] << "\'." << std::endl;
			feature_file[index].load(inputfile[index],feature_map[index]);
		}

		//-------------------------------------------------------------

		// and now ... do the job!

		FeatureMatcher feature_matcher;

		feature_matcher.setParam(param.copy(ini_location,true));

		for ( Size index = 0; index < 2; ++index )
		{
			feature_matcher.setFeatureMap(index,feature_map[index]);
		}

		FeaturePairVector feature_pair_vector;
		feature_matcher.setFeaturePairs(feature_pair_vector);

		GridType grid;
		feature_matcher.setGrid(grid);

		log << ini_location << " Running UnlabeledMatcher." << endl;
	
		feature_matcher.run();

		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
	
		log << ini_location << " Writing feature pairs, `" << pairsfile << "\'." << std::endl;
		log << ini_location << feature_pair_vector.size() << endl;
		
		FeaturePairVectorFile feature_pair_vector_file;
		feature_pair_vector_file.store(pairsfile,feature_pair_vector);

		DGridFile grid_file;
		grid_file.store(gridfilename,feature_matcher.getGrid());

		DataValue fm_p_d_dfi = feature_matcher.getParam().getValue("debug:dump_feature_input");
		if ( !fm_p_d_dfi.isEmpty() )
		{
			std::string dump_filenameprefix = fm_p_d_dfi;
			for ( Size index = 0; index < 2; ++index )
			{
				std::string dump_filename = dump_filenameprefix+'_'+('0'+index);
				std::ofstream dump_file(dump_filename.c_str());
				dump_file << "# " << dump_filename << " generated " << Date::now() << std::endl;
				dump_file << feature_matcher.getFeatureMap(index) << std::endl;
				dump_file << "# " << dump_filename << " EOF " << Date::now() << std::endl;
			}
		}

		//-------------------------------------------------------------
	
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
