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
// $Id: PairMatcher.C,v 1.4 2006/06/09 08:19:41 j-joachim Exp $
// $Author: j-joachim $
// $Maintainer: Jens Joachim $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>

#include <OpenMS/KERNEL/DFeatureMap.h>

#include <TOPPCommon.h>

#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page PairMatcher PairMatcher

	@brief Executes the pair matching algorithm
	for labeled peptides.

	This module identifies pairs of labeled "features" in a LC/MS map.
	By feature, we understand a peptide in a MS sample that
	reveals a characteristic isotope distribution.

  
  <ul>
	<li><b>min_intensity</b> :
	minimum intensity of a seed
	</li>
	<li><b>priority_thr</b> :
	The priority of data point is a function of its intensity and
	its distance from the seed. Data points with a priority below this
	threshold are not included into the feature region.</li>
	<li><b>min_quality</b> :
	minimum quality of feature, if smaller feature will be discarded</li>
	<li><b>intensity_cutoff_factor</b>
	For each data points in the feature region, we compute its probability given
	the model. Data points with a probability below this cutoff are discarded.</li>
  </ul>

	@ingroup TOPP
*/



//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "PairMatcher";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- find pairs of labeled features in LC/MS data" << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [-in <file>] [-out <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]" << endl
			 << "  -in <file>   input file in mzData format (default read from INI file)" << endl
			 << "  -out <file>  output file (default read from INI file)" << endl
			 << "  -vis_all <file>  output file of all pairs "
			 << "for visualisation in TOPPView (default read from INI file)" << endl
			 << "  -vis_best <file>  output file of the best pairs "
			 << "for visualisation in TOPPView (default read from INI file)" << endl
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

void print_helpopts()
{
	cerr << endl
       << tool_name << " -- find pairs of labeled features in LC/MS data" << endl
       << endl
       << "INI options:" << endl
			 << endl
			 << " min_rt : minimum difference in retention time of second peptide to first" << endl
			 << " max_rt : maximum difference in retention time of second peptide to first" << endl
			 << " max_mz : maximum deviation from optimal m/z-difference "
			 << "between the features of a pairs" << endl
			 << endl
			 << "For a detailled description, please have a look at the doxygen documentation." << endl
			 << endl ;
}

//-------------------------------------------------------------
// main program
//-------------------------------------------------------------

int main( int argc, char ** argv )
{
	// instance specific location of settings in INI file (e.g. 'TOPP_Skeleton:1:')
	String ini_location = "PairMatcher";
	// path to the log file
	String logfile = "";
	// debug level
	int debug_level = 0;
	// log filestream
	ofstream log;
	log.open("TOPP.log", ofstream::out | ofstream::app);
	
	// input file to be read
	String inputfile = "";

	// output file to be written
	String outputfile = "";
	String vis_all_outputfile = "";
	String vis_best_outputfile = "";

	//-------------------------------------------------------------
	// command line parsing
	//-------------------------------------------------------------
	
	//list of all the valid options
	map<string,string> valid_options;
	valid_options["-out"] = "out";
	valid_options["-in"] = "in";
	valid_options["-vis_best"] = "vis_best";
	valid_options["-vis_all"] = "vis_all";
	valid_options["-ini"] = "ini";
	valid_options["-log"] = "log";
	valid_options["-n"] = "instance";
	valid_options["-d"] = "debug";
	valid_options["--help"] = "help";
	valid_options["--help-opt"] = "helpopt";
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

	//-------------------------------------------------------------
	// check command line options
	//-------------------------------------------------------------
	
	// '--help' given
	if (!(param.getValue("help").isEmpty()))
	{
		print_usage();
		return OK;
	}

	// '--help-opt' given
	if (!(param.getValue("helpopt").isEmpty()))
	{
		print_helpopts();
		return OK;
	}

	// test if unknown options were given
	if (!param.getValue("unknown").isEmpty())
	{
		log << Date::now() << " " << ini_location << " Unknown option '" << (string)(param.getValue("unknown")) << "' given. Aborting!" << endl;
		print_usage();
		return ILLEGAL_PARAMETERS;
	}
	
	// test if unknown text argument were given (we do not use them)
	if (!param.getValue("misc").isEmpty())
	{
		log << Date::now() << " " << ini_location << " Trailing text argument '" << (string)(param.getValue("misc")) << "' given. Aborting!" << endl;
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
	else
	{
		logfile = "TOPP.log";
	}
	if (debug_level>0) cout << ini_location << " log file: " << logfile << endl;
	log.close();
	log.open(logfile.c_str(), ofstream::out | ofstream::app);
	
	//-------------------------------------------------------------
	// calculations
	//-------------------------------------------------------------

	// determine name of input file
	if (!(param.getValue("in").isEmpty()))
	{
		inputfile = (string)(param.getValue("in"));
	}
	else if (!(param.getValue(ini_location+"in").isEmpty()))
	{
		inputfile = (string)(param.getValue(ini_location+"in"));
	}
	else 
	{
		log << Date::now() << " " << ini_location << " Could not find input file. Aborting!" << endl;
		return INPUT_FILE_NOT_FOUND;
	}

	// determine name ouf output file
	if (!(param.getValue("out").isEmpty()))
	{
		outputfile = (string)(param.getValue("out"));
	}
	else if (!(param.getValue(ini_location+"out").isEmpty()))
	{
		outputfile = (string)(param.getValue(ini_location+"out"));
	}
	else
	{
		log << Date::now() << " " << ini_location << " No output file given. Aborting!" << endl;
		return CANNOT_WRITE_OUTPUT_FILE;
	}


	// determine name ouf visualization output file
	if (!(param.getValue("vis_all").isEmpty()))
	{
		vis_all_outputfile = (string)(param.getValue("vis_all"));
	}
	else if (!(param.getValue(ini_location+"vis_all").isEmpty()))
	{
		vis_all_outputfile = (string)(param.getValue(ini_location+"vis_all"));
	}
	else
	{
		vis_all_outputfile = "";
	}

	// determine name ouf visualization output file
	if (!(param.getValue("vis_best").isEmpty()))
	{
		vis_best_outputfile = (string)(param.getValue("vis_best"));
	}
	else if (!(param.getValue(ini_location+"vis_best").isEmpty()))
	{
		vis_best_outputfile = (string)(param.getValue(ini_location+"vis_best"));
	}
	else
	{
		vis_best_outputfile = "";
	}


	log << Date::now() << " " << ini_location << " Reading input file " << inputfile << endl;

	DFeatureMap<2> features;
	DFeatureMapFile().load(inputfile,features);

	// sort input file
	enum DimensionId
	{
			RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
			MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
	};
	typedef DFeature<2>::NthPositionLess< RT > RTless;
	typedef DFeature<2>::NthPositionLess<MZ> MZless;
	std::sort(features.begin(),features.end(), LexicographicComparator<RTless,MZless>());

	PairMatcher pm(features);
	pm.setParam(param.copy(ini_location+"algorithm:",true));

	//pm.setDebugLevel(debug_level);
	//pm.setDebugStream(&log);
	//pm.setInstanceId(ini_location);

	log << Date::now() << " " << ini_location << " Running PairMatcher..." << endl;


	const DFeaturePairVector<2>* pairs = &pm.run();

	// save pairs in DFeatureMap for visualization in TOPPView
	// (until visualization of DFeaturePairFile is available)
	if (vis_all_outputfile!="")
	{
		DFeatureMap<2> map;
		PairMatcher::fillFeatureMap(map,*pairs);
		DFeatureMapFile().store(vis_all_outputfile,map);
	}

	log << Date::now() << " " << ini_location << "\nAll pairs:\n";
	PairMatcher::printInfo(log,*pairs);

	if (vis_best_outputfile!="")
	{
		DFeatureMap<2> map;
		pairs = &pm.getBestPairs();
		PairMatcher::fillFeatureMap(map,*pairs);
		DFeatureMapFile().store(vis_best_outputfile,map);
	}

	log << Date::now() << " " << ini_location << "\nBest pairs:\n";
	PairMatcher::printInfo(log,*pairs);

	//-------------------------------------------------------------
	// writing files
	//-------------------------------------------------------------
	log << Date::now() << " " << ini_location << " Writing results to " << outputfile << endl;
	DFeaturePairsFile().store(outputfile,*pairs);


	}
	catch(Exception::UnableToCreateFile& e)
	{
		log << Date::now() << " " << ini_location << " Error: Unable to write file (" << e.what() <<")"<< endl;
		return CANNOT_WRITE_OUTPUT_FILE;
	}
	catch(Exception::FileNotFound& e)
	{
		log << Date::now() << " " << ini_location << " Error: File not found (" << e.what() <<")"<< endl;
		return INPUT_FILE_NOT_FOUND;
	}
	catch(Exception::ParseError& e)
	{
		log << Date::now() << " " << ini_location << " Error: Unable to read file (" << e.what() <<")"<< endl;
		return INPUT_FILE_CORRUPT;
	}
	catch(Exception::InvalidValue& e)
	{
		log << Date::now() << " " << ini_location << " Error: Invalid parameter value in ini-file (" << e.what() << ")" << endl;
		return INPUT_FILE_CORRUPT;
	}
	catch(Exception::Base& e)
	{
		log << Date::now() << " " << ini_location << " Error: Unknown error (" << e.what() <<")"<< endl;
		return UNKNOWN_ERROR;
	}
	log.close();

	return OK;
}
