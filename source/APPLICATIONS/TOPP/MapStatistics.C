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
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/DFeatureMap.h>

#include <TOPPCommon.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page MapStatistics MapStatistics
	
	@brief Computes a five-number summary of
	intensities in raw data, picked peak or feature map.
		
	This TOPP module computes a five-number summary
	of the feature intensities and qualities in a map.
	
	The five-number summary consists of median, upper
	and lower quartile, minimum and maximum. These values
	are computed for qualities and intensities. They 
	give a measure of spread and location and are stored
	in a XML format for further processing.
	
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "MapStatistics";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
	cerr << endl
       << tool_name << " -- Computes a five-number summary for " << endl
       << " features / raw data intensities and qualities in a map." << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>        feature or raw data map (default read from INI file)" << endl
			 << "  -in_type <file>   either feat or mzData (default read from INI file)" << endl
			 << "  -out <file>  output file in XML format (default read from INI file)" << endl
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
	valid_options["-in"] = "in";
	valid_options["-ini"] = "ini";
	valid_options["-log"] = "log";
	valid_options["-n"] = "instance";
	valid_options["-d"] = "debug";
	valid_options["--help"] = "help";
	valid_options["-in_type"] = "in_type";
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
		// reading input
		//-------------------------------------------------------------
		String inputfile, outputfile;
		String in_type; // type of inputfile (either .mzData for raw data maps or .feat for feature maps) 
		
		//determine input file name and type
		if (!(param.getValue("in").isEmpty())) //from command line
		{
			inputfile = (string)(param.getValue("in"));
	
			//determine input file type
			if (!(param.getValue("in_type").isEmpty()))
			{
				in_type = (string)(param.getValue("in_type"));
			}
			if (debug_level>1) log << Date::now() << " " << ini_location << " input file was determined from command line!" << endl;
		}
		else if (!(param.getValue(ini_location+"in").isEmpty())) //from INI file
		{
			inputfile = (string)(param.getValue(ini_location+"in"));
			
			//determine input file type
			if (!(param.getValue(ini_location+"in_type").isEmpty()))
			{
				in_type = (string)(param.getValue(ini_location+"in_type"));
			}
			if (debug_level>1) log << Date::now() << " " << ini_location << " input file was determined from INI file!" << endl;
		}
		else 
		{
			log << Date::now() << " " << ini_location << " Could not find input file. Aborting!" << endl;
			return INPUT_FILE_NOT_FOUND;
		}
	
		if (in_type=="")
		{
			in_type = inputfile.suffix('.');
			if (debug_level>1) log << Date::now() << " " << ini_location << " input file type is determined from file extension!" << endl;
		}	
		in_type.toUpper();
				
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
						
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
		
		// We need to distinguish whether we deal with peak
		// or feature data because each type requires slightly
		// different statistics
		
		if (in_type == "FEAT")
		{
			DFeatureMap<2> map;
			DFeatureMapFile map_file;
			map_file.load(inputfile,map);		
			
			unsigned int size = map.size();
			
			typedef DFeatureMap<2>::FeatureType::IntensityType IntensityType;
			typedef DFeatureMap<2>::FeatureType::QualityType QualityType;
		
			IntensityType * intensities = new IntensityType[ size ];
			QualityType * qualities     = new QualityType[ size ];
				
			for (unsigned int i = 0; i < size; 	++i)
			{
				intensities[i] = map.at(i).getIntensity();
				qualities[i]   = map.at(i).getOverallQuality();
			}
			
			gsl_sort(intensities, 1, size);
			gsl_sort(qualities, 1, size);
			
			double mean_int, var_int, max_int, min_int;
			mean_int = gsl_stats_mean(intensities,1,size);
  		var_int  = gsl_stats_variance(intensities,1,size);
  		max_int  = gsl_stats_max(intensities,1,size);
  		min_int  = gsl_stats_min(intensities,1,size);
  		
  		double mean_q, var_q, max_q, min_q;
			mean_q = gsl_stats_mean(qualities,1,size);
  		var_q  = gsl_stats_variance(qualities,1,size);
  		max_q  = gsl_stats_max(qualities,1,size);
  		min_q  = gsl_stats_min(qualities,1,size);
			
			double median_int, upperq_int, lowerq_int;
			median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
	 		upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
  		lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
  		
  		double median_q, upperq_q, lowerq_q;
			median_q = gsl_stats_median_from_sorted_data(qualities,1,size);
	 		upperq_q = gsl_stats_quantile_from_sorted_data(qualities,1,size,0.75);
  		lowerq_q = gsl_stats_quantile_from_sorted_data (qualities,1,size,0.25);
			
			delete [] intensities;
			delete [] qualities;
			
			ofstream out(outputfile.c_str());
			out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
			out << "<mapstatistics>" << endl;
		
			out << "\t<intensities>" << endl;
			out << "\t\t<mean>" << mean_int << "</mean>" << endl;
			out << "\t\t<median>" << median_int << "</median>" << endl;
			out << "\t\t<variance>" << var_int << "</variance>" << endl;
			out << "\t\t<min>" << min_int << "</min>" << endl;
			out << "\t\t<max>" << max_int << "</max>" << endl;
			out << "\t\t<lower_quartile>" << lowerq_int << "</lower_quartile>" << endl;
			out << "\t\t<upper_quartile>" << upperq_int << "</upper_quartile>" << endl;
			out << "\t</intensities>" << endl;
		
			out << "\t<qualities>" << endl;
			out << "\t\t<mean>" << mean_q << "</mean>" << endl;
			out << "\t\t<median>" << median_q << "</median>" << endl;
			out << "\t\t<variance>" << var_q << "</variance>" << endl;
			out << "\t\t<min>" << min_q << "</min>" << endl;
			out << "\t\t<max>" << max_q << "</max>" << endl;
			out << "\t\t<lower_quartile>" << lowerq_q << "</lower_quartile>" << endl;
			out << "\t\t<upper_quartile>" << upperq_q << "</upper_quartile>" << endl;
			out << "\t</qualities>" << endl;
		
			out << "</mapstatistics>" << endl;
		
			out.close(); 
						
		}
		else if (in_type == "MZDATA")
		{
			MSExperiment< DPeak<1> > exp;
			MzDataFile f;
			f.load(inputfile,exp);		
			
			DPeakArray<2> array;
			exp.get2DData(array);
									
			unsigned int size = array.size();
			
			typedef DPeak<1>::IntensityType IntensityType;
					
			IntensityType * intensities = new IntensityType[ size ];
							
			for (unsigned int i = 0; i < size; 	++i)
			{
				intensities[i] = array[i].getIntensity();
			}
			
			gsl_sort(intensities, 1, size);
						
			double mean_int, var_int, max_int, min_int;
			mean_int = gsl_stats_mean(intensities,1,size);
  		var_int  = gsl_stats_variance(intensities,1,size);
  		max_int  = gsl_stats_max(intensities,1,size);
  		min_int  = gsl_stats_min(intensities,1,size);
  		  		
			double median_int, upperq_int, lowerq_int;
			median_int = gsl_stats_median_from_sorted_data(intensities,1,size);
	 		upperq_int = gsl_stats_quantile_from_sorted_data(intensities,1,size,0.75);
  		lowerq_int = gsl_stats_quantile_from_sorted_data (intensities,1,size,0.25);
  				
			delete [] intensities;
			
			ofstream out(outputfile.c_str());
			out << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
			out << "<mapstatistics>" << endl;
		
			out << "\t<intensities>" << endl;
			out << "\t\t<mean>" << mean_int << "</mean>" << endl;
			out << "\t\t<median>" << median_int << "</median>" << endl;
			out << "\t\t<variance>" << var_int << "</variance>" << endl;
			out << "\t\t<min>" << min_int << "</min>" << endl;
			out << "\t\t<max>" << max_int << "</max>" << endl;
			out << "\t\t<lower_quartile>" << lowerq_int << "</lower_quartile>" << endl;
			out << "\t\t<upper_quartile>" << upperq_int << "</upper_quartile>" << endl;
			out << "\t</intensities>" << endl;
		
			out << "</mapstatistics>" << endl;
		
			out.close(); 
					
		}
		else
		{
			log << Date::now() << " " << ini_location << " Unknown file type '" << in_type << "' given. Aborting!" << endl;
			cout << "Unknown file type '" << in_type << "' given. Aborting!" << endl;
			print_usage();
			return ILLEGAL_PARAMETERS;			
		}
				
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




