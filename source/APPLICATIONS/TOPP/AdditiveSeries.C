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
// $Id: AdditiveSeries.C,v 1.5 2006/03/28 08:03:30 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Ole Schulz-Trieglaff$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/Param.h>
#include<OpenMS/FORMAT/DFeatureMapFile.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include<OpenMS/KERNEL/DFeatureMap.h>
#include<OpenMS/KERNEL/DFeature.h>
#include<OpenMS/KERNEL/DPosition.h>
#include<OpenMS/KERNEL/DimensionDescription.h>

#include <TOPPCommon.h>

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include<vector>
#include<algorithm>

#include <gsl/gsl_math.h>

#include <OpenMS/MATH/STATISTICS/LinearRegression.h>

using namespace OpenMS;
using namespace std;

typedef DFeature<2>::CoordinateType CoordinateType;
typedef DFeature<2>::IntensityType IntensityType;

/// Defines the coordinates of peaks / features.
enum DimensionId
{
	RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
	MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
};

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page AdditiveSeries AdditiveSeries
	
	@brief Computes an additive serives to quantify a peptide in a set of samples.
	
	This module computes an additve series for an absolute
	quantification of a peptide in a set of samples. The
	output consisits of a GNUplot script which can be used
	to visualise the results and the results in XML format.
	
	In this version, the application computes the additive
	series as a ratio of the intensities of two different peptides.
	One of these peptides serves as internal standard for
	calibration. For details of the procedure, please have
	a look at the publications:
	
	Groepl at al. (2005) Proc. CompLife pages 151-163 
	and
	Mayr et al. (2006) Journal of Proteome Research (5), pp. 414-421
	
	There are several parameters that influence the behaviour of
	this application. This is an overview of the most important
	ones, for a full description please have a look at the INI file
	in the subdirectory Examples.
	<ul>
	<li><b>write_gnuplot_output</b>:
	If set to true, a file with GNUplot commands is written that
	draw the regression line together with its error bars.
	</li>
	<li><b>mz_tolerance</b>:
	m/z range in which we search for the feature.
	</li>
	<li><b>rt_tolerance</b>:
		m/z range in which we search for the feature.
	</li>
	</ul>
		
	@ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "AdditiveSeries";


// searches for a features with coordinates within the tolerance in this map
bool readMapFile(String filename, vector<double>& intensities, 
                 CoordinateType tol_mz, CoordinateType tol_rt,
                 DPosition<2> fpos1, DPosition<2> fpos2)
{
	DFeatureMapFile map_file;
	DFeatureMap<2> map;
	map_file.load(filename,map);
	
	DFeature<2>* feat1 = 0;
	DFeature<2>* feat2 = 0;
	
	DFeatureMap<2>::iterator iter = map.begin();
	while (iter!= map.end() )
	{
				
		if ( (iter->getPosition()[RT] <  fpos1[RT] + tol_rt) && 
		     (iter->getPosition()[RT] >  fpos1[RT] - tol_rt) &&
		     (iter->getPosition()[MZ] <  fpos1[MZ] + tol_mz) &&
		     (iter->getPosition()[MZ] >  fpos1[MZ] - tol_mz) )
		{
			// feature at correct position found, save intensity
			if (!feat1) feat1 = &(*iter);
					
		}
		
		if ( (iter->getPosition()[RT] <  fpos2[RT] + tol_rt) && 
		     (iter->getPosition()[RT] >  fpos2[RT] - tol_rt) &&
		     (iter->getPosition()[MZ] <  fpos2[MZ] + tol_mz) &&
		     (iter->getPosition()[MZ] >  fpos2[MZ] - tol_mz) )
		{
			// same as above
			if (!feat2) feat2 = &(*iter);
		}
		
		iter++;
	}	// end of while
	
	if (feat1 != 0 && feat2 != 0)
	{
		intensities.push_back( feat1->getIntensity() / feat2->getIntensity());
	
		return true;
	}
	
	return false;
}

/**\brief Computes the linear regression for a series of measurements, the
	 x-axis intercept of the regression line and its confidence interval, and
	 writes a couple of files from which a nice plot of all this can be
	 generated using the gnuplot program.
 */
bool computeRegressionAndWriteGnuplotFiles ( std::vector<double>::const_iterator const conc_vec_begin,
																						 std::vector<double>::const_iterator const conc_vec_end,
																						 std::vector<double>::const_iterator const area_vec_begin,
																						 double const confidence_p,
																						 String const filename_prefix,
																						 String const output_filename,
																						 String const format = "",
																						 bool const write_gnuplot = true
																					 )
{
	
  try 
	{ 
		LinearRegression<vector<double>::const_iterator> linreg;

		linreg.computeInterceptXAxis ( confidence_p, conc_vec_begin, conc_vec_end, area_vec_begin );
		
		if (write_gnuplot)
		{
		
			// the peak data goes here
			String datafilename(filename_prefix);
			datafilename+=String(".dat");
			ofstream dataout(datafilename.c_str());
	
			// the gnuplot commands go here
			String commandfilename(filename_prefix);
			commandfilename+=String(".cmd");
			ofstream cmdout(commandfilename.c_str());
		
			// the error bar for the x-axis intercept goes here
			String errorbarfilename(filename_prefix);
			errorbarfilename+=String(".err");
			ofstream errout(errorbarfilename.c_str());

			// writing the commands
			cmdout <<
				"set ylabel \"ion count\"\n"
				"set xlabel \"concentration\"\n"
				"set key left Left reverse\n"
				"set title \"" << filename_prefix << "\"\n"
				;
			
			if ( ! format.empty() )
			{
				if ( format == "png" )
				{
					cmdout <<
						"set terminal png \n"
						"set output \"" << filename_prefix << ".png\"\n" ;
				}
				else if ( format == "eps" )
				{
					cmdout <<
						"set terminal postscript eps \n"
						"set output \"" << filename_prefix << ".eps\"\n" ;
				}
			
			}
			cmdout <<
				"plot \""  << datafilename <<"\"  w points ps 2 pt 1 lt 8 title \"data\" " // want data on first line of key
				",  " << linreg.getIntercept() << "+" <<  linreg.getSlope() << "*x lt 2 lw 3 title \"linear regression: "
						  << linreg.getIntercept() << " + " <<  linreg.getSlope() << " * x\" "
				", \""  << datafilename <<"\"  w points ps 2 pt 1 lt 8 notitle " // draw data a second time, on top of reg. line
				", \"" << errorbarfilename << "\"  using ($1):(0) w points pt 13 ps 2 lt 1 title \"x-intercept: " << linreg.getXIntercept() << "\" "
				", \"" << errorbarfilename << "\"  w xerrorbars lw 3 lt 1 title \"95% interval: [ " << linreg.getLower() << ", " << linreg.getUpper() << " ]\"\n"
				;
			cmdout.close();
	  
			// writing the x-axis intercept error bar
			errout << linreg.getXIntercept() << " 0 " << linreg.getLower() << " " << linreg.getUpper() << endl;
			errout.close();
	  
			// writing the peak data points
			vector<double>::const_iterator cit = conc_vec_begin;
			vector<double>::const_iterator ait = area_vec_begin;
			dataout.precision(15);
			for ( ;cit != conc_vec_end; ++cit, ++ait )
			{
				dataout << *cit << ' ' << *ait << '\n';
			}
			dataout.close();
		
		} // end if (write_gnuplot)
		
		// write results to XML file
		ofstream results;
		results.open(output_filename.c_str());
		
		results << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>" << endl;
		results << "<results_additiveseries>" << endl; 
		results << "\t<x_intercept>" << linreg.getXIntercept() << "</x_intercept>" << endl;
		results << "\t<confidence_lowerlimit>" << linreg.getLower() << "</confidence_lowerlimit>" << endl;
		results << "\t<confidence_upperlimit>" << linreg.getUpper() << "</confidence_upperlimit>" << endl;
		results << "</results_additiveseries>" << endl; 	

		results.close();
	}
  catch (std::string s)
	{
		cout << s <<  endl;
		return 1;
	} 

  return 0;
}


void print_usage()
{
	cerr << endl
       << tool_name << " -- Computes an additive series to quantify" << endl 
       << "a peptide in a set of samples." << endl
       << "Detailed procedure is described in Groepl et. (2005) Proc. CompLife-05 ." << endl
       << endl
       << "Usage:" << endl
			 << " " << tool_name << " [options]" << endl
			 << endl
			 << "Options are:" << endl
			 << "  -in <file>   input file containing the spiked concentrations (default read from INI file)" << endl
			 << "  -out <file>  output file XML containg regression line and confidence interval (default read from INI file)" << endl
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
       << tool_name << " -- Computes an additive series to quantify a peptide in a set of samples." << endl
       << "Detailed procedure is described in Groepl et. (2005) Proc. CompLife-05 ." << endl
       << endl
       << "INI options (excerpt):" << endl
			 << endl
			 << " -write_gnuplot_output : True => write script with GNUplot commands" << endl
			 << " -mz_tolerance         : m/z range for feature coordinates" << endl
			 << " -rt_tolerance         : rt range for feature coordinates" << endl
			 << endl
			 << "All further parameters are explained in the example INI file" << endl
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
				
		Param add_param = param.copy(ini_location,true);
		
		DPosition<2> feat_pos1;
		DPosition<2> feat_pos2;
		
		if (add_param.getValue("mz_tolerance").isEmpty() || add_param.getValue("rt_tolerance").isEmpty() )
		{
			log << Date::now() << " " << ini_location << "Tolerances not set in INI file. Aborting."<< endl;
			return ILLEGAL_PARAMETERS;
		}

		CoordinateType tol_mz = (CoordinateType) add_param.getValue("mz_tolerance");
		CoordinateType tol_rt = (CoordinateType) add_param.getValue("rt_tolerance");
		
		if (add_param.getValue("in").isEmpty() || add_param.getValue("out").isEmpty() )
		{
			log << Date::now() << " " << ini_location << "Input / outputfile not given. Aborting."<< endl;
			return ILLEGAL_PARAMETERS;
		}
	
		String conc_f = (String) add_param.getValue("in");
		String out_f  = (String) add_param.getValue("out");
		
		if (add_param.getValue("Feature:MZ").isEmpty() || add_param.getValue("Feature:RT").isEmpty() )
		{
			log << Date::now() << " " << ini_location << "Feature coordinates not given. Aborting."<< endl;
			return ILLEGAL_PARAMETERS;
		}		
		feat_pos1[MZ] = (CoordinateType) add_param.getValue("Feature:MZ");
		feat_pos1[RT] = (CoordinateType) add_param.getValue("Feature:RT");
	
		if (add_param.getValue("Standard:MZ").isEmpty() || add_param.getValue("Standard:RT").isEmpty() )
		{
			log << Date::now() << " " << ini_location << "Standard coordinates not given. Aborting."<< endl;
			return ILLEGAL_PARAMETERS;
		}		
	
		feat_pos2[MZ] = (CoordinateType) add_param.getValue("Standard:MZ");
		feat_pos2[RT] = (CoordinateType) add_param.getValue("Standard:RT");
	
		log << Date::now() << " " << ini_location << " Setting tolerances to " << tol_mz << " " << tol_rt << endl;
		log << Date::now() << " " << ini_location << " Feature position 1: " << feat_pos1 << endl;
		log << Date::now() << " " << ini_location << " Feature position 2: " << feat_pos2 << endl;
	
		String title = (String) add_param.getValue("title");
	
		// read the spiked concentrations 		
		ifstream conc_file(conc_f.c_str());
		char line[256];
		vector<double> sp_concentrations;
		while (conc_file.getline(line,256))
		{
			String line_str(line);
			line_str.trim();
			sp_concentrations.push_back(line_str.toDouble());
		}
	
		// introduce a flag for each concetration
		// true => the corresponding feature was found
		vector<bool> flags;
	
		// fetching list of files
		vector<String> files;
		Param file_param = add_param.copy("Files:",true);
	
		Param::const_iterator pit = file_param.begin();
	
		while (pit != file_param.end() )
		{
			files.push_back(pit->second);
			pit++;
		}	

		sort(files.begin(),files.end());
	
		// collect features
		vector<IntensityType> intensities;
		vector<String>::const_iterator cit = files.begin();
		while (cit != files.end())
		{
			cout << "Opening file " << *cit << endl;
		
			if (readMapFile(*cit,intensities,tol_mz,tol_rt,feat_pos1,feat_pos2) )
			{
		 		flags.push_back(true);
			}
			else 
			{
		 		flags.push_back(false);
			}
			cit++;
		}
	
		vector<double> sp_concentrations2;
		for (unsigned int i=0; i<sp_concentrations.size(); i++)
		{
			if (flags.at(i) == true )
			{
			sp_concentrations2.push_back( sp_concentrations.at(i) );
			}
		}
	
	
		if (intensities.size() == 0 || sp_concentrations.size() == 0 )
		{
			log << Date::now() << " " << ini_location << " Did not find any data! Aborting..." << endl;
			return UNKNOWN_ERROR;		
		}
		
		// set prefix of gnuplot output
		String filename_prefix = "gnuplot_";
		filename_prefix += String(title);
		
		DataValue dv = add_param.getValue("write_gnuplot_output");
		if (!dv.isEmpty() && dv.toString() != "false") 
		{
			// compute regression and write GNUplot files
			computeRegressionAndWriteGnuplotFiles(sp_concentrations2.begin(), sp_concentrations2.end(),
																						intensities.begin(), 0.95, filename_prefix, out_f, "eps", true);
		}
		else
		{
			log << Date::now() << " " << ini_location << " No GNUplot output is written..." << endl;	
			// compute regression and write GNUplot files
			computeRegressionAndWriteGnuplotFiles(sp_concentrations2.begin(), sp_concentrations2.end(),
																						intensities.begin(), 0.95, filename_prefix, out_f, "eps", false);
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




