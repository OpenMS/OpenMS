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
// $Id: RawDataResampler.C,v 1.3 2006/05/30 15:46:40 marc_sturm Exp $
// $Author: marc_sturm $
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/Param.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <TOPPCommon.h>

#include <OpenMS/KERNEL/MSExperiment.h>

# include <OpenMS/SYSTEM/StopWatch.h>

#include <OpenMS/FILTERING/TRANSFORMERS/LinearResampler.h>

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
   @page RawDataResampler RawDataResampler

    The RawDataResampler module can be used to generate uniform data from non-uniform raw data (e.g. ESI-TOF or MALDI-TOF experiments).
    Therefore the intensity at every position x in the input raw data is spread to the two
    adjacent resampling points.

    This method preserves the area of the input signal and also the centroid position of a peak.
    Therefore it is recommended for quantitation as well as for identification experiments.

  @note Use this method only for high resoluted data (< 0.1 Th between two adjacent raw data points).
       The resampling rate should be >= the accuracy.

   @ingroup TOPP
*/

//-------------------------------------------------------------
// command line name of this tool
//-------------------------------------------------------------
const char* tool_name = "RawDataResampler";

//-------------------------------------------------------------
// description of the usage of this TOPP tool
//-------------------------------------------------------------

void print_usage()
{
  cerr << endl
  << tool_name << " -- generate equally spaced raw data" << endl
  << "This application implements a linear resampling method" << endl
  << "which preserves the total area of the input data" << endl
  << "as well as the peak's centroids. (The default sampling rate is 0.05Th.)"
  << endl
  << endl
  << "Use this module only for high resoluted data " << endl
  << "(< 0.1 Th between two adjacent raw data points)." << endl
  << endl
  << "Usage:" << endl
  << " " << tool_name << " [-in <file>] [-out <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]" << endl
  << "  -spacing <resampling_rate> spacing of the equally spaced resampled data (default read from INI file)" << endl
  << "  -in <file>   input file in MzData format (default read from INI file)" << endl
  << "  -out <file>  output file in MzData format (default read from INI file)" << endl
  << endl
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
  String ini_location = "RawDataResampler";
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
  valid_options["-spacing"] = "spacing";
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
    else
    {
      logfile = "TOPP.log";
    }
    if (debug_level>0) log << Date::now() << " " << ini_location << " log file: " << logfile << endl;
    log.close();
    log.open(logfile.c_str(), ofstream::out | ofstream::app);

    //-------------------------------------------------------------
    // parsing parameters
    //-------------------------------------------------------------
    //File names and types
    String inputfile, outputfile;
    double spacing = 0.;

    // determine filter type
    if (!(param.getValue("spacing").isEmpty()))//from command line
    {
      spacing = (double)(param.getValue("spacing"));
    }
    else if (!(param.getValue(ini_location+"ResamplingWidth").isEmpty())) //from INI file
    {
      spacing = (double)(param.getValue(ini_location+"ResamplingWidth"));
    }

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

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------
    cout << "SPACING " << spacing << endl;

    StopWatch timer;
    timer.reset();
    timer.start();

    MzDataFile mz_data_file;
    MSExperiment<DRawDataPoint<1> > ms_exp_raw;
    mz_data_file.load(inputfile,ms_exp_raw);

    timer.stop();
    cout << "read end " << timer.getUserTime() << endl;

    if (debug_level>0) log << Date::now() << " " << ini_location << " Number of spectra in input file: " << ms_exp_raw.size() << endl;

    //-------------------------------------------------------------
    // calculations
    //-------------------------------------------------------------
    Param th_param = param.copy(ini_location,true);

    timer.reset();
    timer.start();
    LinearResampler< DRawDataPoint <1> > linear_resampler(th_param);
    linear_resampler.setSpacing(spacing);

    MSExperiment<DRawDataPoint<1> > ms_exp_resampled;
    ms_exp_raw >> linear_resampler(ms_exp_resampled);

    timer.stop();
    cout << "resampling end " << timer.getUserTime() << endl;


    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    timer.reset();
    timer.start();

    if (debug_level>0) log << Date::now() << " " << ini_location << " Number of spectra for writing: " << ms_exp_resampled.size() << endl;

    mz_data_file.store(outputfile,ms_exp_resampled);

    timer.stop();
    cout << "write end " << timer.getUserTime()  << endl;
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




