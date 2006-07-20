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
// $Maintainer: Ole Schulz-Trieglaff $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/Param.h>

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Date.h>

#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>

#include <OpenMS/KERNEL/DFeature.h>
#include <OpenMS/KERNEL/DFeatureMap.h>
#include <OpenMS/KERNEL/DPeakArray.h>

#include "TOPPBase.h"

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
	@page FeatureFinder FeatureFinder
	
	@brief Executes the feature finding algorithm
	as decribed by Groepl et al. (2005) Proc. CompLife-05.
	
	This module identifies "features" in a LC/MS map.
	By feature, we understand a peptide in a MS sample that
	reveals a characteristic isotope distribution. The algorithm
	computes position in rt and m/z dimension and a charge estimate
	of the peptide.The algorithm identifies pronounced regions of raw data points around so-called <tt>seeds</tt>. 
  In the next step, we iteratively fit a model of the isotope profile and the retention time to
  these data points. Data points with a low probability under this model are removed from the
  feature region. The intensity of the feature is then given by the sum of the data points included
  in its regions.<br>
  There are several parameters that can be used to fine-tune the behaviour of the algorithm. Here we
  give only an overview, for a detailed description of all parameters, please have a look at the
  Doxygen Documentation.
  
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

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES

class TOPPFeatureFinder
	: public TOPPBase
{
	public:
		TOPPFeatureFinder()
			: TOPPBase("FeatureFinder")
		{
			
		}
	
	protected:
		void printToolUsage_()
		{
			 cerr << endl
       		 << tool_name_ << " -- detects two-dimensional features in LC/MS data" << endl
       		 << "This application implements an algorithm for peptide feature detection " << endl
       		 << "as described in Groepl et al. (2005) Proc. CompLife 05" << endl
       		 << endl
       		 << "Usage:" << endl
					 << " " << tool_name_ << " [-in <file>] [-out <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]" << endl
					 << "  -in <file>   input file in mzData format" << endl
					 << "  -out <file>  output file in feature format" << endl
					 << endl;
		}
	
		void printToolHelpOpt_()
		{
			cerr << endl
       		 << tool_name_ << " -- find two-dimensional features in LC/MS data" << endl
       		 << "This application implements an algorithm for peptide feature detection " << endl
       		 << "as described in Groepl et al. (2005) Proc. CompLife 05" << endl
       		 << endl
       		 << "INI options:" << endl
					 << endl
					 << " in    input file" << endl 
					 << " out   output file" << endl 
					 << endl
					 << "All other options of the Featurefinder depend on the Seeder, Extender and Modelfitter used." << endl
					 << "For a detailled description, please have a look at the doxygen documentation." << endl
					 << "How the docu can be built is explained in OpenMS/doc/index.html." << endl
					 << endl ;
		}
	
		void setOptionsAndFlags_()
		{
			//list of all the valid options
			options_["-out"] = "out";
			options_["-in"] = "in";
		}
	
		ExitCodes main_(int , char**)
		{
			//input file names and types
			String in = getParamAsString_("in");
			writeDebug_(String("Input file: ") + in, 1);
			
			String out = getParamAsString_("out");
			writeDebug_(String("Output file: ") + in, 1);
									
			writeLog_(String(" Reading input file ") + in);
			
			MzDataFile mzdata_file;
			MSExperiment<DPeak<1> > exp;
			mzdata_file.load(in,exp);

			String ini_location = String(tool_name_) + ":" + String(instance_number_) + ":";
			
			FeatureFinder ff;
			Param feafi_param = getParamCopy_(ini_location,true);
				
			if (feafi_param.empty())
			{
				writeLog_("No params given, aborting.");
				return ILLEGAL_PARAMETERS;
			}
			
			ff.setParam(feafi_param);
			ff.setData(exp);
	
			writeLog_(" Running FeatureFinder...");
		
			DFeatureMap<2> features = ff.run();
	
			//-------------------------------------------------------------
			// writing files
			//-------------------------------------------------------------
	
			writeLog_(String(" Writing results to ") + out);
			DFeatureMapFile map_file;
			map_file.store(out,features);			
			
			return OK;

		}
};

/// @endcond

int main( int argc, char ** argv )
{
	TOPPFeatureFinder tool;
	return tool.main(argc,argv);
}

