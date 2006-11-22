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
#include <OpenMS/KERNEL/MSExperimentExtern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/APPLICATIONS/TOPPBase2.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FeatureFinder FeatureFinder
	
	@brief Executes the feature finding algorithm as decribed by Groepl et al. (2005) Proc. CompLife-05.
	
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

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinder
	: public TOPPBase2
{
 public:
	TOPPFeatureFinder()
		: TOPPBase2("FeatureFinder","detects two-dimensional features in LC/MS data")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in MzData format");
		registerStringOption_("out","<file>","","output file in feature format");
		
		addEmptyLine_();
		addText_("All other options of the Featurefinder depend on the Seeder, Extender and Modelfitter used.\n"
						 "For a detailled description, please have a look at the doxygen documentation.\n"
						 "How the docu can be built is explained in OpenMS/doc/index.html.");
		
		addEmptyLine_();
		addText_("This application implements an algorithm for peptide feature detection\n"
						 "as described in Groepl et al. (2005) Proc. CompLife 05.");
	}
	
	ExitCodes main_(int , char**)
	{
		//input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");
		
		writeLog_(String(" Reading input file ") + in);
			
		MSExperimentExtern<DPeak<1> > exp;
		MzDataFile().load(in,exp);

		FeatureFinder ff;
		Param const& feafi_param = getParam_();

		writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);
		
		if (feafi_param.empty())
		{
			writeLog_("No parameters for FeatureFinder modules given. Aborting!");
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
			
		return EXECUTION_OK;
	}
};


int main( int argc, char ** argv )
{
	TOPPFeatureFinder tool;
	return tool.main(argc,argv);
}

/// @endcond
