// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

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
  in its regions.
  
  <b>Seeding:</b><br>
  Right now only the seeder with the ID <tt>SimpleSeeder</tt> should be used.
  The critical parameters for this seeder are:
  - <b>min_intensity</b> the minimum intensity to consider a data point as seed
  - <b>intensity_perc</b> Minimum intensity percentage (relative to maximum peak) that a peak has to have
       as a seed. Only used when min_intensity is 0.
  
  <b>Extension:</b><br>
  Right now only the extender with the ID <tt>SimpleExtender</tt> should be used.
	- <b>priority_thr</b> The priority of data point is a function of its intensity and its distance from the seed. 
												Data points with a priority below this threshold are not included into the feature region.
  - <b>dist_mz_up</b> determines how far from the seed data points are searched in direction to higher m/z
  - <b>dist_mz_down</b> determines how far from the seed data points are searched in direction to lower m/z
  - <b>dist_rt_up</b> determines how far from the seed data points are searched in direction to higher RT
  - <b>dist_rt_down</b> determines how far from the seed data points are searched in direction to lower RT
  
  <b>Model fitting:</b><br>
	Right now only the model fitter with the ID <tt>SimpleModelFitter</tt> should be used.
	- <b>min_num_peaks:extended</b> the miniumum number of data points that a extended area has to contain
	- <b>min_num_peaks:final</b> the minimum number of data points that a feature has to contain
	- <b>model_type:first</b> min charge to consider for the model
	- <b>model_type:last</b> max charge to consider for the model
	- <b>quality:type</b> how the quality of a feature is measured. Use 'Correlation'. 
	- <b>quality:minimum</b> min quality that a feature has to achieve
	- <b>isotope_model:stdev:first</b> first std deviation for isotope model peaks to try
	- <b>isotope_model:stdev:last</b> last std deviation for isotope model peaks to try
	- <b>isotope_model:stdev:step</b> steps in between first and last std deviation for isotope model peaks
	- <b>intensity_cutoff_factor</b> intensity ratio (compared to seed)
	- <b>tolerance_stdev_bounding_box</b> influence of the width of the bounding box during the fit
	- <b>intensity_cutoff_factor</b> After fitting, the model is used to truncate the feature region and to remove points with low probability under the model. This is the corresponding threshold.
	- <b>mz:interpolation_step</b> Gives the interpolation step size in m/z domain
	- <b>rt:interpolation_step</b> interpolation step size in time domain
	
	@todo Add test with SimpleModelFitter (Ole)
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinder
	: public TOPPBase
{
 public:
	TOPPFeatureFinder()
		: TOPPBase("FeatureFinder","detects two-dimensional features in LC/MS data")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in MzData format");
		registerStringOption_("out","<file>","","output file in feature format");
		registerIntOption_("buffer_size","<size>",1500,"size of the spectrum buffer used internally", false);
		
		addEmptyLine_();
		addText_("This application implements an algorithm for peptide feature detection\n"
						 "as described in Groepl et al. (2005) Proc. CompLife 05.");
		
		addEmptyLine_();
		addText_("All other options of the Featurefinder depend on the Seeder, Extender and Modelfitter used.\n"
						 "They can be given only in the 'algorithm' seciton  of the INI file.\n"
						 "For a detailled description, please have a look at the doxygen documentation.\n"
						 "How the docu can be built is explained in OpenMS/doc/index.html.");	
		
		registerSubsection_("algorithm");
	}
	
	ExitCodes main_(int , char**)
	{
		//input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");

		FeatureFinder ff;
		Param const& feafi_param = getParam_().copy("algorithm:",true);

		writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);
		
		if (feafi_param.empty())
		{
			writeLog_("No parameters for FeatureFinder modules given. Aborting!");
			return ILLEGAL_PARAMETERS;
		}
		
		ff.setParam(feafi_param);
		
		//New scope => exp is deleted as soon as the FeatureFinder has made a copy
		writeLog_(String("Reading input file ") + in);
		{
			MSExperiment<DPeak<1> > exp;
			MzDataFile().load(in,exp);
			ff.setData(exp.begin(),exp.end(),getIntOption_("buffer_size"));
		}
		writeLog_("Running FeatureFinder...");
		
		DFeatureMap<2> features = ff.run();
	
		//-------------------------------------------------------------
		// writing files
		//-------------------------------------------------------------
	
		writeLog_(String("Writing results to ") + out);
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
