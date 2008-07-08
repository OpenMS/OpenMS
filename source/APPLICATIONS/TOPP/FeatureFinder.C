// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Ole Schulz-Trieglaff, Marcel Grunert, Clemens Groepl$
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder_impl.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page FeatureFinder FeatureFinder
	
	@brief The feature detection application (quantitation)
	
	This module identifies "features" in a LC/MS map.
	
	By feature, we understand a peptide in a MS sample that
	reveals a characteristic isotope distribution. The algorithm
	computes positions in rt and m/z dimension and a charge estimate
	of each peptide.
	
	The algorithm identifies pronounced regions of raw data points around so-called <tt>seeds</tt>. 
  In the next step, we iteratively fit a model of the isotope profile and the retention time to
  these data points. Data points with a low probability under this model are removed from the
  feature region. The intensity of the feature is then given by the sum of the data points included
  in its regions.
  
  How to find suitable parameters is described in the TOPP tutorial.
  
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureFinder
	: public TOPPBase
{
 public:
	TOPPFeatureFinder()
		: TOPPBase("FeatureFinder","Detects two-dimensional features in LC-MS data.")
	{
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("mzData"));
		registerOutputFile_("out","<file>","","output feature list ");
		setValidFormats_("out",StringList::create("featureXML"));
		registerStringOption_("type","<name>","","FeatureFinder algorithm type\n",true);
		setValidStrings_("type", Factory<FeatureFinderAlgorithm<Peak1D,Feature> >::registeredProducts());
		addEmptyLine_();
		addText_("All other options of the Featurefinder depend on the algorithm type used.\n"
									 "They are set in the 'algorithm' seciton of the INI file.\n");	

		registerSubsection_("algorithm","Algorithm section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		return FeatureFinder().getParameters(type);
	}

	ExitCodes main_(int , const char**)
	{
		//input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");

		Param feafi_param = getParam_().copy("algorithm:",true);

		writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);
				
		String type = getStringOption_("type");
		
		//setup of FeatureFinder
		FeatureFinder ff;
		ff.setLogType(log_type_);
		
		//reading input data
		MSExperiment<Peak1D> exp;
		MzDataFile f;
		f.setLogType(log_type_);
		//prevent loading of fragment spectra
		PeakFileOptions options;
		options.setMSLevels(vector<Int>(1,1));
		f.getOptions() = options;
		f.load(in,exp);

		exp.updateRanges();
		
		//ouput data
		FeatureMap<> features;

		//running algorithm
		ff.run(type, exp, features, feafi_param);

		//-------------------------------------------------------------
		// writing files
		//-------------------------------------------------------------
		FeatureXMLFile map_file;
		map_file.store(out,features);			
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPFeatureFinder tool;
	return tool.main(argc,argv);
}

/// @endcond
