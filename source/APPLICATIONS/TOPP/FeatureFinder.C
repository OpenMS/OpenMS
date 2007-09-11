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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinder.h>
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
	computes position in rt and m/z dimension and a charge estimate
	of the peptide.The algorithm identifies pronounced regions of raw data points around so-called <tt>seeds</tt>. 
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
		: TOPPBase("FeatureFinder","detects two-dimensional features in LC/MS data")
	{
			
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerStringOption_("in","<file>","","input file in MzData format");
		registerStringOption_("out","<file>","","output file in FeatureXML format");
		registerStringOption_("type","<name>","","FeatureFinder algorithm type ('simple', )");
		
		addEmptyLine_();
		addText_("This application implements an algorithm for peptide feature detection\n"
						 "as described in Groepl et al. (2005) Proc. CompLife 05.");
		
		addEmptyLine_();
		addText_("All other options of the Featurefinder depend on the Seeder, Extender and Modelfitter used.\n"
						 "They can be given only in the 'algorithm' seciton  of the INI file.\n");	

		registerSubsection_("algorithm","Algorithm section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		Param tmp;
		
		FeatureFinder ff;
		try
		{
			tmp.insert("",ff.getParameters(getStringOption_("type")));
		}
		catch(Exception::RequiredParameterNotGiven)
		{
			cout << "Error: Required parameter 'type' not given!" << endl;
			tmp.setValue("algorithm:dummy","value","Here the algorithms of the FeatureFinder are given!",true);
		}
		return tmp;
	}

	ExitCodes main_(int , char**)
	{
		//input file names and types
		String in = getStringOption_("in");	
		String out = getStringOption_("out");

		Param feafi_param = getParam_().copy("algorithm:",true);

		writeDebug_("Parameters passed to FeatureFinder", feafi_param, 3);
				
		String type = getStringOption_("type");
		if (type!="simple")
		{
			writeLog_("Invalid FeatureFinder type given. Aborting!");
			return ILLEGAL_PARAMETERS;
		}
		
		//setup of FeatureFinder
		FeatureFinder ff;
		ff.setLogType(log_type_);
		
		//reading input data
		writeLog_(String("Reading input file ") + in);
		MSExperiment<RawDataPoint1D> exp;
		MzDataFile f;
		f.setLogType(log_type_);
		f.load(in,exp);
		exp.updateRanges();
		
		//ouput data
		FeatureMap<> features;

		//running algorithm
		writeLog_("Running FeatureFinder...");
		
		ff.run(type, exp, features, feafi_param);

		//-------------------------------------------------------------
		// writing files
		//-------------------------------------------------------------
	
		writeLog_(String("Writing results to ") + out);
		FeatureXMLFile map_file;
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
