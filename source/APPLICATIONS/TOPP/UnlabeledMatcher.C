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
// $Maintainer: Clemens Groepl, Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DGeomHashPairwiseMapMatcher.h> // the new one
#include <OpenMS/ANALYSIS/MAPMATCHING/DSimpleFeatureMatcher.h> // the old one
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/DGridFile.h>

#include "TOPPBase.h"

#include <sstream>


using namespace OpenMS;
using namespace std;

typedef DFeature<2,KernelTraits> Feature;
typedef DFeatureMap<2,Feature> FeatureMap;
typedef DFeatureMapFile FeatureMapFile;
typedef DFeaturePair<2,Feature> FeaturePair;
typedef DFeaturePairVector<2,Feature> FeaturePairVector;
typedef DFeaturePairsFile FeaturePairVectorFile;
typedef DSimpleFeatureMatcher<2,KernelTraits,Feature> SimpleFeatureMatcherType;
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

// We do not want this class to show up in the docu -> @cond
/// @cond TOPPCLASSES

class TOPPUnlabeledMatcher
      : public TOPPBase
{
 public:
  TOPPUnlabeledMatcher()
		: TOPPBase("UnlabeledMatcher")
  {}

 protected:
  void printToolUsage_()
  {
    cerr << endl
				 << tool_name_ << " -- match common two-dimensional features of two LC/MS data sets\n"
			"\n"
			"Usage:\n"
			"  " << tool_name_ << " [-in1 <file>] [-in2 <file>] [-grid <file>] [-pairs <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]\n\n"
			"Options are:\n"
			"  -in1 <file>   input feature file 1\n"
			"  -in2 <file>   input feature file 2\n"
			"  -pairs <file> output file: XML formatted list of feature pairs\n"
			"  -grid <file>  output file: grid covering the feature map\n"
				 << endl;
  }

  void printToolHelpOpt_()
  {
    cerr << "\n"
				 << tool_name_ << "\n"
			"\n"
			"INI options:\n"
			"  in1    input feature file 1\n"
			"  in2 	  input feature file 2\n"
			"  pairs  output file: XML formatted list of feature pairs\n"
			"  grid   output file: grid covering the feature map\n"
			"\n"
			"INI File example section:\n"
			"  <ITEM name=\"in1\" value=\"input_1.mzData\" type=\"string\"/>\n"
			"  <ITEM name=\"in2\" value=\"input_2.mzData\" type=\"string\"/>\n"
			"  <ITEM name=\"pairs\" value=\"pairs.xml\" type=\"string\"/>\n"
			"  <ITEM name=\"grid\" value=\"grid.xml\" type=\"string\"/>\n"
			"Note: many more parameters can be set in the INI File.\n"
			;
  }

  void setOptionsAndFlags_()
  {
    options_["-in1"] = "in1";
    options_["-in2"] = "in2";
    options_["-grid"] = "grid";
    options_["-pairs"] = "pairs";
  }

  ExitCodes main_(int , char**)
  {

		writeDebug_("--------------------------------------------------",1);
		writeDebug_("Running UnlabeledMatcher.",1);

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

		String param_path = tool_name_ + ':' + String(instance_number_) + ':';
		Param param = getParamCopy_(param_path,true);
		{
			std::string bla;
			std::stringstream blastream(bla);
			blastream << param;
			writeDebug_(String("Parameters from `")+param_path+"' are:\n"
									+blastream.str(),2);
		}
		const int inheritance_steps_max = 15;
		int inheritance_steps = 0;
		for (	String param_inherit_path = param.getValue("inherit");
					!param_inherit_path.empty();
					param_inherit_path = param.getValue("inherit")
				)
		{
			if ( ++inheritance_steps > inheritance_steps_max )
			{
				writeLog_("Too many inheritance steps ("+String(inheritance_steps_max)+" allowed).  Perhaps there is a cycle?");
				return ILLEGAL_PARAMETERS;
			}
			param.remove("inherit");
			param.setDefaults(getParamCopy_(param_inherit_path,true),"",false);
			{
				std::string bla;
				std::stringstream blastream(bla);
				blastream << param;
				writeDebug_(String("Parameters after inheriting from `")+param_inherit_path+"' are:\n"
										+blastream.str(),2);
			}
		}

    // input files to be read
    String inputfile[2];

    /// determine names of input files
    for ( Size index = 0; index < 2; ++index )
    {
      string const inputfile_key = string("in") + ('1'+index);
      inputfile[index] = param.getValue(inputfile_key);
      writeDebug_(String("Input file: ") + String(index) + ' ' + inputfile_key, 1);
    }

    // determine name ouf grid file
    String gridfilename = param.getValue("grid");
    String pairsfile = param.getValue("pairs");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // read input files
    FeatureMapFile feature_file[2];
    FeatureMap feature_map[2];
    for ( Size index = 0; index < 2; ++index )
    {
      writeLog_(String("Reading input file ") + String(index+1) + String(", `") + inputfile[index]+'\'');
      feature_file[index].load(inputfile[index],feature_map[index]);
    }

    //-------------------------------------------------------------

		// polymorphism, the hard way...
		enum Algorithm
			{
				ALGORITHM_NONE = 0,
				ALGORITHM_SIMPLE = 1,
				ALGORITHM_GEOMHASH_SHIFT = 2,
				// ...
				ALGORITHM_OLD = 99
			} algorithm = ALGORITHM_NONE;

		if      ( param.getValue("algorithm") == "simple" )
		{
			writeDebug_("Using algorithm \"simple\".",1);
			algorithm = ALGORITHM_SIMPLE;
		}
		else if ( param.getValue("algorithm") == "geomhash_shift" )
		{
			writeDebug_("Using algorithm \"geomhash_shift\".",1);
			algorithm = ALGORITHM_GEOMHASH_SHIFT;
		}
		else if ( param.getValue("algorithm") == "old" )
		{
			writeDebug_("Using algorithm \"old\".",1);
			algorithm = ALGORITHM_OLD;
		}
		else
		{
			writeLog_("Error: No existing algorithm specified!  Stop.");
			// We do not provide a default, so just return.
			return ILLEGAL_PARAMETERS;
		}
		
		writeDebug_(String("(Useless fact: This algorithm has number ")+String(algorithm)+" in class TOPPUnlabeledMatcher (see UnlabeledMatcher.C))", 13);
			
    //-------------------------------------------------------------


		// the resulting feature pairs go here
		FeaturePairVector feature_pair_vector;


		// the resulting grid goes here
		GridType grid;


		switch ( algorithm )
		{
			// --------------------------------------------------
			// We can still distinguish among them by the "algorithm" parameter setting.
		case ALGORITHM_SIMPLE:
		case ALGORITHM_GEOMHASH_SHIFT:
			{

				DGeomHashPairwiseMapMatcher<2> geomhash_feature_matcher;
				
				geomhash_feature_matcher.setParam(param);

				for ( Size index = 0; index < 2; ++index )
				{
					geomhash_feature_matcher.setFeatureMap(index,feature_map[index]);
				}

				geomhash_feature_matcher.setFeaturePairs(feature_pair_vector);

				writeDebug_("Running algorithm.",1);

				geomhash_feature_matcher.run();
				
				writeDebug_("Running algorithm...done.",1);

				// this does a deep copy
				grid = geomhash_feature_matcher.getGrid();

				// this does a deep copy
				feature_pair_vector = geomhash_feature_matcher.getFeaturePairs();

				break;
			}
			// --------------------------------------------------
		case ALGORITHM_OLD:
			{
				// This is the old algorithm.  It shall be replaced soon.

				SimpleFeatureMatcherType simple_feature_matcher;

				simple_feature_matcher.setParam(param);

				for ( Size index = 0; index < 2; ++index )
				{
					simple_feature_matcher.setFeatureMap(index,feature_map[index]);
				}

				simple_feature_matcher.setFeaturePairs(feature_pair_vector);

				simple_feature_matcher.setGrid(grid);

				writeDebug_("Running.",1);

				simple_feature_matcher.run();

				writeDebug_("Running...done.",1);

				break;
			}
			// --------------------------------------------------
		default:
			{
				writeLog_("Unrecognized algorithm number in switch statement - I am confused!");
				return UNKNOWN_ERROR;
				break;
			}
		} // switch (algorithm)

		writeDebug_(String("Number of feature pairs: ") + String(feature_pair_vector.size()),1);
		writeDebug_(String("Writing feature pairs file `") + pairsfile + String("'."),1);
	
		FeaturePairVectorFile feature_pair_vector_file;
		feature_pair_vector_file.store(pairsfile,feature_pair_vector);
	
		writeDebug_(String("Writing grid file `") + gridfilename + String("'."),1);
	
		DGridFile grid_file;
		grid_file.store(gridfilename,grid);

		writeDebug_("Running UnlabeledMatcher...done.",1);
 
		return OK;
	
	} // main_()

}; // TOPPUnlabeledMatcher


int main( int argc, char ** argv )
{
  TOPPUnlabeledMatcher tool;
  return tool.main(argc,argv);
}

// @endcond

