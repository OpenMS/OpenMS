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
// $Id: UnlabeledMatcher.C,v 1.5 2006/05/02 09:39:07 ole_st Exp $
// $Author: ole_st $
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/DSimpleFeatureMatcher.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/DGridFile.h>

//#include <map>
//#include <string>

#include "TOPPBase.h"

using namespace OpenMS;
using namespace std;

typedef DFeature<2,KernelTraits> Feature;
typedef DFeatureMap<2,KernelTraits,Feature> FeatureMap;
typedef DFeatureMapFile FeatureMapFile;
typedef DFeaturePair<2,Feature> FeaturePair;
typedef DFeaturePairVector<2,Feature> FeaturePairVector;
typedef DFeaturePairsFile FeaturePairVectorFile;
typedef DSimpleFeatureMatcher<2,KernelTraits,Feature> FeatureMatcher;
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
    "  " << tool_name_ << " [options]" << endl
    << endl
    << "Options are: " << endl
    << " [-in1 <file>] [-in2 <file>] [-grid <file>] [-pairs <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]\n"
    "  -in1 <file>   input file 1 in xml format\n"
    "  -in2 <file>   input file 2 in xml format\n"
    "  -pairs <file> XML formatted list of feature pairs\n"
    "  -grid <file>  grid covering the feature map\n"
    << endl;
  }

  void printToolHelpOpt_()
  {
    cerr << endl
    << tool_name_ << endl
    << endl
    << "INI options:" << endl
    << "  in1    input file 1 in xml format" << endl
    << "  in2 	 input file 2 in xml format" << endl
    << "  pairs	 XML formatted list of feature pairs)" << endl
    << "  grid   grid covering the feature map" << endl
    << endl
    << "INI File example section:" << endl
    << "  <ITEM name=\"in1\" value=\"input_1.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"in2\" value=\"input_2.mzData\" type=\"string\"/>" << endl
    << "  <ITEM name=\"pairs\" value=\"pairs.xml\" type=\"string\"/>" << endl
    << "  <ITEM name=\"grid\" value=\"grid.xml\" type=\"string\"/>" << endl;
  }

  void setOptionsAndFlags_()
  {
    options_["--help"] = "help";
    options_["-d"] = "debug";
    options_["-in1"] = "in1";
    options_["-in2"] = "in2";
    options_["-ini"] = "ini";
    options_["-log"] = "log";
    options_["-n"] = "instance";
    options_["-grid"] = "grid";
    options_["-pairs"] = "pairs";
    //for debugging the parameters
    options_["unknown"] = "unknown";
    options_["misc"] = "misc";

  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    // input files to be read
    String inputfile[2];

    /// determine names of input files
    for ( Size index = 0; index < 2; ++index )
    {
      string const inputfile_key = string("in") + ('1'+index);
      inputfile[index] = getParamAsString_(inputfile_key);
      writeDebug_(String("Input file: ") + String(index) + ' ' + inputfile_key, 1);
    }

    // determine name ouf grid file
    String gridfilename = getParamAsString_("grid");
    String pairsfile = getParamAsString_("pairs");

    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------

    // read input files
    FeatureMapFile feature_file[2];
    FeatureMap feature_map[2];
    for ( Size index = 0; index < 2; ++index )
    {
      writeLog_(String(" Reading input file ") + String(index+1) + String(", `") + inputfile[index]);
      feature_file[index].load(inputfile[index],feature_map[index]);
    }

    //-------------------------------------------------------------

    // and now ... do the job!

    FeatureMatcher feature_matcher;
    String ini_location = String(tool_name_) + ":" + getParamAsString_("instance") + ":";
    feature_matcher.setParam(getParamCopy_(ini_location));

    for ( Size index = 0; index < 2; ++index )
    {
      feature_matcher.setFeatureMap(index,feature_map[index]);
    }

    FeaturePairVector feature_pair_vector;
    feature_matcher.setFeaturePairs(feature_pair_vector);

    GridType grid;
    feature_matcher.setGrid(grid);

    //log << ini_location << " Running UnlabeledMatcher." << endl;

    feature_matcher.run();

    //-------------------------------------------------------------
    // writing output
    //-------------------------------------------------------------

    writeDebug_(String(" Writing feature pairs, ") + pairsfile + String("\'."),1);
    writeDebug_(String("Number of feature pairs: " + feature_pair_vector.size()),1);

    FeaturePairVectorFile feature_pair_vector_file;
    feature_pair_vector_file.store(pairsfile,feature_pair_vector);

    DGridFile grid_file;
    grid_file.store(gridfilename,feature_matcher.getGrid());

//     DataValue fm_p_d_dfi = feature_matcher.getParam().getValue("debug:dump_feature_input");
//     if ( !fm_p_d_dfi.isEmpty() )
//     {
//       std::string dump_filenameprefix = fm_p_d_dfi;
//       for ( Size index = 0; index < 2; ++index )
//       {
//         std::string dump_filename = dump_filenameprefix+'_'+('0'+index);
//         std::ofstream dump_file(dump_filename.c_str());
//         dump_file << "# " << dump_filename << " generated " << Date::now() << std::endl;
//         dump_file << feature_matcher.getFeatureMap(index) << std::endl;
//         dump_file << "# " << dump_filename << " EOF " << Date::now() << std::endl;
//       }
//     }

    return OK;
  }
};

///@endcond

int main( int argc, char ** argv )
{
  TOPPUnlabeledMatcher tool;
  return tool.main(argc,argv);
}
