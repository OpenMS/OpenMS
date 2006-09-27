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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/DFeatureMap.h>

#include "TOPPBase.h"

#include <iostream>
#include <fstream>
#include <string>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page LabeledMatcher LabeledMatcher
 
	@brief Executes the pair matching algorithm for labeled peptides.
 
	This module identifies pairs of labeled "features" in a LC/MS map.
	By feature, we understand a peptide in a MS sample that
	reveals a characteristic isotope distribution.
 
  
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


class TOPPLabeledMatcher
      : public TOPPBase
{
public:
  TOPPLabeledMatcher()
      : TOPPBase("LabeledMatcher")
  {}

protected:
  void printToolUsage_()
  {
    cerr << endl
    << tool_name_ << " -- find pairs of labeled features in LC/MS data" << endl
    << endl
    << "Usage:" << endl
    << " " << tool_name_ << " [-in <file>] [-out <file>] [-ini <file>] [-log <file>] [-n <int>] [-d <level>]" << endl
    << "  -in  <file>       input file" << endl
    << "  -out <file>  	    output file" << endl
    << "  -vis_all <file>   output file of all pairs for visualisation in TOPPView"
    << "  -vis_best <file>  output file of the best pairs for visualisation in TOPPView"
    << ""
    << "The boundaries on m/z and RT dimension can be defined only in the INI file!"
    << endl;
  }

   void printToolHelpOpt_()
  {
    cerr << endl
    << tool_name_ << endl
    << endl
    << "INI options:" << endl
    << "in        input file in mzData format (default read from INI file)" << endl
    << "out  	  output file" << endl
    << "vis_all   output file of all pairs "
    << "vis_best  output file of the best pairs "
    << endl
    << "INI File example section:" << endl
    << "  <ITEM name=\"in\" value=\"input.xml\" type=\"string\"/>" << endl
    << "  <ITEM name=\"out\" value=\"output.xml\" type=\"string\"/>" << endl
    << "  <ITEM name=\"vis_all\" value=\"output_all_pairs.xml\" type=\"string\"/>" << endl
    << "  <ITEM name=\"vis_best\" value=\"output_vis_best.xml\" type=\"string\"/>" << endl
		<< "  <NODE name=\"algorithm\">" << endl
		<< "    <ITEM name=\"rt_pair_dist\" value=\"0.5\" type=\"float\"/>" << endl
		<< "    <ITEM name=\"rt_stdev_low\" value=\"0.22\" type=\"float\"/>" << endl
		<< "    <ITEM name=\"rt_stdev_high\" value=\"0.65\" type=\"float\"/>" << endl
		<< "    <ITEM name=\"mz_pair_dist\" value=\"4.0\" type=\"float\"/>" << endl
		<< "    <ITEM name=\"mz_stdev\" value=\"0.025\" type=\"float\"/>" << endl
		<< "  </NODE>" << endl;
  }

  void setOptionsAndFlags_()
  {
    //list of all the valid options
    options_["-out"] = "out";
    options_["-in"] = "in";
    options_["-vis_best"] = "vis_best";
    options_["-vis_all"] = "vis_all";
  }

  ExitCodes main_(int , char**)
  {

    //-------------------------------------------------------------
    // parameter handling
    //-------------------------------------------------------------

    // input file to be read
    String inputfile = "";

    // output file to be written
    String outputfile = "";
    String vis_all_outputfile = "";
    String vis_best_outputfile = "";

    // determine name of input file
    inputfile = getParamAsString_("in");
    outputfile = getParamAsString_("out");

    // determine name ouf visualization output file
    vis_all_outputfile = getParamAsString_("vis_all");
    vis_best_outputfile = getParamAsString_("vis_best");



    //-------------------------------------------------------------
    // reading input
    //-------------------------------------------------------------


    DFeatureMap<2> features;
    DFeatureMapFile().load(inputfile,features);

    // sort input file
    enum DimensionId
    {
      RT = DimensionDescription < DimensionDescriptionTagLCMS >::RT,
      MZ = DimensionDescription < DimensionDescriptionTagLCMS >::MZ
    };

    typedef DFeature<2>::NthPositionLess< RT > RTless;
    typedef DFeature<2>::NthPositionLess<MZ> MZless;
    std::sort(features.begin(),features.end(), LexicographicComparator<RTless,MZless>());

    PairMatcher pm(features);

    String ini_location = String(tool_name_) + ":" + String(instance_number_) + ":";
    pm.setParam(getParamCopy_(ini_location+"algorithm:"));

    //pm.setDebugLevel(debug_level);
    //pm.setDebugStream(&log);
    //pm.setInstanceId(ini_location);

    writeLog_(" Running LabeledMatcher...");

    const DFeaturePairVector<2>* pairs = &pm.run();

    // save pairs in DFeatureMap for visualization in TOPPView
    // (until visualization of DFeaturePairFile is available)
    if (vis_all_outputfile!="")
    {
      DFeatureMap<2> map;
      PairMatcher::fillFeatureMap(map,*pairs);
      DFeatureMapFile().store(vis_all_outputfile,map);
    }

    //     log << Date::now() << " " << ini_location << "\nAll pairs:\n";
    //     PairMatcher::printInfo(log,*pairs);

    if (vis_best_outputfile!="")
    {
      DFeatureMap<2> map;
      pairs = &pm.getBestPairs();
      PairMatcher::fillFeatureMap(map,*pairs);
      DFeatureMapFile().store(vis_best_outputfile,map);
    }

    //     log << Date::now() << " " << ini_location << "\nBest pairs:\n";
    //     PairMatcher::printInfo(log,*pairs);

    //-------------------------------------------------------------
    // writing files
    //-------------------------------------------------------------

    writeLog_(String(" Writing results to ") + outputfile);
    DFeaturePairsFile().store(outputfile,*pairs);

    return OK;
  }
};


int main( int argc, char ** argv )
{
  TOPPLabeledMatcher tool;
  return tool.main(argc,argv);
}

/// @endcond
