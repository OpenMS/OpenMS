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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/DFeaturePairsFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/PairMatcher.h>
#include <OpenMS/KERNEL/ComparatorUtils.h>
#include <OpenMS/KERNEL/DFeatureMap.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

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
 
	This module identifies pairs of isotope-labeled features in a LC/MS features map.
	By feature, we understand a peptide in a MS sample that reveals a characteristic isotope distribution.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPLabeledMatcher
      : public TOPPBase
{
	public:
		
	  TOPPLabeledMatcher()
	    : TOPPBase("LabeledMatcher","find pairs of labeled features in LC/MS data")
	  {
	  }

	protected:
	
	  void registerOptionsAndFlags_()
	  {
			registerStringOption_("in","<file>","","input file in FeatureMap format");
			registerStringOption_("out","<file>","","output file in FeaturePairs format");
			registerStringOption_("best","<file>","","output file of the best pairs in FeaturePairs format", false);
	  	addEmptyLine_();
	  	addText_("RT and m/z shifts and ranges can currently only be given in the 'algorithm' part of INI file:\n"
							 "  <NODE name=\"algorithm\">\n"
							 "    <ITEM name=\"rt_pair_dist\" value=\"0.5\" type=\"float\"/>\n"
							 "    <ITEM name=\"rt_stdev_low\" value=\"0.22\" type=\"float\"/>\n"
							 "    <ITEM name=\"rt_stdev_high\" value=\"0.65\" type=\"float\"/>\n"
							 "    <ITEM name=\"mz_pair_dist\" value=\"4.0\" type=\"float\"/>\n"
							 "    <ITEM name=\"mz_stdev\" value=\"0.025\" type=\"float\"/>\n" 
							 "  </NODE>");
			addEmptyLine_();
			addText_("Note: The mz_pair_dist is added while the rt_pair_dist is substracted when searching pairs.\n"
			         "      This is due to the fact, that the heavier peptide normally elutes earlier!");
			registerSubsection_("algorithm");
	  }
	
	  ExitCodes main_(int , char**)
	  {
	
	    //-------------------------------------------------------------
	    // parameter handling
	    //-------------------------------------------------------------
	
	    // input file to be read
	    String inputfile = "";
	
	    // determine name of input file
	    inputfile = getStringOption_("in");
	    String outputfile = getStringOption_("out");
	
	    // determine name ouf visualization output file
	    String best_outputfile = getStringOption_("best");
	
	    //-------------------------------------------------------------
	    // reading input
	    //-------------------------------------------------------------
	
	
	    DFeatureMap<2> features;
	    DFeatureMapFile().load(inputfile,features);
	
	    // sort input file
	    enum DimensionId
	    {
	      RT = DimensionDescription < LCMS_Tag >::RT,
	      MZ = DimensionDescription < LCMS_Tag >::MZ
	    };
	
	    typedef DFeature<2>::NthPositionLess< RT > RTless;
	    typedef DFeature<2>::NthPositionLess<MZ> MZless;
	    sort(features.begin(),features.end(), LexicographicComparator<RTless,MZless>());
	
	    PairMatcher pm(features);
	
	    Param pm_param = getParam_().copy("algorithm:",true);
	    writeDebug_("Parameters passed to PairMatcher", pm_param, 3);
	    pm.setParameters(pm_param);
	
	    writeDebug_(" Running LabeledMatcher...",1);
	
	    const DFeaturePairVector<2>* pairs = &pm.run();
	
	    //-------------------------------------------------------------
	    // writing files
	    //-------------------------------------------------------------
	
	    writeDebug_(String(" Writing results to ") + outputfile, 1 );
	    DFeaturePairsFile().store(outputfile,*pairs);
			
			writeDebug_(String(" Writing results to ") + best_outputfile, 1 );
			if (best_outputfile!="")
	    {
	    	DFeaturePairsFile().store(best_outputfile,pm.getBestPairs());
	    }
			
	    return EXECUTION_OK;
	  }
};


int main( int argc, char ** argv )
{
  TOPPLabeledMatcher tool;
  return tool.main(argc,argv);
}

/// @endcond
