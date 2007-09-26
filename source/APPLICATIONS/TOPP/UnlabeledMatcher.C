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
// $Maintainer: Clemens Groepl, Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FeaturePairsXMLFile.h>
#include <OpenMS/FORMAT/GridFile.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <sstream>


using namespace OpenMS;
using namespace std;

typedef ElementPair<Feature> FeaturePair;
typedef std::vector<FeaturePair> FeaturePairVector;
typedef FeaturePairsXMLFile FeaturePairVectorFile;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
  @page UnlabeledMatcher UnlabeledMatcher
  
  @brief For each feature in a given map, this module tries to find its partner in the second map.
  
  This module is the first step in the map matching workflow. It identifies pairs of features in two
  feature map. Currently, two different approaches can be used: if there is only a slight shift
  between the feature positions in the two maps, a simple pairwise matching procedure can be used.
  For more complex situations, an algorithm based on pose clustering can be used to estimate
  a transform and to compute feature pairs based on this transform.
  
  The next steps in the map matching workflow are done by MapMatcher and MapDewarper.
*/

// We do not want this class to show up in the docu, thus:
/// @cond TOPPCLASSES

class TOPPUnlabeledMatcher
	: public TOPPBase
{
  public:
    TOPPUnlabeledMatcher()
        : TOPPBase("UnlabeledMatcher","matches common two-dimensional features/peaks of two LC/MS maps")
    {
    }

  protected:

    void registerOptionsAndFlags_()
    {
			registerStringOption_("in1","<file>","","input FeatureXML file 1");
			registerStringOption_("in2","<file>","","input FeatureXML file 2");
			registerStringOption_("pairs","<file>","","output file in FeaturePairsXML format");
			registerStringOption_("grid","<file>","","output file: grid covering the feature map");

			addEmptyLine_();
			addText_("All other options can be given only in the 'algorithm' section  of the INI file.\n");
    	registerSubsection_("algorithm","Clustering parameters section");
    }

		Param getSubsectionDefaults_(const String& /*section*/) const
		{
			return PoseClusteringPairwiseMapMatcher<>().getDefaults();
		}

    ExitCodes main_(int , const char**)
    {
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------

      // determine name of grid file
      String gridfilename = getStringOption_("grid");

      // determine name of pairs file
      String pairsfile = getStringOption_("pairs");

      // input files to be read
      String inputfile[2];

      // read input files
      FeatureXMLFile feature_file[2];
      FeatureMap< > feature_map[2];

      // determine names of input files
      for ( UInt index = 0; index < 2; ++index )
      {
        inputfile[index] = getStringOption_(String("in") + (index+1));
        writeLog_(String("Reading input file ") + (index+1) + ", `" + inputfile[index] + '\'');
        feature_file[index].load(inputfile[index],feature_map[index]);
      }
			
      //-------------------------------------------------------------


      // the resulting feature pairs go here
      FeaturePairVector feature_pair_vector;

      PoseClusteringPairwiseMapMatcher<> poseclust_feature_matcher;
			
			Param param_alg = getParam_().copy("algorithm:",true);
			writeDebug_("Parameters passed to PoseClusteringMapMatcher", param_alg,3);

      poseclust_feature_matcher.setParameters(param_alg);

      for ( UInt index = 0; index < 2; ++index )
      {
        poseclust_feature_matcher.setElementMap(index,feature_map[index]);
      }

      writeDebug_("Running algorithm.",1);

      poseclust_feature_matcher.run();

      writeDebug_("Running algorithm...done.",1);

      writeDebug_(String("Number of feature pairs: ") + String(poseclust_feature_matcher.getElementPairs().size()),1);
      writeDebug_(String("Writing feature pairs file `") + pairsfile + "'.",1);

      FeaturePairVectorFile feature_pair_vector_file;
      feature_pair_vector_file.store(pairsfile,poseclust_feature_matcher.getElementPairs());

      writeDebug_(String("Writing grid file `") + gridfilename + "'.",1);

      GridFile grid_file;
      grid_file.store(gridfilename,poseclust_feature_matcher.getGrid());

      writeDebug_("Running UnlabeledMatcher...done.",1);

      return EXECUTION_OK;
    }

};


int main( int argc, const char** argv )
{
  TOPPUnlabeledMatcher tool;
  return tool.main(argc,argv);
}

/// @endcond

