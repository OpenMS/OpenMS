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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzDataFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusPeak.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/StarAlignment.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

typedef DFeatureMap< 2, Feature > FeatureMapType;
typedef DPeakArray< 2, Peak2D > PeakArrayType;
typedef ConsensusFeature< FeatureMapType > ConsensusFeatureType;
typedef ConsensusPeak< PeakArrayType > ConsensusPeakType;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page MapAlignment MapAlignment
 
   @brief Aligns multiple element maps (e.g. feature or peak maps) to one consensus map.
   
   @todo Type (mzData or featureXML) map_type abfangen (try)
   
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignment
      : public TOPPBase
{

  public:
    TOPPMapAlignment()
        : TOPPBase("MapAlignment")
    {}

  protected:
    void printToolUsage_() const
    {
      cerr << endl
      << getToolName() << " -- aligns multiple element (e.g. feature or peak) maps" << endl
      << "Version: " << VersionInfo::getVersion() << endl
      << endl
      << "Usage:" << endl
      << " " << getToolName() << " [options]" << endl
      << endl
      << "Options are:" << endl
      << "  -out <file>       output consensusXML file name" << endl;
    }

    void printToolHelpOpt_() const
    {
      cerr << endl
      << getToolName() << endl
      << endl;
      //         << "INI options:" << endl
      //         << "  optimize_peaks   flag that turns on for the optimization of peak parameters" << endl
      //         << "  in <file>        input mzData file name" << endl
      //         << "  out <file>       output mzData file name" << endl
      //         << endl
      //         << "INI File example section:" << endl
      //         << "  <ITEM name=\"in\" value=\"input.mzData\" type=\"string\"/>" << endl
      //         << "  <ITEM name=\"out\" value=\"output.mzData\" type=\"string\"/>" << endl
      //         << "  <ITEM name=\"optimize_peaks\" value=\"\" type=\"string\"/>" << endl;
    }

    void setOptionsAndFlags_()
    {
      options_["-out"] = "out";
    }

    ExitCodes main_(int , char**)
    {
      //output file name
      String out = getParamAsString_("out");
      
      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      Param mapali_param = getParam_();
      writeDebug_("Parameters:", mapali_param, 2);

      Param files_param = mapali_param.copy("Files:",true);
      writeDebug_("Files parameters:", files_param, 2);
      Param::ConstIterator pit = files_param.begin();
      
      String map_type = mapali_param.getValue("map_type");

      //-------------------------------------------------------------
      // loading input and initialize the alignment object
      //-------------------------------------------------------------
      if (map_type == "feature_map")
      {
        StarAlignment< ConsensusFeatureType > alignment;
        alignment.setParam(mapali_param);
        DFeatureMapFile feature_file;
        std::vector< String > file_names;
        // Vector for the feature maps
        std::vector< FeatureMapType > feature_maps(distance(pit,files_param.end()));

        // Reference to the map vector of the alignment object
        std::vector< FeatureMapType* >& map_vector = alignment.getElementMapVector();
        unsigned int i=0;
        while (pit != files_param.end())
        {
          file_names.push_back(pit->second);
          // load the feature file into a feature_map
          feature_file.load(pit->second, feature_maps[i]);
          map_vector.push_back(&(feature_maps[i]));
          pit++;
          ++i;
        }
        alignment.setFileNames(file_names);
        //-------------------------------------------------------------
        // align
        //-------------------------------------------------------------
        alignment.run();
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
        ConsensusXMLFile cons_file;
        cons_file.store(out,alignment);
      }
      // peak maps
      else
      {
        StarAlignment< ConsensusPeakType > alignment;
        alignment.setParam(mapali_param);
        MzDataFile mzdata_file;
        std::vector< String > file_names;
        // Vector for the feature maps
        std::vector< PeakArrayType > peak_maps(distance(pit,files_param.end()));

        // Reference to the map vector of the alignment object
        std::vector< PeakArrayType* >& map_vector = alignment.getElementMapVector();
        unsigned int i=0;
        while (pit != files_param.end())
        {
          file_names.push_back(pit->second);
          // load the feature file into a feature_maps
          PeakMap ms_exp;
          mzdata_file.load(pit->second, ms_exp);
          ms_exp.get2DData(peak_maps[i]);
          map_vector.push_back(&(peak_maps[i]));
          pit++;
          ++i;
        }
        alignment.setFileNames(file_names);
        //-------------------------------------------------------------
        // align
        //-------------------------------------------------------------
        alignment.run();
        //-------------------------------------------------------------
        // writing output
        //-------------------------------------------------------------
//         ConsensusXMLFile cons_file;
//         cons_file.store(out,alignment);
      }

      return EXECUTION_OK;
    }
};


int main( int argc, char ** argv )
{
  TOPPMapAlignment tool;
  return tool.main(argc,argv);
}

/// @endcond
