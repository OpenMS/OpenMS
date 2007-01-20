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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/MzDataFile.h>

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/ConsensusPeak.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/DFeatureMapFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

typedef DFeatureMap< 2, Feature > FeatureMapType;
typedef DPeakArray< 2, Peak2D > PeakArrayType;
typedef ConsensusFeature< FeatureMapType > ConsensusFeatureType;
typedef ConsensusPeak< PeakArrayType > ConsensusPeakType;
typedef ConsensusMap< ConsensusFeatureType > ConsensusMapType;


//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
   @page MapAlignment MapAlignment
 
   @brief Aligns multiple element maps to one consensus map.
   
   This application implements an algorithm for the alignment of mulitple maps.
   It accepts feature maps (in featureXML), peak maps (in mzData) or consensus maps (in ConsensusXML).
   This tool requires an INI file with at least the names of the input files and the map_type.
   Parameters for the alignment algorithm can be given only in the 'algorithm' seciton  of the INI file.
   
   @Note If you use consensus maps, the consensus elements are used as normal elements and you will
         loose the former consensus information.
      
   @ingroup TOPP
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPMapAlignment
      : public TOPPBase
{

  public:
    TOPPMapAlignment()
        : TOPPBase("MapAlignment","aligns multiple feature, peak or consensus maps")
    {
    }

  protected: 
    void registerOptionsAndFlags_()
    {
      registerStringOption_("out","<file>","","output consensusXML file name");
      
      addEmptyLine_();
      addText_("This application implements an algorithm for the alignment of mulitple maps.\n"
               "It accepts feature maps (in featureXML), peak maps (in mzData) or consensus maps (in ConsensusXML)\n"
               "Note: If you use consensus maps , the consensus elements are used as normal elements and you will\n"
               "loose the former consensus information.");

      addEmptyLine_();
      addText_("This tool requires an INI file with at least the names of the input files and the map_type.\n"
              	"Parameters for the alignment algorithm can be given only in the 'algorithm' seciton  of the INI file:\n"
								"  <NODE name=\"file_names\">\n"
								"    <ITEM name=\"1\" value=\"file1.xml\" type=\"string\"/>\n"
								"    <ITEM name=\"2\" value=\"file2.xml\" type=\"string\"/>\n"
								"    <ITEM name=\"3\" value=\"file3.xml\" type=\"string\"/>\n"
								"  </NODE>\n"
								"  <NODE name=\"algorithm\">\n"
								"    <ITEM name=\"map_type\" value=\"feature_map\" type=\"string\"/>\n"
								"    ...\n"
								"  </NODE>");
			
      registerSubsection_("algorithm");
      registerSubsection_("file_names");
    }


    ExitCodes main_(int , char**)
    {
      //output file name
      String out = getStringOption_("out");

      //-------------------------------------------------------------
      // parameter handling
      //-------------------------------------------------------------
      Param const& mapali_param = getParam_().copy("algorithm:",true);
      writeDebug_("Parameters:", mapali_param, 2);
      
      Param files_param = getParam_().copy("file_names:",true);
      writeDebug_("Files parameters:", files_param, 2);
      Param::ConstIterator pit = files_param.begin();

      String map_type = getParam_().getValue("algorithm:map_type");
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
          try
          {
            feature_file.load(pit->second, feature_maps[i]);
          }
          catch(Exception::FileNotFound& e)
          {
            writeLog_(String("File not found '") + (String)pit->second + "'. Aborting!");
            return INPUT_FILE_NOT_FOUND;
          }
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
        if (map_type == "peak_map")
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

            try
            {
              mzdata_file.load(pit->second, ms_exp);
            }
            catch(Exception::FileNotFound& e)
            {
              writeLog_(String("File not found '") + (String)pit->second + "'. Aborting!");
              return INPUT_FILE_NOT_FOUND;
            }
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
          ConsensusXMLFile cons_file;
          cons_file.store(out,alignment);
        }
        else if (map_type == "consensus_map")
        {
          StarAlignment< ConsensusFeature< ConsensusMapType > > alignment;
          alignment.setParam(mapali_param);

          ConsensusXMLFile cons_file;
          std::vector< String > file_names;
          // Vector for the feature maps
          std::vector< ConsensusMapType > cons_maps(distance(pit,files_param.end()));

          // Reference to the map vector of the alignment object
          std::vector< ConsensusMapType* >& map_vector = alignment.getElementMapVector();
          unsigned int i=0;
          while (pit != files_param.end())
          {
            file_names.push_back(pit->second);
            // load the feature file into a feature_map
            try
            {
              cons_file.load(pit->second, cons_maps[i]);
            }
            catch(Exception::FileNotFound& e)
            {
              writeLog_(String("File not found '") + (String)pit->second + "'. Aborting!");
              return INPUT_FILE_NOT_FOUND;
            }
            map_vector.push_back(&(cons_maps[i]));
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
          cons_file.store(out,alignment);
        }
        else
        {
          writeLog_(String("Unknown map type '") + map_type + "' (valid map types are 'feature_map', 'peak_map' and 'consensus_map'. Aborting!");
          return ILLEGAL_PARAMETERS;
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
