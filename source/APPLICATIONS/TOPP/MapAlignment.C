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
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/GridFile.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

typedef FeatureMap< Feature > FeatureMapType;
typedef DPeakArray<Peak2D > PeakArrayType;
typedef ConsensusFeature< FeatureMapType > ConsensusFeatureType;
typedef ConsensusPeak< PeakArrayType > ConsensusPeakType;
typedef ConsensusMap< ConsensusFeatureType > ConsensusMapType;

#define DEBUG_CONSENSUS
#undef DEBUG_CONSENSUS


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
   
   @todo Talk through the concept of MapDewarper, MapMatcher, MapAlignment (Marc, Eva)
	 @todo document parameters! (Eva)   
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
		registerStringOption_("out","<file>","Consensus.xml","output consensusXML file name",false);
      
		addEmptyLine_();
		addText_("This application implements an algorithm for the alignment of multiple maps.\n"
						 "It accepts feature maps (in featureXML), peak maps (in mzData) or consensus maps (in ConsensusXML)\n"
						 "The output of the MapAlignment tool depends on the type of the input maps. \n"
             "In case of peak maps it returns the warping functions that map each input map onto the reference map along with the dewarped maps itself.\n"
             "The alignment of feature or consensus maps result in a consensus map, which contains all grouped elements.\n"  
             "Note: If you use consensus maps , the consensus elements are used as normal elements and you will\n"
						 "loose the former consensus information.");

		addEmptyLine_();
		addText_("This tool requires an INI file with at least the names of the input files and the map_type.\n"
						 "Parameters for the alignment algorithm can be given only in the 'algorithm' seciton  of the INI file.\n");			
		registerSubsection_("algorithm","Algorithm parameters section");
		registerSubsection_("file_names","Input file name section");
	}
	
	 Param getSubsectionDefaults_(const String& section) const
    {
      Param tmp;
		
			if (section == "algorithm")
				{
					tmp.setValue("map_type","feature_map");
					tmp.setValue("number_buckets:RT",1);
      		tmp.setValue("number_buckets:MZ",1);
					tmp.setValue("matching_algorithm:type","poseclustering_pairwise");
					tmp.setValue("matching_algorithm:superimposer:type","poseclustering_affine");
					tmp.insert("matching_algorithm:superimposer",Factory<BaseSuperimposer<> >::create("poseclustering_affine")->getDefaults());
					tmp.setValue("matching_algorithm:pairfinder:type","DelaunayPairFinder");
					tmp.insert("matching_algorithm:pairfinder",Factory<BasePairFinder<> >::create("DelaunayPairFinder")->getDefaults());
					tmp.insert("consensus_algorithm",Factory<BasePairFinder<> >::create("DelaunayPairFinder")->getDefaults());
				}
				if (section == "file_names")
				{
					tmp.setValue("1","feature_map_1.xml");
					tmp.setValue("2","feature_map_2.xml");
					tmp.setValue("3","feature_map_3.xml");
					tmp.setValue("4","feature_map_4.xml");
					tmp.setValue("5","feature_map_5.xml");
				}
			return tmp;
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
		Param::ParamIterator pit = files_param.begin();

		String map_type = getParam_().getValue("algorithm:map_type");
		//-------------------------------------------------------------
		// loading input and initialize the alignment object
		//-------------------------------------------------------------
		if (map_type == "feature_map")
      {
				if (out == "") 
					{
						writeLog_("No name for the output consensus map is given! Please specify the \"out\" option. Aborting!");
						return MISSING_PARAMETERS;
					}
        StarAlignment< ConsensusFeatureType > alignment;
        alignment.setParameters(mapali_param);
        FeatureXMLFile feature_file;
        std::vector< String > file_names;
        // Vector for the feature maps
        std::vector< FeatureMapType > feature_maps(files_param.size());

        // Reference to the map vector of the alignment object
        std::vector< FeatureMapType* >& map_vector = alignment.getElementMapVector();
        unsigned int i=0;
        while (pit != files_param.end())
					{
						file_names.push_back(pit->value);
						// load the feature file into a feature_map
						try
							{
								feature_file.load(pit->value, feature_maps[i]);
							}
						catch(Exception::FileNotFound& e)
							{
								writeLog_(String("File not found '") + String(pit->value) + "'. Aborting!");
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
        
#ifdef DEBUG_CONSENSUS
        std::ofstream out_pairs("MapAlignment_pairs.dat",std::ios::out);
        const ConsensusMap<ConsensusFeatureType>& final_consensus_map_(alignment.getFinalConsensusMap());
        for (UInt i = 0; i < final_consensus_map_.size(); ++i)
        {
          bool ref = false;
          std::vector<const ConsensusFeatureType::ElementType*> features(1);

          const ConsensusFeatureType* c = &(final_consensus_map_[i]);
          for (ConsensusFeatureType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
          {
            if (it->getMapIndex() == alignment.getReferenceMapIndex())
            {
              ref = true;
              features[0] = &(it->getElement());
            }
            else
            {
              features.push_back(&(it->getElement()));
            }
          }
          if (ref)
          {
            out_pairs
                << features[0]->getIntensity() << ' '
                << features[0]->getRT() << ' '
                << features[0]->getMZ() << ' ';

            UInt j=1;
            for (; j < features.size(); ++j)
            {
              out_pairs
                  << features[j]->getIntensity() << ' '
                  << features[j]->getRT() << ' '
                  << features[j]->getMZ() << ' ';
            }
            for (;j < alignment.getElementMapVector().size(); ++j)
            {
              out_pairs
                  << 0 << ' '
                  << 0 << ' '
                  << 0 << ' ';
            }
            out_pairs << std::endl;
          }
        }
#endif       
        
        alignment.merge();

#ifdef DEBUG_CONSENSUS       
        std::ofstream out_pairs_2("MapAlignment_pairs_merged.dat",std::ios::out);
        const ConsensusMap<ConsensusFeatureType>& final_consensus_map_2(alignment.getFinalConsensusMap());
        for (UInt i = 0; i < final_consensus_map_2.size(); ++i)
        {
          bool ref = false;
          std::vector<const ConsensusFeatureType::ElementType*> features(1);

          const ConsensusFeatureType* c = &(final_consensus_map_2[i]);
          for (ConsensusFeatureType::Group::const_iterator it = c->begin(); it != c->end(); ++it)
          {
            if (it->getMapIndex() == alignment.getReferenceMapIndex())
            {
              ref = true;
              features[0] = &(it->getElement());
            }
            else
            {
              features.push_back(&(it->getElement()));
            }
          }
          if (ref)
          {
            out_pairs_2
                << features[0]->getIntensity() << ' '
                << features[0]->getRT() << ' '
                << features[0]->getMZ() << ' ';

            UInt j=1;
            for (; j < features.size(); ++j)
            {
              out_pairs_2
                  << features[j]->getIntensity() << ' '
                  << features[j]->getRT() << ' '
                  << features[j]->getMZ() << ' ';
            }
            for (;j < alignment.getElementMapVector().size(); ++j)
            {
              out_pairs_2
                  << 0 << ' '
                  << 0 << ' '
                  << 0 << ' ';
            }
            out_pairs_2 << std::endl;
          }
        }
#endif

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
					alignment.setParameters(mapali_param);
					MzDataFile mzdata_file;
					mzdata_file.setLogType(log_type_);
					std::vector< String > file_names;
          // Vector for the feature maps
					std::vector< PeakArrayType > peak_maps(files_param.size());

          // Reference to the map vector of the alignment object
					std::vector< PeakArrayType* >& map_vector = alignment.getElementMapVector();
					unsigned int i=0;
					while (pit != files_param.end())
						{
							file_names.push_back(pit->value);
							// load the feature file into a feature_maps
							PeakMap ms_exp;

							try
								{
									mzdata_file.load(pit->value, ms_exp);
								}
							catch(Exception::FileNotFound& e)
								{
									writeLog_(String("File not found '") + (String)pit->value + "'. Aborting!");
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

					UInt ref_index = alignment.getReferenceMapIndex();
					writeLog_("File " + String(file_names[ref_index]) + " is the reference map of the starwise alignment.");
					GridFile grid_file;
					for (UInt m = 0; m < i; ++m)
						{
							PeakArrayType& dewarped_map = peak_maps[m];
							if (m != ref_index)
								{
									// store the transformation
									String file_name(file_names[m]);
									file_name.trim();
									std::vector< String > substrings;
									file_name.split('.',substrings);
									file_name.implode(substrings.begin(),substrings.end()-1);
									String file_name_grid(file_name + ".grid");
									String file_name_dewarped(file_name + "_dewarped.mzData");
									writeLog_("Store the transformation, which maps " + file_name_dewarped + " onto the reference map in " + file_name_grid + '.');
									grid_file.store(file_name_grid,alignment.getTransformationVector()[m]);

									// iterate over all Elements...
									UInt n = map_vector[m]->size();
									for (UInt j = 0; j < n; ++j)
										{
											// Test in which cell this element is included
											// and apply the corresponding transformation
											Grid::const_iterator grid_it = (alignment.getTransformationVector()[m]).begin();
											while (grid_it != (alignment.getTransformationVector()[m]).end() )
												{
													if (grid_it->encloses(dewarped_map[j].getPosition()) )
														{
															LinearMapping* mapping_rt = dynamic_cast<LinearMapping* >(grid_it->getMappings()[RawDataPoint2D::RT]);
															LinearMapping* mapping_mz = dynamic_cast<LinearMapping* >(grid_it->getMappings()[RawDataPoint2D::MZ]);

															DPosition<2> pos = dewarped_map[j].getPosition();

															mapping_rt->apply(pos[RawDataPoint2D::RT]);
															mapping_mz->apply(pos[RawDataPoint2D::MZ]);

															dewarped_map[j].setPosition(pos);
														}
													grid_it++;
												} // end while (grid)
										} // end for
            
									writeLog_("Write dewarped map to " +  file_name_dewarped + '.');
									PeakMap ms_exp;
									ms_exp.set2DData(dewarped_map);
									mzdata_file.store(file_name_dewarped,ms_exp);
								}
						}
				}
      else if (map_type == "consensus_map")
				{
					if (out == "") 
						{
							writeLog_("No name for the output consensus map is given! Please specify the \"out\" option. Aborting!");
							return MISSING_PARAMETERS;
						}
					StarAlignment< ConsensusFeature< ConsensusMapType > > alignment;
					alignment.setParameters(mapali_param);
					
					ConsensusXMLFile cons_file;
					std::vector< String > file_names;
          // Vector for the feature maps
					std::vector< ConsensusMapType > cons_maps(files_param.size());

          // Reference to the map vector of the alignment object
					std::vector< ConsensusMapType* >& map_vector = alignment.getElementMapVector();
					unsigned int i=0;
					while (pit != files_param.end())
						{
							file_names.push_back(pit->value);
							// load the feature file into a feature_map
							try
								{
									cons_file.load(pit->value, cons_maps[i], false);
								}
							catch(Exception::FileNotFound& e)
								{
									writeLog_(String("File not found '") + (String)pit->value + "'. Aborting!");
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
