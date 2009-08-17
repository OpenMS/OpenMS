// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/Map.h>
#include <fstream>
#include <sstream>

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
@page UTILS_CaapConvert CaapConvert

	@brief Convert a CAAP ground truth file into consensusXML format.

	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_CaapConvert.cli
*/


// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

int main( int argc, const char** argv )
{
	if ( argc != 3 && argc != 4 )
	{
		std::cout <<
			"Usage:  " << argv[0] << " input output [-v]\n"
			"\n"
			"where:\n"
			"  input    is a ground truth file as described on the CAAP web page\n"
			"  output   is the result in consensusXML format as described in the OpenMS docu.\n"
			"  [supply optional third argument -v for verbose output]\n"
			"\n"
			"See the paper:\n"
			"\"Critical assessment of alignment procedures for LC-MS proteomics and metabolomics measurements\"\n"
			"Eva Lange, Ralf Tautenhahn, Steffen Neumann, Clemens Groepl\n"
			"BMC Bioinformatics 2008, 9:375.\n"
			"http://dx.doi.org/10.1186/1471-2105-9-375\n"
			;
		return 1;
	}

	int verbose = 0;
	if ( argc == 4 )
	{
		if ( std::string(argv[3]) =="-v" )
		{
			verbose = 1;
		}
		else
		{
			verbose = 2;
		}
	}

#define VERBOSEMSG1(bla) if ( verbose >= 1 ) { std::cout << bla; }
#define VERBOSEMSG2(bla) if ( verbose >= 2 ) { std::cout << bla; }

	const char * const argv_input  = argv[1];
	const char * const argv_output = argv[2];

	std::fstream input(argv_input);
	std::string line;
	std::string map_id_str;
	OpenMS::Size map_id_num = 0;
	typedef std::map<std::string,OpenMS::Size> MapType;
	MapType map_filename_to_map_index;
	double score;
	double intensity;
	double retention_time;
	double mass_to_charge;
	OpenMS::ConsensusMap consensus_map;
	OpenMS::ConsensusFeature consensus_feature;
	OpenMS::FeatureHandle feature_handle;
	for (;;)
	{
		line.clear();
		std::getline(input,line);
		if (!input) break;
		VERBOSEMSG2("line: " << line << std::endl);
		std::stringstream linestream(line);
		consensus_feature.clear();
		for (;;)
		{
			linestream >> map_id_str >> score >> intensity >> retention_time >> mass_to_charge;
			if (!linestream) break;
			MapType::const_iterator map_iter = map_filename_to_map_index.find(map_id_str);
			if ( map_iter != map_filename_to_map_index.end() )
			{
				map_id_num = map_iter->second;
			}
			else
			{
				map_filename_to_map_index.insert(std::make_pair(map_id_str,(map_id_num=map_filename_to_map_index.size())));
			}
			VERBOSEMSG2
				(
				 "CE: " <<
				 map_id_str << " " <<
				 map_id_num << " " <<
				 score << " " <<
				 intensity << " " <<
				 retention_time << " " <<
				 mass_to_charge <<
				 std::endl
				);
			feature_handle.setMapIndex(map_id_num);
			// We currently do not trace the element indices back to the original feature maps.
			feature_handle.setElementIndex(0);
			feature_handle.setIntensity(intensity);
			feature_handle.setRT(retention_time);
			feature_handle.setMZ(mass_to_charge);
#if 1
			// We need to cast away the overloaded insert that checks for reusage of
			// element indices because we want to cowardly ignore a few duplicates.
			bool is_no_duplicate = static_cast<OpenMS::ConsensusFeature::HandleSetType&>(consensus_feature).insert(feature_handle).second;
			if (!is_no_duplicate)
			{
				VERBOSEMSG1
					(
					 "\nNote: cowardly ignoring a duplicate feature_handle:\n" << feature_handle <<
					 "---------- /FeatureHandle ----------------\n"
					 "In this line:   " << line << "\n\n"
					);
			}
#else
			// Fussy mode ... ;)  (will throw exceptions)
			consensus_feature.insert(feature_handle);
#endif
		}
		consensus_feature.computeConsensus();
		consensus_feature.setQuality(1);
		consensus_map.push_back(consensus_feature);
	}

	OpenMS::ConsensusMap::FileDescription file_description;
	VERBOSEMSG2("map_id_numbers:");
	for ( MapType::const_iterator map_iter = map_filename_to_map_index.begin();
				map_iter != map_filename_to_map_index.end();
				++map_iter
			)
	{
		VERBOSEMSG2(" " << map_iter->first << ':' << map_iter->second);
		file_description.filename = map_iter->first;
		file_description.label = "";
		file_description.size = 1; // element_index is always 0
		consensus_map.getFileDescriptions()[map_iter->second] = file_description;
	}
	VERBOSEMSG2(std::endl);

	OpenMS::ConsensusXMLFile consensus_xml_file;
	consensus_xml_file.store(argv_output,consensus_map);

}

/// @endcond
