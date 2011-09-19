// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Marc Sturm, Clemens Groepl, Steffen Sass $
// --------------------------------------------------------------------------
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithm.h>


#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_FeatureLinkerBase FeatureLinkerBase

	@brief Base class for different FeatureLinker tools.

*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinkerBase
  : public TOPPBase
{

public:
	TOPPFeatureLinkerBase(String name, String description)
		: TOPPBase(name, description)
	{
	}

protected:
	void registerOptionsAndFlags_() // only for "unlabeled" algorithms!
	{
		registerInputFileList_("in", "<files>", StringList(), "input files separated by blanks", true);
		setValidFormats_("in", StringList::create("featureXML,consensusXML"));
		registerOutputFile_("out", "<file>", "", "Output file", true);
		setValidFormats_("out", StringList::create("consensusXML"));
		addEmptyLine_();
		// addText_("Additional parameters for consensusXML input:");
		registerFlag_("keep_subelements", "For consensusXML input only: If set, the sub-features of the inputs are transferred to the output.");
	}


	ExitCodes common_main_(FeatureGroupingAlgorithm* algorithm,
												 bool labeled=false)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		StringList ins;
		if (labeled) ins << getStringOption_("in");
		else ins = getStringList_("in");
		String out = getStringOption_("out");

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		// check if all input files have the correct type
		FileTypes::Type file_type = FileHandler::getType(ins[0]);
		for (Size i = 0; i < ins.size(); ++i)
		{
			if (FileHandler::getType(ins[i]) != file_type)
			{
				writeLog_("Error: All input files must be of the same type!");
				return ILLEGAL_PARAMETERS;
			}
		}

    //-------------------------------------------------------------
    // set up algorithm
    //-------------------------------------------------------------
		Param algorithm_param = getParam_().copy("algorithm:", true);
		writeDebug_("Used algorithm parameters", algorithm_param, 3);
		algorithm->setParameters(algorithm_param);

    //-------------------------------------------------------------
    // perform grouping
    //-------------------------------------------------------------
		// load input
		ConsensusMap out_map;
		if (file_type == FileTypes::FEATUREXML)
		{
			std::vector< FeatureMap<> > maps(ins.size());
			FeatureXMLFile f;
			for (Size i = 0; i < ins.size(); ++i)
			{
				f.load(ins[i], maps[i]);
			}
			for (Size i = 0; i < ins.size(); ++i)
			{
				out_map.getFileDescriptions()[i].filename = ins[i];
				out_map.getFileDescriptions()[i].size = maps[i].size();
				out_map.getFileDescriptions()[i].unique_id = maps[i].getUniqueId();
			}
			// exception for "labeled" algorithms: copy file descriptions
			if (labeled)
			{
				out_map.getFileDescriptions()[1] = out_map.getFileDescriptions()[0];
				out_map.getFileDescriptions()[0].label = "light";
				out_map.getFileDescriptions()[1].label = "heavy";
			}
			// group
			algorithm->group(maps, out_map);
		}
		else
		{
			std::vector<ConsensusMap> maps(ins.size());
			ConsensusXMLFile f;
			for (Size i = 0; i < ins.size(); ++i)
			{
				f.load(ins[i], maps[i]);
			}
			// group
			algorithm->group(maps, out_map);

			// set file descriptions:
			bool keep_subelements = getFlag_("keep_subelements");
			if (!keep_subelements)
			{
				for (Size i = 0; i < ins.size(); ++i)
				{
					out_map.getFileDescriptions()[i].filename = ins[i];
					out_map.getFileDescriptions()[i].size = maps[i].size();
					out_map.getFileDescriptions()[i].unique_id = maps[i].getUniqueId();
				}
			}
			else
			{
				// components of the output map are not the input maps themselves, but
				// the components of the input maps:
				algorithm->transferSubelements(maps, out_map);
			}
		}

		// assign unique ids
		out_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

		// annotate output with data processing info
		addDataProcessing_(out_map, getProcessingInfo_(DataProcessing::FEATURE_GROUPING));

		// write output
		ConsensusXMLFile().store(out, out_map);

		// some statistics
		map<Size, UInt> num_consfeat_of_size;
		for (ConsensusMap::const_iterator cmit = out_map.begin(); cmit != out_map.end(); ++cmit)
		{
			++num_consfeat_of_size[cmit->size()];
		}

		LOG_INFO << "Number of consensus features:" << endl;
		for (map<Size, UInt>::reverse_iterator i = num_consfeat_of_size.rbegin(); i != num_consfeat_of_size.rend(); ++i)
		{
			LOG_INFO << "  of size " << setw(2) << i->first << ": " << setw(6) << i->second << endl;
		}
		LOG_INFO << "  total:      " << setw(6) << out_map.size() << endl;

		return EXECUTION_OK;
	}
};

/// @endcond
