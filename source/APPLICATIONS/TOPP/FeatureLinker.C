// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
	@page TOPP_FeatureLinker FeatureLinker

	@brief Groups corresponding features in one map or across maps.

<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=4> \f$ \longrightarrow \f$ FeatureLinker \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_FeatureFinder </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=2> @ref TOPP_MapAligner </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_SeedListGenerator </td>
		</tr>
	</table>
</CENTER>


	This tool provides algorithms for grouping corresponding features in @ref OpenMS::FeatureGroupingAlgorithmLabeled "isotope-labeled" and @ref OpenMS::FeatureGroupingAlgorithmUnlabeled "label-free" experiments. (Click on the links for detailed information including algorithm-specific parameters.)

	FeatureLinker takes one or several feature maps (featureXML files) and stores the corresponding features in a consensus map (consensusXML file). Feature maps can be created from MS experiments (peak data) using the @ref TOPP_FeatureFinder.

	It is assumed that major retention time distortions are corrected before applying this tool. Use @ref TOPP_MapAligner to do that on the peak or feature level.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_FeatureLinker.cli

	For the parameters of the algorithm section see the algorithms documentation: @n
		@ref OpenMS::FeatureGroupingAlgorithmUnlabeled "algorithm unlabeled" @n
		@ref OpenMS::FeatureGroupingAlgorithmLabeled "algorithm labeled" @n
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPFeatureLinker
  : public TOPPBase
{

public:
	TOPPFeatureLinker()
		: TOPPBase("FeatureLinker","Groups corresponding features in one map or across maps.")
	{
	}

protected:
	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>",StringList(),"input files separated by blanks",true);
		setValidFormats_("in",StringList::create("featureXML,consensusXML"));
		registerOutputFile_("out","<file>","","Output file",true);
		setValidFormats_("out",StringList::create("consensusXML"));
		registerStringOption_("type","<name>","","Feature grouping algorithm type",true);
		setValidStrings_("type", getToolList()[toolName_()]);

		registerSubsection_("algorithm","Algorithm parameters section");
	}

	Param getSubsectionDefaults_(const String& /*section*/) const
	{
		String type = getStringOption_("type");
		FeatureGroupingAlgorithm* algo = Factory<FeatureGroupingAlgorithm>::create(type);
		Param p = algo->getParameters();
		delete algo;
		return p;
	}

	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
		StringList ins = getStringList_("in");

		String out = getStringOption_("out");

		String type = getStringOption_("type");

		//-------------------------------------------------------------
		// check for valid input
		//-------------------------------------------------------------
		//check if all input files have the correct type
		FileTypes::Type file_type= FileHandler::getType(ins[0]);
		if (type=="unlabeled_qt")
		{
			for (Size i=0;i<ins.size();++i)
			{
				if (FileHandler::getType(ins[i])!=file_type)
				{
					writeLog_("Error: All input files must be of same type!");
					return ILLEGAL_PARAMETERS;
				}
			}
		}
		else
		{
			for (Size i=0;i<ins.size();++i)
			{
				if (FileHandler::getType(ins[i])!=FileTypes::FEATUREXML)
				{
					writeLog_("Error: All input files must be of type FeatureXML!");
					return ILLEGAL_PARAMETERS;
				}
			}
		}

    //-------------------------------------------------------------
    // set up algorithm
    //-------------------------------------------------------------
    FeatureGroupingAlgorithm* algorithm = Factory<FeatureGroupingAlgorithm>::create(type);
		Param algorithm_param = getParam_().copy("algorithm:", true);
		writeDebug_("Used algorithm parameters",algorithm_param, 3);
		algorithm->setParameters(algorithm_param);

    //-------------------------------------------------------------
    // perform grouping
    //-------------------------------------------------------------
		//load input

		ConsensusMap out_map;
		if (file_type==FileTypes::FEATUREXML)
		{
			std::vector< FeatureMap<> > maps(ins.size());
			FeatureXMLFile f;
			for (Size i=0; i<ins.size(); ++i)
			{
				f.load(ins[i], maps[i]);
			}
			for (Size i=0; i<ins.size(); ++i)
			{
				out_map.getFileDescriptions()[i].filename = ins[i];
				out_map.getFileDescriptions()[i].size = maps[i].size();
				out_map.getFileDescriptions()[i].unique_id = maps[i].getUniqueId();
			}
			//Exception for 'labeled' algorithms: copy file descriptions
			if (type=="labeled")
			{
				out_map.getFileDescriptions()[1] = out_map.getFileDescriptions()[0];
				out_map.getFileDescriptions()[0].label = "light";
				out_map.getFileDescriptions()[1].label = "heavy";
			}
			// group
			algorithm->group(maps,out_map);
		}
		else
		{
			std::vector<ConsensusMap> maps(ins.size());
			ConsensusXMLFile f;
			for (Size i=0; i<ins.size(); ++i)
			{
				f.load(ins[i], maps[i]);
			}
			if (out_map.getFileDescriptions().empty())
			{
				for (Size i=0; i<ins.size(); ++i)
				{
					out_map.getFileDescriptions()[i].filename = ins[i];
					out_map.getFileDescriptions()[i].size = maps[i].size();
					out_map.getFileDescriptions()[i].unique_id = maps[i].getUniqueId();
				}
			}
			// group
			algorithm->group(maps,out_map);
		}

		//set file names

		// assign unique ids
		out_map.applyMemberFunction(&UniqueIdInterface::setUniqueId);

		// annotate output with data processing info
		addDataProcessing_(out_map, getProcessingInfo_(DataProcessing::FEATURE_GROUPING));

		// write output
		ConsensusXMLFile().store(out,out_map);

		// some statistics
		map<Size,UInt> num_consfeat_of_size;
		for ( ConsensusMap::const_iterator cmit = out_map.begin(); cmit != out_map.end(); ++cmit )
		{
			++num_consfeat_of_size[cmit->size()];
		}

		LOG_INFO << "Number of consensus features:" << endl;
		for ( map<Size,UInt>::reverse_iterator i = num_consfeat_of_size.rbegin(); i != num_consfeat_of_size.rend(); ++i )
		{
			LOG_INFO << "  of size " << setw(2) << i->first << ": " << setw(6) << i->second << endl;
		}
		LOG_INFO << "  total:      " << setw(6) << out_map.size() << endl;


		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
  TOPPFeatureLinker tool;
  return tool.main(argc,argv);
}

/// @endcond
