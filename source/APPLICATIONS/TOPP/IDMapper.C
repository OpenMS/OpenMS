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
// $Maintainer: Andreas Bertsch $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/config.h>

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/ANALYSIS/ID/IDMapper.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDMapper IDMapper

	Assigns protein/peptide identifications to feature or consensus features.

	This tool is typically used before @ref TOPP_ConsensusID.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDMapper.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMapper
	: public TOPPBase
{

	public:

		TOPPIDMapper()
			: TOPPBase("IDMapper", "Assigns protein/peptide identifications to feature or consensus features.")
		{
		}

	protected:

		void registerOptionsAndFlags_()
		{
			registerInputFile_("id", "<file>", "", "Protein/peptide identifications file");
			setValidFormats_("id",StringList::create("idXML"));
			registerInputFile_("in", "<file>", "", "Protein/peptide identifications file");
			setValidFormats_("in",StringList::create("featureXML,consensusXML"));
			registerOutputFile_("out", "<file>", "", "Output file (the format depends on the input file format).");
			setValidFormats_("out",StringList::create("featureXML,consensusXML"));

			addEmptyLine_();
			IDMapper mapper;
			Param p = mapper.getParameters();
			registerDoubleOption_("rt_delta","<value>",p.getValue("rt_delta"), "Maximum allowed RT delta (seconds) between identification and peak/feature.", false);
			setMinFloat_("rt_delta",0.0);
			registerDoubleOption_("mz_delta","<value>",p.getValue("mz_delta"), "Maximum allowed m/z delta (ppm or Da) between identification and peak/feature.", false);
			setMinFloat_("mz_delta",0.0);
			registerStringOption_("mz_measure","<String>",p.getEntry("mz_measure").valid_strings[0],"Unit of mz_delta", false);
			setValidStrings_("mz_measure", p.getEntry("mz_measure").valid_strings);
			registerStringOption_("mz_reference","<String>",p.getEntry("mz_reference").valid_strings[1],"Method to determine m/z of identification", false);
			setValidStrings_("mz_reference", p.getEntry("mz_reference").valid_strings);
			registerFlag_("use_centroids","[FeatureMap only] use RT&MZ coordinate instead of convex hull");
			registerFlag_("use_subelements","[ConsensusMap only] use RT&MZ coordinate of sub-features instead of consensus RT&MZ");
		}

		ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in");
			FileTypes::Type in_type = FileHandler::getType(in);
			String out = getStringOption_("out");

			//----------------------------------------------------------------
			// load idXML
			//----------------------------------------------------------------
			vector<ProteinIdentification> protein_ids;
			vector<PeptideIdentification> peptide_ids;
			String document_id;
			IdXMLFile().load(getStringOption_("id"),protein_ids,peptide_ids, document_id);

			//----------------------------------------------------------------
			//create mapper
			//----------------------------------------------------------------
			IDMapper mapper;
			Param p = mapper.getParameters();
			p.setValue("rt_delta", getDoubleOption_("rt_delta"));
			p.setValue("mz_delta", getDoubleOption_("mz_delta"));
			p.setValue("mz_measure",getStringOption_("mz_measure"));
			p.setValue("mz_reference",getStringOption_("mz_reference"));
			mapper.setParameters(p);

			//----------------------------------------------------------------
			// consensusXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::CONSENSUSXML)
			{
				ConsensusXMLFile file;
				ConsensusMap map;
				file.load(in,map);
				
				bool measure_from_subelements=getFlag_("use_subelements");
				
				mapper.annotate(map,peptide_ids,protein_ids,measure_from_subelements);
				
				//annotate output with data processing info
				addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));
				
				file.store(out,map);
			}

			//----------------------------------------------------------------
			// featureXML
			//----------------------------------------------------------------
			if (in_type == FileTypes::FEATUREXML)
			{
				FeatureMap<> map;
				FeatureXMLFile file;
				file.load(in,map);

				bool measure_from_centroids=getFlag_("use_centroids");
				
				mapper.annotate(map,peptide_ids,protein_ids,measure_from_centroids);
				
				//annotate output with data processing info
				addDataProcessing_(map, getProcessingInfo_(DataProcessing::IDENTIFICATION_MAPPING));
				
				file.store(out,map);
			}

			return EXECUTION_OK;
		}

};


int main( int argc, const char** argv )
{
	TOPPIDMapper tool;
	return tool.main(argc, argv);
}

/// @endcond
