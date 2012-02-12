// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page UTILS_IDSplitter IDSplitter
	
	@brief Splits protein/peptide identifications off of annotated data files.

	This performs the reverse operation as IDMapper.
  	
	<B>The command line parameters of this tool are:</B>
	@verbinclude UTILS_IDSplitter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDSplitter
	: public TOPPBase
{
public:

	TOPPIDSplitter()
		: TOPPBase("IDSplitter", "Splits protein/peptide identifications off of annotated data files", false)
		{
		}


protected:

	void removeDuplicates_(vector<PeptideIdentification>& peptides)
		{
			// there is no "PeptideIdentification::operator<", so we can't use a set
			// or sort + unique to filter out duplicates...
			// just use the naive O(n²) algorithm
			vector<PeptideIdentification> unique;
			for (vector<PeptideIdentification>::iterator in_it = peptides.begin();
					 in_it != peptides.end(); ++in_it)
			{
				bool duplicate = false;
				for (vector<PeptideIdentification>::iterator out_it = unique.begin();
						 out_it != unique.end(); ++out_it)
				{
					if (*in_it == *out_it)
					{
						duplicate = true;
						break;
					}
				}
				if (!duplicate) unique.push_back(*in_it);
			}
			peptides.swap(unique);
		}


	void registerOptionsAndFlags_()
		{
			registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
			setValidFormats_("in", StringList::create("mzML,featureXML,consensusXML"));
			registerOutputFile_("out", "<file>", "", "Output file (data without identifications)", false);
			setValidFormats_("out", StringList::create("mzML,featureXML,consensusXML"));
			registerOutputFile_("id_out", "<file>", "", "Output file (identifications)", false);
			setValidFormats_("id_out", StringList::create("idXML"));
			addEmptyLine_();
			addText_("Either 'out' or 'id_out' are required. They can be used together.");
		}


	ExitCodes main_(int , const char**)
		{
			String in = getStringOption_("in"), out = getStringOption_("out"),
				id_out = getStringOption_("id_out");

			if (out.empty() && id_out.empty())
			{
				throw Exception::RequiredParameterNotGiven(__FILE__, __LINE__,
																									 __PRETTY_FUNCTION__, 
																									 "out/id_out");
			}

			vector<ProteinIdentification> proteins;
			vector<PeptideIdentification> peptides;

			FileTypes::Type in_type = FileHandler::getType(in);

			if (in_type == FileTypes::MZML)
			{
				MSExperiment<> experiment;
				MzMLFile().load(in, experiment);
				// what about unassigned peptide IDs?
				for (MSExperiment<>::Iterator exp_it = experiment.begin();
						 exp_it != experiment.end(); ++exp_it)
				{
					peptides.insert(peptides.end(),
													exp_it->getPeptideIdentifications().begin(),
													exp_it->getPeptideIdentifications().end());
					exp_it->getPeptideIdentifications().clear();
				}
				experiment.getProteinIdentifications().swap(proteins);
				if (!out.empty()) 
				{
					addDataProcessing_(experiment, 
														 getProcessingInfo_(DataProcessing::FILTERING));
					MzMLFile().store(out, experiment);
				}
			}

			else if (in_type == FileTypes::FEATUREXML)
			{
				FeatureMap<> features;
				FeatureXMLFile().load(in, features);
				features.getUnassignedPeptideIdentifications().swap(peptides);
				for (FeatureMap<>::Iterator feat_it = features.begin();
						 feat_it != features.end(); ++feat_it)
				{
					peptides.insert(peptides.end(),
													feat_it->getPeptideIdentifications().begin(),
													feat_it->getPeptideIdentifications().end());
					feat_it->getPeptideIdentifications().clear();
				}
				features.getProteinIdentifications().swap(proteins);
				if (!out.empty())
				{
					addDataProcessing_(features, 
														 getProcessingInfo_(DataProcessing::FILTERING));
					FeatureXMLFile().store(out, features);
				}
			}

			else // consensusXML
			{
				ConsensusMap consensus;
				ConsensusXMLFile().load(in, consensus);
				consensus.getUnassignedPeptideIdentifications().swap(peptides);
				for (ConsensusMap::Iterator cons_it = consensus.begin();
						 cons_it != consensus.end(); ++cons_it)
				{
					peptides.insert(peptides.end(),
													cons_it->getPeptideIdentifications().begin(),
													cons_it->getPeptideIdentifications().end());
					cons_it->getPeptideIdentifications().clear();
				}
				consensus.getProteinIdentifications().swap(proteins);
				if (!out.empty()) 
				{
					addDataProcessing_(consensus, 
														 getProcessingInfo_(DataProcessing::FILTERING));
					ConsensusXMLFile().store(out, consensus);
				}
			}

			if (!id_out.empty())
			{
				// IDMapper can match a peptide ID to several overlapping features,
				// resulting in duplicates; this shouldn't be the case for peak data
				if (in_type != FileTypes::MZML) removeDuplicates_(peptides);
				IdXMLFile().store(id_out, proteins, peptides);
			}

			return EXECUTION_OK;
		}
};


int main( int argc, const char** argv )
{
	TOPPIDSplitter tool;
	return tool.main(argc, argv);
}
  
/// @endcond
