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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona $
// --------------------------------------------------------------------------

#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <algorithm>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDConflictResolver IDConflictResolver
	
	@brief Resolves ambiguous annotations of features with peptide identifications.
 
	<CENTER>
	<table>
		<tr>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential predecessor tools </td>
			<td VALIGN="middle" ROWSPAN=3> \f$ \longrightarrow \f$ IDConflictResolver \f$ \longrightarrow \f$</td>
			<td ALIGN = "center" BGCOLOR="#EBEBEB"> potential successor tools </td>
		</tr>
		<tr>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_IDMapper </td>
			<td VALIGN="middle" ALIGN = "center" ROWSPAN=1> @ref TOPP_TextExporter </td>
		</tr>
		<tr>
		  <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_FeatureLinkerUnlabeled @n (or another feature grouping algorithm) </td>
		  <td VALIGN="middle" ALIGN="center" ROWSPAN=1> @ref TOPP_ProteinQuantifier </td>
		</tr>
	</table>
	</CENTER>
 	
	The peptide identifications are filtered so that only one identification with a single hit (with the best score) is associated to each feature. (If two IDs have the same best score, either one of them may be selected.)

	This step may be useful before applying @ref TOPP_ProteinQuantifier "ProteinQuantifier", because features with ambiguous annotation are not considered for the quantification.
  	
	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDConflictResolver.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDConflictResolver
	: public TOPPBase
{
public:

	TOPPIDConflictResolver()
		: TOPPBase("IDConflictResolver", "Resolves ambiguous annotations of features with peptide identifications")
	{
	}

protected:

	// compare peptide IDs by score of best hit (hits must be sorted first!)
	// (note to self: the "static" is necessary to avoid cryptic "no matching
	// function" errors from gcc when the comparator is used below)
	static bool compareIDs_(const PeptideIdentification& left, 
													const PeptideIdentification& right)
	{
		if (left.getHits()[0].getScore() < right.getHits()[0].getScore())
		{
			return true;
		}
		return false;
	}


	void resolveConflict_(vector<PeptideIdentification>& peptides)
	{
		if (peptides.empty()) return;
		for (vector<PeptideIdentification>::iterator pep_id = peptides.begin();
				 pep_id != peptides.end(); ++pep_id)
		{
			pep_id->sort();
		}
		vector<PeptideIdentification>::iterator pos;
		if (peptides[0].isHigherScoreBetter()) // find highest-scoring ID
		{
			pos = max_element(peptides.begin(), peptides.end(), compareIDs_);
		}
		else // find lowest-scoring ID
		{
			pos = min_element(peptides.begin(), peptides.end(), compareIDs_);
		}
		peptides[0] = *pos;
		peptides.resize(1);
		// remove all but the best hit:
		vector<PeptideHit> best_hit(1, peptides[0].getHits()[0]);
		peptides[0].setHits(best_hit);
	}


	void registerOptionsAndFlags_()
	{
		registerInputFile_("in", "<file>", "", "Input file (data annotated with identifications)");
		setValidFormats_("in", StringList::create("featureXML,consensusXML"));
		registerOutputFile_("out", "<file>", "", "Output file (data with one peptide identification per feature)");
		setValidFormats_("out", StringList::create("featureXML,consensusXML"));
	}


	ExitCodes main_(int , const char**)
	{
		String in = getStringOption_("in"), out = getStringOption_("out");
		FileTypes::Type in_type = FileHandler::getType(in);
		if (in_type == FileTypes::FEATUREXML)
		{
			FeatureMap<> features;
			FeatureXMLFile().load(in, features);
			for (FeatureMap<>::Iterator feat_it = features.begin();
					 feat_it != features.end(); ++feat_it)
			{
				resolveConflict_(feat_it->getPeptideIdentifications());
			}
			addDataProcessing_(features, 
												 getProcessingInfo_(DataProcessing::FILTERING));
			FeatureXMLFile().store(out, features);
		}
		else // consensusXML
		{
			ConsensusMap consensus;
			ConsensusXMLFile().load(in, consensus);
			for (ConsensusMap::Iterator cons_it = consensus.begin();
					 cons_it != consensus.end(); ++cons_it)
			{
				resolveConflict_(cons_it->getPeptideIdentifications());
			}
			addDataProcessing_(consensus, 
												 getProcessingInfo_(DataProcessing::FILTERING));
			ConsensusXMLFile().store(out, consensus);
		}
		return EXECUTION_OK;
	}
};


int main(int argc, const char** argv)
{
	TOPPIDConflictResolver tool;
	return tool.main(argc, argv);
}
  
/// @endcond
