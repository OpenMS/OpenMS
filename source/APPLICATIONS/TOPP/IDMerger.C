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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/APPLICATIONS/TOPPBase.h>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDMerger IDMerger
	
	@brief Merges several idXML files into one idXML file.
	
	You can merge an unlimited number of files into one idXML file.

	With the @p -pepxml_protxml option, results from corresponding PeptideProphet and ProteinProphet runs can be combined. In this case, exactly two idXML files are expected as input: one containing data from a pepXML file, and the other containing data from a protXML file that was created based on the pepXML (meaningful results can only be obtained for matching files!). The @ref TOPP_IDFileConverter can be used to convert pepXML or protXML to idXML.
	
	This tool is typically applied before @ref TOPP_ConsensusID or @ref TOPP_IDMapper.

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDMerger.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPIDMerger
	: public TOPPBase
{
 public:
	TOPPIDMerger()
		: TOPPBase("IDMerger","Merges several protein/peptide identification files into one file.")
	{
			
	}
	
 protected:
	void mergePepXMLProtXML_(StringList filenames, vector<ProteinIdentification>&
													 proteins, vector<PeptideIdentification>& peptides)
		{
			IdXMLFile idxml;
			idxml.load(filenames[0], proteins, peptides);
			vector<ProteinIdentification> pepxml_proteins, protxml_proteins;
			vector<PeptideIdentification> pepxml_peptides, protxml_peptides;

			if (proteins[0].getProteinGroups().empty())
			{	// first idXML contains data from the pepXML
				proteins.swap(pepxml_proteins);
				peptides.swap(pepxml_peptides);
				idxml.load(filenames[1], protxml_proteins, protxml_peptides);
				if (protxml_proteins[0].getProteinGroups().empty())
				{
					throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "None of the input files seems to be derived from a protXML file (information about protein groups is missing).");
				}
			}
			else
			{ // first idXML contains data from the protXML
				proteins.swap(protxml_proteins);
				peptides.swap(protxml_peptides);
				idxml.load(filenames[1], pepxml_proteins, pepxml_peptides);
			}

			if ((protxml_peptides.size() > 1) || (protxml_proteins.size() > 1))
			{
				throw Exception::InvalidParameter(__FILE__, __LINE__, __PRETTY_FUNCTION__, "The idXML derived from a protXML file should contain only one 'ProteinIdentification' and one 'PeptideIdentification' instance.");
			}

			// peptide information comes from the pepXML (additional information in
			// the protXML - adapted peptide hit score, "is_unique", "is_contributing"
			// - is not transferred):
			peptides.swap(pepxml_peptides);

			// prepare scores and coverage values of protein hits from the protXML:
			map<String, pair<DoubleReal, DoubleReal> > hit_values;
			ProteinIdentification& protein = protxml_proteins[0];
			for (vector<ProteinHit>::iterator hit_it = protein.getHits().begin();
					 hit_it != protein.getHits().end(); ++hit_it)
			{
				hit_values[hit_it->getAccession()] = make_pair(hit_it->getScore(),
																											 hit_it->getCoverage());
			}
			
			// merge protein information:
			proteins.swap(pepxml_proteins);
			for (vector<ProteinIdentification>::iterator prot_it = proteins.begin();
					 prot_it != proteins.end(); ++prot_it)
			{
				prot_it->getProteinGroups() = protein.getProteinGroups();
				prot_it->getIndistinguishableProteins() = 
					protein.getIndistinguishableProteins();
				// TODO: since a protXML file can integrate data from several protein
				// identification runs, the protein groups/indistinguishable proteins
				// that we write to one identification run could contain references to
				// proteins that are not observed in this run, but in others; also, some
				// protein hits without enough evidence may not occur in the protXML
				// (thus also not in the protein groups) - clean this up?

				prot_it->setScoreType(protein.getScoreType());
				prot_it->setHigherScoreBetter(protein.isHigherScoreBetter());
				prot_it->setSignificanceThreshold(protein.getSignificanceThreshold());

				for (vector<ProteinHit>::iterator hit_it = prot_it->getHits().begin();
						 hit_it != prot_it->getHits().end(); ++hit_it)
				{
					map<String, pair<DoubleReal, DoubleReal> >::const_iterator pos = 
						hit_values.find(hit_it->getAccession());
					if (pos == hit_values.end())
					{
						hit_it->setScore(-1);
					}
					else
					{
						hit_it->setScore(pos->second.first);
						hit_it->setCoverage(pos->second.second);
					}
				}
			}
		}

	void registerOptionsAndFlags_()
	{
		registerInputFileList_("in","<files>",StringList(),"two or more input files separated by blanks");
		setValidFormats_("in",StringList::create("idXML"));
		registerOutputFile_("out","<file>","","output file ");
		setValidFormats_("out",StringList::create("idXML"));
		registerFlag_("pepxml_protxml", "Merge idXML files derived from a pepXML and corresponding protXML file.\nExactly two input files are expected in this case.");
	}
	
	ExitCodes main_(int , const char**)
	{
		//-------------------------------------------------------------
		// parameter handling
		//-------------------------------------------------------------
	
		StringList file_names = getStringList_("in");
		String out = getStringOption_("out");
		
		if (file_names.size() < 2)
		{
			writeLog_("Less than two filenames given. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}

		bool pepxml_protxml = getFlag_("pepxml_protxml");
		if (pepxml_protxml && (file_names.size() != 2))
		{
			writeLog_("Exactly two filenames expected for option 'pepxml_protxml'. Aborting!");
			printUsage_();
			return ILLEGAL_PARAMETERS;
		}
				
		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------

		vector<ProteinIdentification> protein_identifications;
		vector<PeptideIdentification> identifications;

		if (pepxml_protxml)
		{
			mergePepXMLProtXML_(file_names, protein_identifications, identifications);
		}
		else
		{
			IdXMLFile().load(file_names[0], protein_identifications, identifications);

			vector<String> used_ids;
			for (Size i=1; i<file_names.size(); ++i)
			{
				vector<ProteinIdentification> additional_protein_identifications;
				vector<PeptideIdentification> additional_identifications;
				IdXMLFile().load(file_names[i], additional_protein_identifications, additional_identifications);
			
				for (Size i=0; i<additional_protein_identifications.size();++i)
				{
					if (find(used_ids.begin(), used_ids.end(), additional_protein_identifications[i].getIdentifier())!=used_ids.end())
					{
						writeLog_(String("Error: The identifier '") + additional_protein_identifications[i].getIdentifier() + "' was used before!");
						return INCOMPATIBLE_INPUT_DATA;
					}
					used_ids.push_back(additional_protein_identifications[i].getIdentifier());
				}
			
				protein_identifications.insert(protein_identifications.end(), additional_protein_identifications.begin(), additional_protein_identifications.end());
				identifications.insert(identifications.end(), additional_identifications.begin(), additional_identifications.end());
			}
		}
															
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
			
		IdXMLFile().store(out, protein_identifications, identifications);
			
		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDMerger tool;
	return tool.main(argc,argv);
}

/// @endcond
