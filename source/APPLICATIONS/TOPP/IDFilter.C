#// -*- mode: C++; tab-width: 2; -*-
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
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/APPLICATIONS/TOPPBase.h>

#include <limits>
#include <cmath>
#include <set>

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
//Doxygen docu
//-------------------------------------------------------------

/**
	@page TOPP_IDFilter IDFilter
	
	@brief Filters ProteinIdentification engine results by different criteria.
	
	This tool is used to filter the identifications found by
	a peptide/protein identification tool like Mascot&copy;. The identifications 
	can be filtered by different criteria.
	
	<ul>
		<li> 
			<b>peptide significance threshold</b>:<br> This parameter 
			specifies which amount of the significance threshold should 
			be reached by a peptide to be kept. If for example a peptide
			has score 30 and the significance threshold is 40, the 
			peptide will only be kept by the filter if the significance 
			threshold fraction is set to 0.75 or lower.
		</li> 
		<li> 
			<b>protein significance threshold</b>:<br> This parameter 
			behaves in the same way as the peptide significance threshold 
			fraction parameter. The only difference is that it is used
			to filter protein hits.
		</li>
		<li> 
			<b>peptide score</b>:<br> This parameter 
			specifies which score a peptide should have to be kept.
		</li> 
		<li> 
			<b>protein score</b>:<br> This parameter 
			specifies which score a protein should have to be kept.
		</li>
		<li>
			<b>peptide seqences</b>:<br> If you know which proteins
			are in the measured sample you can specify a FASTA file 
			which contains the protein sequences of those proteins. All 
			peptides which are not a substring of a protein contained
			in the sequences file will be filtered out. 
		</li>
		<li>
			<b>rt_filtering</b>:<br> To filter identifications according to their 
			predicted retention times you have to set this flag. You can set the used significance level
			by setting the 'p_value' parameter.<br>  
			This filter can only be applied to IdXML files produced by @ref TOPP_RTPredict.
		</li>
		<li>
			<b>exclusion peptides</b>:<br> For this option you specify an IdXML file.
			All peptides that are present in both files (in-file and exclusion peptides
			file) will be dropped.
		</li>
		<li>
			<b>best hits only</b>:<br> Only the best hit of a spectrum is kept.
			If there is more than one hit for a spectrum with the maximum score, then
			none of the hits will be kept.
		</li>
		<li>
			<b>best_n_peptide_hits</b>:<br> Only the best n peptide hits of a spectrum are kept.
		</li>
		<li>
			<b>best_n_protein_hits</b>:<br> Only the best n protein hits of a spectrum are kept.
		</li>
	</ul>

	<B>The command line parameters of this tool are:</B>
	@verbinclude TOPP_IDFilter.cli
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES


class TOPPIDFilter
	: public TOPPBase
{
 public:
	TOPPIDFilter()
		: TOPPBase("IDFilter","Filters results from protein or peptide identification engines based on different criteria.")
	{
		
	}
	
 protected:
	void registerOptionsAndFlags_()
	{
		registerInputFile_("in","<file>","","input file ");
		setValidFormats_("in",StringList::create("idXML"));
		registerOutputFile_("out","<file>","","output file ");
	  setValidFormats_("out",StringList::create("idXML"));
		registerInputFile_("sequences_file","<file>","","filename of a fasta file containing protein sequences.\n"
																											 "All peptides that are not a substring of a sequence in this file are filtered out",false);
		registerFlag_("no_protein_identifiers_in_seq_filter","If this flag is set the protein identifiers will not be used for sequence filtering");
		registerInputFile_("exclusion_peptides_file","<file>","","Peptides having the same sequence as any peptide in this file will be filtered out\n",false);
		setValidFormats_("exclusion_peptides_file",StringList::create("idXML"));
		registerDoubleOption_("pep_fraction","<fraction>",0.0,"the fraction of the peptide significance threshold that should be reached by a peptide hit",false);	
		registerDoubleOption_("prot_fraction","<fraction>",0.0,"the fraction of the protein significance threshold that should be reached by a protein hit",false);
		registerDoubleOption_("pep_score","<score>", 0,"the score which should be reached by a peptide hit to be kept",false);	
		registerDoubleOption_("prot_score","<score>", 0,"the score which should be reached by a protein hit to be kept",false);
		registerDoubleOption_("p_value","<significance>",0.05,"The probability of a correct ProteinIdentification having a deviation between observed and predicted rt equal or bigger than allowed",false);	
		registerIntOption_("best_n_peptide_hits","<score>", 0, "If this value is set only the n highest scoring peptide hits are kept per spectrum.", false);
		setMinInt_("best_n_peptide_hits", 1);
		registerIntOption_("best_n_protein_hits","<score>", 0, "If this value is set only the n highest scoring protein hits are kept.", false);
		setMinInt_("best_n_protein_hits", 1);
		registerIntOption_("min_length","<property>", 6, "If this value is set only peptide hits with a length greater or equal this value are kept.", false);
		setMinInt_("min_length", 1);
		registerFlag_("best_hits", "If this flag is set only the highest scoring hit is kept.\n"
															"If there are two or more highest scoring hits, none are kept.");
		registerFlag_("rt_filtering","If this flag is set rt filtering will be pursued.");
		registerFlag_("first_dim_rt","If this flag is set rt filtering will be pursued for first_dim.");
		registerFlag_("unique","If this flag is set and a peptide hit occurs more than once, only one instance is kept.");
		registerFlag_("unique_per_protein","If this flag is set, only peptides matching exactly one protein are kept.");
	}

	ExitCodes main_(int , const char**)
	{
	
		//-------------------------------------------------------------
		// varaibles
		//-------------------------------------------------------------
			
		IDFilter filter;
		IdXMLFile IdXML_file;
		vector<ProteinIdentification> protein_identifications;
		vector<PeptideIdentification> identifications;
		vector<PeptideIdentification> identifications_exclusion;
		vector<PeptideIdentification> filtered_peptide_identifications;
		vector<ProteinIdentification> filtered_protein_identifications;
		PeptideIdentification filtered_identification;
		ProteinIdentification filtered_protein_identification;
		vector<UInt> charges;
		vector< FASTAFile::FASTAEntry > sequences;
		set<String> exclusion_peptides;
		bool rt_filtering = false;
		bool first_dim_rt = false;
		DoubleReal p_value = 0.05;
		UInt min_length = 1;
		bool unique = false;
		
		
		//-------------------------------------------------------------
		// parsing parameters
		//-------------------------------------------------------------
			
		String inputfile_name = getStringOption_("in");			
		String outputfile_name = getStringOption_("out");
		
		DoubleReal peptide_significance_threshold_fraction = getDoubleOption_("pep_fraction");
		DoubleReal protein_significance_threshold_fraction = getDoubleOption_("prot_fraction");
		DoubleReal peptide_threshold_score = getDoubleOption_("pep_score");
		DoubleReal protein_threshold_score = getDoubleOption_("prot_score");
		
		Int best_n_peptide_hits = getIntOption_("best_n_peptide_hits");
		Int best_n_protein_hits = getIntOption_("best_n_protein_hits");
		min_length = getIntOption_("min_length");
		
		String sequences_file_name = getStringOption_("sequences_file");
		String exclusion_peptides_file_name = getStringOption_("exclusion_peptides_file");
		
		p_value = getDoubleOption_("p_value");
		rt_filtering = getFlag_("rt_filtering");
		first_dim_rt = getFlag_("first_dim_rt");
		unique = getFlag_("unique");

		bool strict = getFlag_("best_hits");
		bool no_protein_identifiers = getFlag_("no_protein_identifiers_in_seq_filter");
		bool unique_per_protein = getFlag_("unique_per_protein");
	
		//-------------------------------------------------------------
		// reading input
		//-------------------------------------------------------------
	

		if (sequences_file_name != "")
		{
			FASTAFile().load(sequences_file_name,sequences);				
		}
			
		if (exclusion_peptides_file_name  != "")
		{
			String document_id;
			IdXML_file.load(exclusion_peptides_file_name, protein_identifications, identifications_exclusion, document_id);
			for (Size i = 0; i < identifications_exclusion.size(); i++)
			{
				for(vector<PeptideHit>::const_iterator it = identifications_exclusion[i].getHits().begin();
						it != identifications_exclusion[i].getHits().end();
						it++)
				{
					exclusion_peptides.insert(it->getSequence().toString());
				}
			} 
		}												 
		String document_id;
		IdXML_file.load(inputfile_name, protein_identifications, identifications, document_id);

		//-------------------------------------------------------------
		// calculations
		//-------------------------------------------------------------
						
		// Filtering peptide identifications	according to set criteria
		for (Size i = 0; i < identifications.size(); i++)
		{
			if (unique_per_protein)
			{
				vector<PeptideHit> hits;
				for (vector<PeptideHit>::const_iterator it = identifications[i].getHits().begin(); it != identifications[i].getHits().end(); ++it)
				{
					if (!it->metaValueExists("protein_references"))
					{
						writeLog_("IDFilter: Warning, filtering with 'unique_per_protein' can only be done after indexing the file with 'PeptideIndexer' first.");
					}
					if (it->metaValueExists("protein_references") && (String)it->getMetaValue("protein_references") == "unique")
					{
						hits.push_back(*it);
					}
				}
				identifications[i].setHits(hits);
			}

			if (fabs(peptide_significance_threshold_fraction - 0) < 0.00001)
			{
				filtered_identification = identifications[i];
			}
			else
			{
				filter.filterIdentificationsByThreshold(identifications[i], peptide_significance_threshold_fraction, filtered_identification);
			}
			if (sequences_file_name != "")
			{
				PeptideIdentification temp_identification = filtered_identification;				
				filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_identification, no_protein_identifiers);
			}

			if (rt_filtering)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByRTPValues(temp_identification, filtered_identification, p_value);																																
			}

			if (first_dim_rt)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByRTFirstDimPValues(temp_identification, filtered_identification, p_value);																																
			}

			if (exclusion_peptides_file_name != "")
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByExclusionPeptides(temp_identification, exclusion_peptides, filtered_identification); 				
			}
			
			if (unique)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsUnique(temp_identification, filtered_identification); 				
			}
			
			if (strict)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByBestHits(temp_identification, filtered_identification, strict); 				
			}
			
			if (setByUser_("min_length"))
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByLength(temp_identification,
																						 min_length,
																						 filtered_identification); 								
			}

			if (setByUser_("pep_score"))
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByScore(temp_identification, peptide_threshold_score, filtered_identification); 				
			}

			if (best_n_peptide_hits != 0)
			{
				PeptideIdentification temp_identification = filtered_identification;
				filter.filterIdentificationsByBestNHits(temp_identification, best_n_peptide_hits, filtered_identification); 				
			}

			if(!filtered_identification.getHits().empty())
			{
			  PeptideIdentification tmp;
			  tmp = filtered_identification;
			  tmp.setMetaValue("RT", identifications[i].getMetaValue("RT"));
			  tmp.setMetaValue("MZ", identifications[i].getMetaValue("MZ"));
				filtered_peptide_identifications.push_back(tmp);
			}
		}
						
		// Filtering protein identifications	according to set criteria
		for (Size i = 0; i < protein_identifications.size(); i++)
		{
			if (!protein_identifications[i].getHits().empty())
			{
				if (!setByUser_("prot_fraction"))
				{       
					filtered_protein_identification = protein_identifications[i];
				}
				else
				{
					filter.filterIdentificationsByThreshold(protein_identifications[i], protein_significance_threshold_fraction, filtered_protein_identification);
				}
			
				if (sequences_file_name != "" && !no_protein_identifiers)
				{
					ProteinIdentification temp_identification = filtered_protein_identification;				
					filter.filterIdentificationsByProteins(temp_identification, sequences, filtered_protein_identification);
				}

				if (setByUser_("prot_score"))
				{
					ProteinIdentification temp_identification = filtered_protein_identification;
					filter.filterIdentificationsByScore(temp_identification, protein_threshold_score, filtered_protein_identification); 				
				}

				if (best_n_protein_hits > 0)
				{
					ProteinIdentification temp_identification = filtered_protein_identification;
					filter.filterIdentificationsByBestNHits(temp_identification, best_n_protein_hits, filtered_protein_identification); 				
				}
				
				ProteinIdentification temp_identification = filtered_protein_identification;
				filter.removeUnreferencedProteinHits(temp_identification, filtered_peptide_identifications, filtered_protein_identification); 								

				if(!(filtered_protein_identification.getHits().empty()))
				{
					filtered_protein_identifications.push_back(filtered_protein_identification);
				}
			}
			else
			{
				// copy the identifiers to the filtered protein ids
				filtered_protein_identifications.push_back(protein_identifications[i]);
			}
		}

		// check whether for each peptide identification identifier a corresponding protein id exists, if not add an empty one from the input file
		set<String> identifiers;
		for (vector<PeptideIdentification>::const_iterator it = filtered_peptide_identifications.begin(); it != filtered_peptide_identifications.end(); ++it)
		{
			identifiers.insert(it->getIdentifier());
		}

		for (set<String>::const_iterator it = identifiers.begin(); it != identifiers.end(); ++it)
		{
			// search for this identifier in filtered protein ids
			bool found(false);
			for (vector<ProteinIdentification>::const_iterator pit = filtered_protein_identifications.begin(); pit != filtered_protein_identifications.end(); ++pit)
			{
				if (*it == pit->getIdentifier())
				{
					found = true;
					break;
				}
			}

			if (!found)
			{
				// search this identifier in the protein id input
				found = false;
				ProteinIdentification new_prot_id;
				for (vector<ProteinIdentification>::const_iterator pit = protein_identifications.begin(); pit != protein_identifications.end(); ++pit)
				{
					if (*it == pit->getIdentifier())
					{
						new_prot_id = *pit;
						found = true;
						break;
					}
				}

				if (!found)
				{
					// this case means that the input file was not standard compatible
					writeLog_("Error: the identification run '" + *it + "' has no corresponding protein identification object!");
				}
				else
				{
					// just through away the protein hits
					new_prot_id.setHits(vector<ProteinHit>());
					filtered_protein_identifications.push_back(new_prot_id);
				}
			}
		}
		
		
		//-------------------------------------------------------------
		// writing output
		//-------------------------------------------------------------
	
		IdXML_file.store(outputfile_name, filtered_protein_identifications, filtered_peptide_identifications);

		return EXECUTION_OK;
	}
};


int main( int argc, const char** argv )
{
	TOPPIDFilter tool;

	return tool.main(argc,argv);
}

/// @endcond
