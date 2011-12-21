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
// $Maintainer: Nico Pfeifer $
// $Authors: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/ID/IDFilter.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <cmath>

using namespace std;

namespace OpenMS 
{
  IDFilter::IDFilter()
  {
  }
   
  IDFilter::~IDFilter()
  {
  }
  
  void IDFilter::filterIdentificationsUnique(const PeptideIdentification& identification,
  																					 PeptideIdentification& 			filtered_identification)
 	{
 		vector<PeptideHit> hits;
		filtered_identification = identification;		
 		vector<PeptideHit> temp_hits = identification.getHits();
 		
 		for(vector<PeptideHit>::iterator it = temp_hits.begin();
 				it != temp_hits.end();
 				++it)
 		{
 			if (find(hits.begin(), hits.end(), *it) == hits.end())
 			{
 				hits.push_back(*it);
 			}
 		}
 		filtered_identification.setHits(hits); 		
 	}  

	void IDFilter::filterIdentificationsByBestHits(const PeptideIdentification& identification,
																								 PeptideIdentification& filtered_identification, 
																								 bool strict)
	{
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		vector<Size> new_peptide_indices;		

		filtered_identification = identification;		
		filtered_identification.setHits(vector<PeptideHit>());
		
    if ( !identification.getHits().empty() )
		{
			Real optimal_value = identification.getHits()[0].getScore();
			new_peptide_indices.push_back(0);
			
			// searching for peptide(s) with maximal score			
			for (Size i = 1; i < identification.getHits().size(); i++)
			{
				Real temp_score = identification.getHits()[i].getScore();
        bool new_leader = false;
				if (   ( identification.isHigherScoreBetter() && (temp_score > optimal_value))
            || (!identification.isHigherScoreBetter() && (temp_score < optimal_value)) ) new_leader=true;
          
        if (new_leader)
				{
					optimal_value = temp_score;
					new_peptide_indices.clear();
					new_peptide_indices.push_back(i);
				}				
				else if (temp_score == optimal_value)
				{
					new_peptide_indices.push_back(i);
				}
			}						
			if (!strict || new_peptide_indices.size() == 1)
			{
				for (Size i = 0; i < new_peptide_indices.size(); i++)
				{
					filtered_peptide_hits.push_back(identification.getHits()[new_peptide_indices[i]]);
				}
			}
		}

    if ( !filtered_peptide_hits.empty() )
		{
  		filtered_identification.setHits(filtered_peptide_hits);
			filtered_identification.assignRanks();  																								
		}
	}

	void IDFilter::filterIdentificationsByLength(const PeptideIdentification& 	identification,
																							 Size            								min_length,
																							 PeptideIdentification& 				filtered_identification)
	{
		vector<Size> new_peptide_indices;		
		vector<PeptideHit> filtered_peptide_hits;
		
		filtered_identification = identification;		
		filtered_identification.setHits(vector<PeptideHit>());

		const vector<PeptideHit>& temp_peptide_hits = identification.getHits();

		for (Size i = 0; i < temp_peptide_hits.size(); i++)
		{
	  	if (temp_peptide_hits[i].getSequence().size() >= min_length)
	  	{
	  		new_peptide_indices.push_back(i);
			}				
		}		
		
		for (Size i = 0; i < new_peptide_indices.size(); i++)
		{
			filtered_peptide_hits.push_back(identification.getHits()[new_peptide_indices[i]]);
		}
    if ( !filtered_peptide_hits.empty() )
		{
  		filtered_identification.setHits(filtered_peptide_hits);
			filtered_identification.assignRanks();  																								
		}		
	}

	void IDFilter::filterIdentificationsByProteins(const PeptideIdentification& identification, 
																								 const vector< FASTAFile::FASTAEntry >& proteins,
																								 PeptideIdentification& filtered_identification,
																								 bool no_protein_identifiers)
	{
    // TODO: this is highly inefficient! the Protein-Index should be build once for all peptide-identifications instead of
    //       doing this once for every ID. Furthermore the index itself is inefficient (use seqan instead)
		String protein_sequences;
		String accession_sequences;
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification = identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for (Size i = 0; i < proteins.size(); i++)
		{
			if (proteins[i].identifier!="")
			{
				accession_sequences.append("*" + proteins[i].identifier);
			}
			if (proteins[i].sequence!="")
			{
				protein_sequences.append("*" + proteins[i].sequence);
			}
		}
		accession_sequences.append("*");
		protein_sequences.append("*");
		
		for (Size i = 0; i < identification.getHits().size(); i++)
		{
			if (no_protein_identifiers || accession_sequences=="*")
			{ // filter by sequence alone if no protein accesssions are available
		  	if (protein_sequences.find(identification.getHits()[i].getSequence().toUnmodifiedString()) != String::npos)
		  	{
		  		filtered_peptide_hits.push_back(identification.getHits()[i]);
		  	}
			}
			else
			{ // filter by protein accessions
				for(vector<String>::const_iterator ac_it = identification.getHits()[i].getProteinAccessions().begin();
						ac_it != identification.getHits()[i].getProteinAccessions().end();
						++ac_it)
				{
		  		if (accession_sequences.find("*" + *ac_it) != String::npos)
		  		{
		  			filtered_peptide_hits.push_back(identification.getHits()[i]);
            break; // we found a matching protein, the peptide is valid -> exit
		  		}
		  	}
			}
		}

    filtered_identification.setHits(filtered_peptide_hits);
		filtered_identification.assignRanks();																			
	}
	
	void IDFilter::filterIdentificationsByProteins(const ProteinIdentification& identification, 
																								 const vector< FASTAFile::FASTAEntry >& proteins,
                                                 ProteinIdentification& filtered_identification)
	{
		String protein_sequences;
		String accession_sequences;
		vector<ProteinHit> filtered_protein_hits;
		ProteinHit temp_protein_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<ProteinHit>());
		
		for (Size i = 0; i < proteins.size(); i++)
		{
			accession_sequences.append("*" + proteins[i].identifier);
		}
		accession_sequences.append("*");
		
		for (Size i = 0; i < identification.getHits().size(); i++)
		{
	  	if (accession_sequences.find("*" + identification.getHits()[i].getAccession()) != String::npos)
	  	{
	  		filtered_protein_hits.push_back(identification.getHits()[i]);
	  	}
		}

    filtered_identification.setHits(filtered_protein_hits);
		filtered_identification.assignRanks();											
	}
	
	void IDFilter::filterIdentificationsByExclusionPeptides(const PeptideIdentification& identification,
																													const set<String>& 				peptides,
																													PeptideIdentification& filtered_identification)
	{
		String protein_sequences;
		String accession_sequences;
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for (Size i = 0; i < identification.getHits().size(); i++)
		{
	  	if (find(peptides.begin(), peptides.end(), identification.getHits()[i].getSequence().toString()) == peptides.end())
	  	{
	  		filtered_peptide_hits.push_back(identification.getHits()[i]);
	  	}
		}
    if ( !filtered_peptide_hits.empty() )
		{
  		filtered_identification.setHits(filtered_peptide_hits);
			filtered_identification.assignRanks();
		}
	}
	
	void IDFilter::filterIdentificationsByRTFirstDimPValues(const PeptideIdentification& 	identification,
																										 			PeptideIdentification& 				filtered_identification,
																										 			DoubleReal 										p_value)
	{
		DoubleReal border = 1 - p_value;
		vector< Size > new_peptide_indices;		
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
    Size missing_meta_value=0;

    for ( Size i = 0; i < identification.getHits().size(); ++i )
    {
			if (identification.getHits()[i].metaValueExists("predicted_RT_p_value_first_dim"))
      {
			  if ((DoubleReal)(identification.getHits()[i].getMetaValue("predicted_RT_p_value_first_dim")) <= border )
			  {
		  	  filtered_peptide_hits.push_back(identification.getHits()[i]);
			  }
      }
      else ++missing_meta_value;
		}
    if (missing_meta_value>0) LOG_WARN << "Filtering identifications by p-value did not work on " << missing_meta_value << " of " << identification.getHits().size() << " hits. Your data is missing a meta-value ('predicted_RT_p_value_first_dim') from RTPredict!\n";

    if ( !filtered_peptide_hits.empty() )
		{
  		filtered_identification.setHits(filtered_peptide_hits);		
  		filtered_identification.assignRanks();											
		}
	}																		 																										 

	void IDFilter::filterIdentificationsByRTPValues(const PeptideIdentification& 	identification,
																						 			PeptideIdentification& 				filtered_identification,
																						 			DoubleReal 										p_value)
	{
		DoubleReal border = 1 - p_value;
		vector< Size > new_peptide_indices;		
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
    Size missing_meta_value=0;

		for (Size i = 0; i < identification.getHits().size(); i++)
		{
			if (identification.getHits()[i].metaValueExists("predicted_RT_p_value"))
      {
			  if ((DoubleReal)(identification.getHits()[i].getMetaValue("predicted_RT_p_value")) <= border )
			  {
		  	  filtered_peptide_hits.push_back(identification.getHits()[i]);
			  }		
      }
      else ++missing_meta_value;
		}
    if (missing_meta_value>0) LOG_WARN << "Filtering identifications by p-value did not work on " << missing_meta_value << " of " << identification.getHits().size() << " hits. Your data is missing a meta-value ('predicted_RT_p_value') from RTPredict!\n";

    if ( !filtered_peptide_hits.empty() )
		{
  		filtered_identification.setHits(filtered_peptide_hits);		
  		filtered_identification.assignRanks();											
		}
	}
	
	void IDFilter::removeUnreferencedProteinHits(const ProteinIdentification& 	identification, 
																							const vector<PeptideIdentification> peptide_identifications, 
																							ProteinIdentification& 	filtered_identification)
	{
		vector<ProteinHit> filtered_protein_hits;
		const vector<ProteinHit>& temp_protein_hits = identification.getHits();
		vector<PeptideHit> temp_peptide_hits;

		filtered_identification=identification;
		filtered_identification.setHits(vector<ProteinHit>());
		String identifier = identification.getIdentifier();
		
		Size i = 0;		
		for (Size j = 0; j < temp_protein_hits.size(); ++j)
		{
			bool found = false;
			i = 0;
 			while(i < peptide_identifications.size() && !found) 
			{
				if (identifier == peptide_identifications[i].getIdentifier())
				{
					temp_peptide_hits.clear();
					peptide_identifications[i].getReferencingHits(temp_protein_hits[j].getAccession(), temp_peptide_hits);
          if ( !temp_peptide_hits.empty() )
					{
						filtered_protein_hits.push_back(temp_protein_hits[j]);
						found = true;
					}
				}
				++i;
			}
		}
		filtered_identification.setHits(filtered_protein_hits);
	}																		 																										 
		
} // namespace OpenMS
