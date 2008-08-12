// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/ID/IDFilter.h>

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
  

	void IDFilter::filterIdentificationsByBestHits(const PeptideIdentification& 	identification,
																								 PeptideIdentification& 				filtered_identification, 
																								 bool 									strict)
	{
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		vector< UInt > new_peptide_indices;		

		filtered_identification = identification;		
		filtered_identification.setHits(vector<PeptideHit>());
		
		if (identification.getHits().size() > 0)
		{
			Real max_value = identification.getHits()[0].getScore();
			new_peptide_indices.push_back(0);
			
			// searching for peptide(s) with maximal score			
			for(UInt i = 1; i < identification.getHits().size(); i++)
			{
				Real temp_score = identification.getHits()[i].getScore();
				if (temp_score > max_value)
				{
					max_value = temp_score;
					new_peptide_indices.clear();
					new_peptide_indices.push_back(i);
				}				
				else if (temp_score == max_value)
				{
					new_peptide_indices.push_back(i);
				}
			}						
			if (!strict || new_peptide_indices.size() == 1)
			{
				for(UInt i = 0; i < new_peptide_indices.size(); i++)
				{
					filtered_peptide_hits.push_back(identification.getHits()[new_peptide_indices[i]]);
				}
			}
		}

		if (filtered_peptide_hits.size() > 0)
		{
  		filtered_identification.setHits(filtered_peptide_hits);
			filtered_identification.assignRanks();  																								
		}
	}

	void IDFilter::filterIdentificationsByProteins(const PeptideIdentification& identification, 
																								 const vector< FASTAFile::FASTAEntry >& proteins,
																								 PeptideIdentification& filtered_identification)
	{
		String protein_sequences;
		String accession_sequences;
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification = identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for(UInt i = 0; i < proteins.size(); i++)
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
		
		for(UInt i = 0; i < identification.getHits().size(); i++)
		{
			if (accession_sequences=="*")
			{
		  	if (protein_sequences.find(identification.getHits()[i].getSequence().toUnmodifiedString()) != string::npos)
		  	{
		  		filtered_peptide_hits.push_back(identification.getHits()[i]);
		  	}
			}
			else
			{
				for(vector<String>::const_iterator ac_it = identification.getHits()[i].getProteinAccessions().begin();
						ac_it != identification.getHits()[i].getProteinAccessions().end();
						++ac_it)
				{
		  		if (accession_sequences.find("*" + *ac_it) != string::npos)
		  		{
		  			filtered_peptide_hits.push_back(identification.getHits()[i]);
		  		}
		  	}
			}
		}
		if (filtered_peptide_hits.size() > 0)
		{
  		filtered_identification.setHits(filtered_peptide_hits);
  		filtered_identification.assignRanks();																			
		}
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
		
		for(UInt i = 0; i < proteins.size(); i++)
		{
			accession_sequences.append("*" + proteins[i].identifier);
		}
		accession_sequences.append("*");
		
		for(UInt i = 0; i < identification.getHits().size(); i++)
		{
	  	if (accession_sequences.find("*" + identification.getHits()[i].getAccession()) != string::npos)
	  	{
	  		filtered_protein_hits.push_back(identification.getHits()[i]);
	  	}
		}
		if (filtered_protein_hits.size() > 0)
		{
      filtered_identification.setHits(filtered_protein_hits);
			filtered_identification.assignRanks();											
		}
	}
	
	void IDFilter::filterIdentificationsByExclusionPeptides(const PeptideIdentification& identification,
																													vector<String> 				peptides,
																													PeptideIdentification& filtered_identification)
	{
		String protein_sequences;
		String accession_sequences;
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for(UInt i = 0; i < identification.getHits().size(); i++)
		{
	  	if (find(peptides.begin(), peptides.end(), identification.getHits()[i].getSequence().toUnmodifiedString()) == peptides.end())
	  	{
	  		filtered_peptide_hits.push_back(identification.getHits()[i]);
	  	}
		}
		if (filtered_peptide_hits.size() > 0)
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
		vector< UInt > new_peptide_indices;		
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for(UInt i = 0; i < identification.getHits().size(); i++)
		{
			if (identification.getHits()[i].metaValueExists("predicted_RT_p_value_first_dim") 
			    && (DoubleReal)(identification.getHits()[i].getMetaValue("predicted_RT_p_value_first_dim")) <= border )
			{
		  	filtered_peptide_hits.push_back(identification.getHits()[i]);
			}		
		}
		if (filtered_peptide_hits.size() > 0)
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
		vector< UInt > new_peptide_indices;		
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		
		filtered_identification=identification;
		filtered_identification.setHits(vector<PeptideHit>());
		
		for(UInt i = 0; i < identification.getHits().size(); i++)
		{
			if (identification.getHits()[i].metaValueExists("predicted_RT_p_value") 
			    && (DoubleReal)(identification.getHits()[i].getMetaValue("predicted_RT_p_value")) <= border )
			{
		  	filtered_peptide_hits.push_back(identification.getHits()[i]);
			}		
		}
		if (filtered_peptide_hits.size() > 0)
		{
  		filtered_identification.setHits(filtered_peptide_hits);		
  		filtered_identification.assignRanks();											
		}
	}																		 																										 
		
} // namespace OpenMS
