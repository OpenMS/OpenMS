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
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/Map.h>

namespace OpenMS
{
	
	ProteinInference::ProteinInference()
	{
	}

	ProteinInference::ProteinInference(const ProteinInference& /*cp*/)
	{
	}

	ProteinInference& ProteinInference::operator = (const ProteinInference& /*rhs*/)
	{
		return *this;
	}
	
	void ProteinInference::infer(ConsensusMap& consensus_map, const UInt reference_map)
	{
		// we infer Proteins for every IdentificationRun seperately. If you want this combined, then
		// do that before calling this function
		// Each ProteinIdentification will be augmented with the quantification (where possible)
		for (size_t i = 0;
				 i < consensus_map.getProteinIdentifications().size();
				 ++i)
		{
			infer_(consensus_map, i, reference_map);
		}
	}
	
	void ProteinInference::infer_(ConsensusMap& consensus_map, 
																const size_t protein_idenfication_index, 
																const UInt reference_map)
	{
		
		ProteinIdentification& protein_ident = consensus_map.getProteinIdentifications()[protein_idenfication_index];
		for (size_t i=0; i< protein_ident.getHits().size(); ++i)
		{
			// Protein Accession
			String accession = protein_ident.getHits()[i].getAccession();

			// consensusfeature -> peptideHit
			Map < size_t, PeptideHit> consensus_to_peptide;
			
			// search for it in consensusElements:
			for (size_t i_cm=0; i_cm < consensus_map.size(); ++i_cm)
			{
				std::vector< PeptideHit > peptide_hits;
				for (std::vector< PeptideIdentification >::iterator it_pepid = consensus_map[i_cm].getPeptideIdentifications().begin();
						 it_pepid != consensus_map[i_cm].getPeptideIdentifications().end();
						 ++it_pepid)
				{
					// are Protein- and PeptideIdentification from the same search engine run?
					if (it_pepid->getIdentifier() != protein_ident.getIdentifier()) continue;
					
					std::vector< PeptideHit > peptide_hits_local;
					
					it_pepid->getReferencingHits(accession, peptide_hits_local);
					
					if (peptide_hits_local.empty()) continue;
					
					if (sortByUnique_(peptide_hits_local, it_pepid->isHigherScoreBetter() ))
					{ // we found a unique peptide
						peptide_hits.push_back(peptide_hits_local[0]);
					}
					
				}
				
				// if several PeptideIdentifications (==Spectra) were assigned to current ConsensusElement
				// --> take the best (as above), e.g. in SILAC this could happen
				// TODO better idea?
				if (peptide_hits.size()>0)
				{
					if (sortByUnique_(peptide_hits, consensus_map[i_cm].getPeptideIdentifications()[0].isHigherScoreBetter() ))
					{ //found a unique peptide for current ConsensusElement
						consensus_to_peptide[i_cm] = peptide_hits[0];
						#ifdef DEBUG_INFERENCE
					std::cout << "assign peptide " <<  peptide_hits[0].getSequence() << " to Protein " << accession << std::endl;
						#endif
					}
				}
				
			} // ! ConsensusMap loop
			
			// no peptides found that match current Protein
			if (consensus_to_peptide.size() == 0) continue;
			
			// Use all matching ConsensusElements to derive a quantitation for current protein
			// build up ratios for every map vs reference
			double coverage = 0;
			Map < Size, std::vector < IntensityType > > ratios;
			
			// number of unique peptides pointing to current protein
			UInt coverage_count = (UInt)consensus_to_peptide.size();
			
			for (Map < size_t, PeptideHit>::iterator it_pephits = consensus_to_peptide.begin();
					 it_pephits != consensus_to_peptide.end();
					 ++it_pephits)
			{
				coverage += it_pephits->second.getSequence().size();
				const ConsensusFeature::HandleSetType& handles = consensus_map[it_pephits->first].getFeatures();
				//search if reference is present
				ConsensusFeature::HandleSetType::const_iterator it_ref = handles.end();
				for (ConsensusFeature::HandleSetType::const_iterator it=handles.begin();
						 it != handles.end();
						 ++it)
				{
					if (it->getMapIndex() == reference_map)
					{
						it_ref = it;
						break;
					}
				}
				
				// did not find a reference
				// TODO assume intensity==0 instead??
				if (it_ref == handles.end()) continue; 
				
				for (ConsensusFeature::HandleSetType::const_iterator it=handles.begin();
						 it != handles.end();
						 ++it)
				{
					ratios[it->getMapIndex()].push_back(it->getIntensity() / it_ref->getIntensity());
				}
				
			}
			
			// sort ratios map-wise and take median
			for (ConsensusMap::FileDescriptions::const_iterator it_file = consensus_map.getFileDescriptions().begin();
					 it_file != consensus_map.getFileDescriptions().end();
					 ++it_file)
			{
				if (ratios.has(it_file->first))
				{
					//sort intensity ratios for map #it_file->first
					std::sort(ratios[it_file->first].begin(), ratios[it_file->first].end());
					//take median
					IntensityType protein_ratio = ratios[it_file->first][ratios[it_file->first].size()/2];
					
					//TODO if ratios have high variance emit a warning!
					
					protein_ident.getHits()[i].setMetaValue(String("ratio_") + String(it_file->first), protein_ratio);
				}
				
			} // ! map loop
			
			// % coverage of protein by peptides
			coverage /= DoubleReal(protein_ident.getHits()[i].getSequence().size()) / 100;
			
			protein_ident.getHits()[i].setMetaValue("coverage", coverage);
			protein_ident.getHits()[i].setMetaValue("hits", coverage_count);
			
		} // ! Protein loop
		
		
		
		// protein_to_peptides now contains the Protein -> Peptides mapping
		// lets estimate the
		
	}
	

	bool ProteinInference::sortByUnique_(std::vector< PeptideHit >& peptide_hits_local, const bool is_higher_score_better )
	{
		if (peptide_hits_local.empty()) return false;

		// several peptideHits from (the same) spectrum point to current Protein
		// -> take the best
		if (peptide_hits_local.size() > 1)
		{
			std::sort(peptide_hits_local.begin(), peptide_hits_local.end(), PeptideHit::ScoreLess() );
			if (is_higher_score_better)
			{
				peptide_hits_local[0] = peptide_hits_local[peptide_hits_local.size()-1];
			}
		}
		
		//-> lets see if its unique:
		if (peptide_hits_local[0].getProteinAccessions().size()!=1)
		{
			// this is a shared peptide --> do not use it
			return false;
		}
		else
		{
			return true;
		}
		// the first element now contains the best peptideHit

	}

}
 
