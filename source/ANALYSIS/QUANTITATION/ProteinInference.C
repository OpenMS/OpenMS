// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
        // TODO: better idea?
        if ( !peptide_hits.empty() )
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
			if (consensus_to_peptide.empty()) continue;
			
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
 
