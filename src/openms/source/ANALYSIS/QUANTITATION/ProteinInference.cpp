// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/ANALYSIS/QUANTITATION/ProteinInference.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <OpenMS/KERNEL/ConsensusMap.h>

#include <map>

namespace OpenMS
{

  ProteinInference::ProteinInference() = default;

  ProteinInference::ProteinInference(const ProteinInference& /*cp*/) = default;

  ProteinInference& ProteinInference::operator=(const ProteinInference& /*rhs*/) = default;

  void ProteinInference::infer(ConsensusMap& consensus_map, const UInt reference_map)
  {
    // we infer Proteins for every IdentificationRun separately. If you want this combined, then
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
    for (size_t i = 0; i < protein_ident.getHits().size(); ++i)
    {
      // Protein Accession
      String accession = protein_ident.getHits()[i].getAccession();

      // consensus feature -> peptide hit
      std::map<size_t, PeptideHit> consensus_to_peptide;

      // search for it in consensus elements:
      for (size_t i_cm = 0; i_cm < consensus_map.size(); ++i_cm)
      {
        std::vector<PeptideHit> peptide_hits;
        for (std::vector<PeptideIdentification>::iterator it_pepid = consensus_map[i_cm].getPeptideIdentifications().begin();
             it_pepid != consensus_map[i_cm].getPeptideIdentifications().end();
             ++it_pepid)
        {
          // are Protein- and PeptideIdentification from the same search engine run?
          if (it_pepid->getIdentifier() != protein_ident.getIdentifier())
            continue;

          std::set<String> accessions;
          accessions.insert(accession);
          std::vector<PeptideHit> peptide_hits_local = PeptideIdentification::getReferencingHits(it_pepid->getHits(), accessions);

          if (peptide_hits_local.empty())
          {
            continue;
          }

          if (sortByUnique_(peptide_hits_local, it_pepid->isHigherScoreBetter())) // we found a unique peptide
          {
            peptide_hits.push_back(peptide_hits_local[0]);
          }

        }

        // if several PeptideIdentifications (==Spectra) were assigned to current ConsensusElement
        // --> take the best (as above), e.g. in SILAC this could happen
        // TODO: better idea?
        if (!peptide_hits.empty())
        {
          if (sortByUnique_(peptide_hits, consensus_map[i_cm].getPeptideIdentifications()[0].isHigherScoreBetter())) //found a unique peptide for current ConsensusElement
          {
            consensus_to_peptide[i_cm] = peptide_hits[0];
#ifdef DEBUG_INFERENCE
            std::cout << "assign peptide " <<  peptide_hits[0].getSequence() << " to Protein " << accession << std::endl;
#endif
          }
        }

      } // ! ConsensusMap loop

      // no peptides found that match current Protein
      if (consensus_to_peptide.empty())
        continue;

      // Use all matching ConsensusElements to derive a quantitation for current protein
      // build up ratios for every map vs reference
      double coverage = 0;
      std::map<Size, std::vector<IntensityType> > ratios;

      // number of unique peptides pointing to current protein
      UInt coverage_count = (UInt)consensus_to_peptide.size();

      for (std::map<size_t, PeptideHit>::iterator it_pephits = consensus_to_peptide.begin();
           it_pephits != consensus_to_peptide.end();
           ++it_pephits)
      {
        coverage += it_pephits->second.getSequence().size();
        const ConsensusFeature::HandleSetType& handles = consensus_map[it_pephits->first].getFeatures();
        //search if reference is present
        ConsensusFeature::HandleSetType::const_iterator it_ref = handles.end();
        for (ConsensusFeature::HandleSetType::const_iterator it = handles.begin();
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
        if (it_ref == handles.end())
          continue;

        for (ConsensusFeature::HandleSetType::const_iterator it = handles.begin();
             it != handles.end();
             ++it)
        {
          ratios[it->getMapIndex()].push_back(it->getIntensity() / it_ref->getIntensity());
        }

      }

      // sort ratios map-wise and take median
      for (ConsensusMap::ColumnHeaders::const_iterator it_file = consensus_map.getColumnHeaders().begin();
           it_file != consensus_map.getColumnHeaders().end();
           ++it_file)
      {
        if (ratios.find(it_file->first) != ratios.end())
        {
          //sort intensity ratios for map #it_file->first
          std::sort(ratios[it_file->first].begin(), ratios[it_file->first].end());
          //take median
          IntensityType protein_ratio = ratios[it_file->first][ratios[it_file->first].size() / 2];

          //TODO if ratios have high variance emit a warning!

          protein_ident.getHits()[i].setMetaValue(String("ratio_") + String(it_file->first), protein_ratio);
        }

      } // ! map loop

      // % coverage of protein by peptides
      coverage /= double(protein_ident.getHits()[i].getSequence().size()) / 100;

      protein_ident.getHits()[i].setMetaValue("coverage", coverage);
      protein_ident.getHits()[i].setMetaValue("hits", coverage_count);

    } // ! Protein loop



    // protein_to_peptides now contains the Protein -> Peptides mapping
    // lets estimate the

  }

  bool ProteinInference::sortByUnique_(std::vector<PeptideHit>& peptide_hits_local, const bool is_higher_score_better)
  {
    if (peptide_hits_local.empty())
      return false;

    // several peptideHits from (the same) spectrum point to current Protein
    // -> take the best
    if (peptide_hits_local.size() > 1)
    {
      std::sort(peptide_hits_local.begin(), peptide_hits_local.end(), PeptideHit::ScoreLess());
      if (is_higher_score_better)
      {
        peptide_hits_local[0] = peptide_hits_local[peptide_hits_local.size() - 1];
      }
    }

    //-> lets see if its unique:
    std::set<String> protein_accessions = peptide_hits_local[0].extractProteinAccessionsSet();
    if (protein_accessions.size() != 1)
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
