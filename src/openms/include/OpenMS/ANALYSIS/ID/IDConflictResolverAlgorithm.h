// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser, Lucia Espona, Moritz Freidank $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------



namespace OpenMS
{

/**
    @brief Resolves ambiguous annotations of features with peptide identifications.

    The peptide identifications are filtered so that only one identification
    with a single hit (with the best score) is associated to each feature.
    (If two IDs have the same best score, either one of them may be selected.)
*/
class OPENMS_DLLAPI IDConflictResolverAlgorithm
{
public:
  /** @brief Resolves ambiguous annotations of features with peptide identifications.
    The the filtered identifications are added to the vector of unassigned peptides
    and also reduced to a single best hit.

    @param features Features to work on
    @param keep_matching Keeps all IDs that match the modified sequence of the best
    hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
  **/
  static void resolve(FeatureMap& features, bool keep_matching = false);

  /** @brief Resolves ambiguous annotations of consensus features with peptide identifications.
    The the filtered identifications are added to the vector of unassigned peptides
    and also reduced to a single best hit.
    
    @param features Features to work on
    @param keep_matching Keeps all IDs that match the modified sequence of the best
    hit in the feature (e.g. keeps all IDs in a ConsensusMap if id'd same across multiple runs)
  **/
  static void resolve(ConsensusMap& features, bool keep_matching = false);

  /** @brief In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear.
   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
   appears only once, i.e. no multiplicities.
   **/
  static void resolveBetweenFeatures(FeatureMap& features);
  
  /** @brief In a single (feature/consensus) map, features with the same (possibly modified) sequence and charge state may appear.
   This filter removes the peptide sequence annotations from features, if a higher-intensity feature with the same (charge, sequence)
   combination exists in the map. The total number of features remains unchanged. In the final output, each (charge, sequence) combination
   appears only once, i.e. no multiplicities.
   **/
  static void resolveBetweenFeatures(ConsensusMap& features);
  
protected:

  template<class T>
  static void resolveConflict_(T& map, bool keep_matching)
  {
    // annotate as not part of the resolution
    for (PeptideIdentification& p : map.getUnassignedPeptideIdentifications())
    {
      p.setMetaValue("feature_id", "not mapped"); // not mapped to a feature
    }

    for (auto& c : map)
    {
      c.setMetaValue("feature_id", String(c.getUniqueId()));
      if (!keep_matching)
      {
        resolveConflict_(c.getPeptideIdentifications(),
                         map.getUnassignedPeptideIdentifications(),
                         c.getUniqueId());
      }
      else
      {
        resolveConflictKeepMatching_(c.getPeptideIdentifications(),
                         map.getUnassignedPeptideIdentifications(),
                         c.getUniqueId());
      }
    }
  }
  
  // compare peptide IDs by score of best hit (hits must be sorted first!)
  // (note to self: the "static" is necessary to avoid cryptic "no matching
  // function" errors from gcc when the comparator is used below)
  static bool compareIDsSmallerScores_(const PeptideIdentification & left,
                          const PeptideIdentification & right);

  static void resolveConflict_(
    std::vector<PeptideIdentification> & peptides,
    std::vector<PeptideIdentification> & removed,
    UInt64 uid);

  static void resolveConflictKeepMatching_(
      std::vector<PeptideIdentification> & peptides,
      std::vector<PeptideIdentification> & removed,
      UInt64 uid);
  
  template<class T>
  static void resolveBetweenFeatures_(T & map)
  {
    // unassigned peptide identifications in this map
    std::vector<PeptideIdentification>& unassigned = map.getUnassignedPeptideIdentifications();
    
    // A std::map tracking the set of unique features.
    // Uniqueness criterion/key is a pair <charge, sequence> for each feature. The peptide sequence may be modified, i.e. is not stripped.
    typedef std::map<std::pair<Int, AASequence>, typename T::value_type*> FeatureSet;
    FeatureSet feature_set;
    
    // Create a std::map `feature_set` mapping pairs <charge, sequence> to a pointer to
    // the feature with the highest intensity for this sequence.
    for (typename T::value_type& element : map)
    {
      std::vector<PeptideIdentification>& pep_ids = element.getPeptideIdentifications();
      
      if (!pep_ids.empty())
      {
        if (pep_ids.size() != 1)
        {
          // Should never happen. In IDConflictResolverAlgorithm TOPP tool
          // IDConflictResolverAlgorithm::resolve() is called before IDConflictResolverAlgorithm::resolveBetweenFeatures().
          throw OpenMS::Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__, "Feature does contain multiple identifications.");
        }
        
        // Make sure best hit is in front, i.e. sort hits first.
        pep_ids.front().sort();
        const std::vector<PeptideHit>& hits = pep_ids.front().getHits();
        
        if (!hits.empty())
        {
          const PeptideHit& highest_score_hit = hits.front();
          
          // Pair <charge, sequence> of charge of the new feature and the sequence of its highest scoring peptide hit.
          std::pair<Int, AASequence> pair = std::make_pair(element.getCharge(), highest_score_hit.getSequence());
          
          // If a <charge, sequence> pair is not yet in the FeatureSet or new feature `feature_in_set`
          // has higher intensity than its counterpart `feature_set[<charge, sequence>]`
          // store a pointer to `feature_in_set` in `feature_set`.
          typename FeatureSet::iterator feature_in_set = feature_set.find(pair);
          if (feature_in_set != feature_set.end())
          {
            // Identical (charge, sequence) key found. Remove annotations from either the old or new feature.
            
            if (feature_in_set->second->getIntensity() < element.getIntensity())
            {
              // Remove annotations from the old low-intensity feature. But only after moving these annotations to the unassigned list.
              std::vector<PeptideIdentification>& obsolete = feature_in_set->second->getPeptideIdentifications();
              unassigned.insert(unassigned.end(), obsolete.begin(), obsolete.end());
              std::vector<PeptideIdentification> pep_ids_empty;
              feature_in_set->second->setPeptideIdentifications(pep_ids_empty);
              
              // Replace feature in the set.
              feature_in_set->second = &(element);
            }
            else
            {
              // Remove annotations from the new low-intensity feature. But only after moving these annotations to the unassigned list.
              std::vector<PeptideIdentification>& obsolete = element.getPeptideIdentifications();
              unassigned.insert(unassigned.end(), obsolete.begin(), obsolete.end());
              std::vector<PeptideIdentification> pep_ids_empty;
              element.setPeptideIdentifications(pep_ids_empty);
            }
          }
          else
          {
            // Feature is not yet in our set -- add it.
            feature_set[pair] = &(element);
          }
        }
      }
    }
  }
  
};

}// namespace OpenMS

